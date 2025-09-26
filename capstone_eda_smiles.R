# Read your tab-delimited TXT file
df <- read.delim("C:/Users/Isaac/Downloads/Repurposing_Hub_exportidkbro.txt", sep = "\t", header = TRUE)

# Write it out as CSV
df_new<-read.csv("C:/Users/Isaac/Downloads/Repurposing_Hub_export3.csv")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)

df_neuro <- df_new %>%
  filter(Disease.Area == "neurology/psychiatry" | Disease.Area == "")%>%
  mutate(rowid = row_number())
summary(df_neuro)

split_dedupe_first <- function(x) {
  xs <- str_split(as.character(x %||% ""), "\\s*,\\s*")[[1]]
  xs <- xs[nzchar(xs)]
  xs <- unique(xs)
  if (length(xs) == 0) NA_character_ else xs[1]
}

dupe_count <- function(x) {
  xs <- str_split(as.character(x %||% ""), "\\s*,\\s*")[[1]]
  length(unique(xs[nzchar(xs)]))
}

df_neuro <- df_neuro %>%
  mutate(SMILES = map_chr(SMILES, split_dedupe_first),
         smiles_uniques_in_cell = map_int(SMILES, dupe_count))
table(df_neuro$smiles_uniques_in_cell)  # sanity check

# --- class balance within neuro: Indication
ind_counts_neuro <- df_neuro %>%
  count(Indication, sort = TRUE) %>%
  mutate(pct = n / sum(n))
print(ind_counts_neuro %>% slice_head(n = 30))  # top 30

ggplot(ind_counts_neuro %>% slice_head(n = 50),
       aes(x = fct_reorder(Indication, n), y = n)) +
  geom_col() + coord_flip() +
  labs(title = "Neurology/Psychiatry: Top 25 Indications",
       x = "Indication", y = "Count")

# --- MOA / Target frequencies (neuro only)
moa_counts_neuro <- df_neuro %>% count(MOA, sort = TRUE)
target_counts_neuro <- df_neuro %>% count(Target, sort = TRUE)
print(moa_counts_neuro %>% slice_head(n = 20))
print(target_counts_neuro %>% slice_head(n = 20))

ggplot(moa_counts_neuro %>% slice_head(n = 30),
       aes(x = fct_reorder(MOA, n), y = n)) +
  geom_col() + coord_flip() +
  labs(title = "Neurology/Psychiatry: Top 25 MOA", x = "MOA", y = "Count")

ggplot(target_counts_neuro %>% slice_head(n = 25),
       aes(x = fct_reorder(Target, n), y = n)) +
  geom_col() + coord_flip() +
  labs(title = "Neurology/Psychiatry: Top 25 Targets", x = "Target", y = "Count")

# --- Phase distribution (neuro only)
df_neuro_phase <- df_neuro %>% mutate(Phase = as.factor(Phase))
phase_counts_neuro <- df_neuro_phase %>% count(Phase, sort = TRUE)
print(phase_counts_neuro)

ggplot(phase_counts_neuro, aes(x = Phase, y = n)) +
  geom_col() +
  labs(title = "Neurology/Psychiatry: Phase Distribution", x = "Phase", y = "Count")

# --- co-occurrence heatmaps: MOA × Indication and Target × Indication (limited)
moa_top_neuro <- moa_counts_neuro %>% slice_head(n = 40) %>% pull(MOA)
ind_top_neuro <- ind_counts_neuro %>% slice_head(n = 40) %>% pull(Indication)

df_neuro %>%
  filter(MOA %in% moa_top_neuro, Indication %in% ind_top_neuro) %>%
  count(MOA, Indication) %>%
  ggplot(aes(x = Indication, y = MOA, fill = n)) +
  geom_tile() + coord_fixed() +
  labs(title = "Neurology/Psychiatry: MOA × Indication (Top 30×30)",
       x = "Indication", y = "MOA", fill = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

target_top_neuro <- target_counts_neuro %>% slice_head(n = 20) %>% pull(Target)

df_neuro %>%
  filter(Target %in% target_top_neuro, Indication %in% ind_top_neuro) %>%
  count(Target, Indication) %>%
  ggplot(aes(x = Indication, y = Target, fill = n)) +
  geom_tile() + coord_fixed() +
  labs(title = "Neurology/Psychiatry: Target × Indication (Top 20×20)",
       x = "Indication", y = "Target", fill = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/Isaac/AppData/Local/r-miniconda/envs/r-rdkit/python.exe")

rd     <- import("rdkit")
Chem   <- import("rdkit.Chem")
Descs  <- import("rdkit.Chem.Descriptors")                 # <-- add this
rdmops <- import("rdkit.Chem.rdmolops")
rdms   <- import("rdkit.Chem.MolStandardize.rdMolStandardize")
rfg    <- import("rdkit.Chem.rdFingerprintGenerator")
DS     <- import("rdkit.DataStructs")                      # <-- add this
np     <- import("numpy")                                   # <-- add this

rd$RDLogger$DisableLog("rdApp.error"); rd$RDLogger$DisableLog("rdApp.warning")

# clean + strip CXSMILES
clean_smiles <- function(x) {
  x <- gsub("\\|.*\\|\\s*$", "", x, perl = TRUE)            # <-- strip CXSMILES payload
  x <- gsub("[\uFEFF\u200B\u00A0]", "", x)
  x <- gsub("[\u2013\u2014]", "-", x)
  x <- gsub("[\"’]", "", x)
  x <- trimws(x)
  x[nchar(x) == 0 | tolower(x) %in% c("na","n/a","nan","null")] <- NA
  x
}
# allow common SMILES chars only; reject anything else
looks_like_smiles <- function(s) {
  if (is.na(s) || !nzchar(s)) return(FALSE)
  # Put '-' at the END of the class (or escape it) and use perl=TRUE
  grepl("^[A-Za-z0-9@+#=()/\\\\\\[\\].%*.:\\-]+$", s, perl = TRUE)
  #                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ allowed chars only
}

safe_mol <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  s <- clean_smiles(s)
  if (!looks_like_smiles(s)) return(NULL)
  
  m <- tryCatch(Chem$MolFromSmiles(s), error = function(e) NULL)
  if (is.null(m)) {
    m <- tryCatch(Chem$MolFromSmiles(s, sanitize = FALSE), error = function(e) NULL)
    if (is.null(m)) return(NULL)
  }
  
  m <- tryCatch(rdms$Cleanup(m), error = function(e) m)
  m <- tryCatch(rdms$LargestFragmentChooser()$choose(m), error = function(e) m)
  m <- tryCatch(rdms$Uncharger()$uncharge(m), error = function(e) m)
  m <- tryCatch(rdms$TautomerEnumerator()$Canonicalize(m), error = function(e) m)
  
  m <- tryCatch({ rdmops$SanitizeMol(m); m }, error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  smi <- tryCatch(Chem$MolToSmiles(m, canonical = TRUE), error = function(e) NA_character_)
  if (is.na(smi) || !nzchar(smi)) return(NULL)
  tryCatch(Chem$MolFromSmiles(smi), error = function(e) NULL)
}

stopifnot("SMILES" %in% names(df_neuro))
smiles_vec <- clean_smiles(as.character(df_neuro$SMILES))
mols <- lapply(smiles_vec, safe_mol)
ok   <- !vapply(mols, is.null, logical(1))
cat("Parsed OK:", sum(ok), "/", length(mols), "\n")


# Morgan FP (use NumPy conversion — no 'environment' error)
nbits  <- 1024L
gen    <- rfg$GetMorganGenerator(radius = as.integer(2), fpSize = as.integer(nbits))

fp_to_bits <- function(m) {
  fp  <- gen$GetFingerprint(m)
  arr <- np$zeros(as.integer(nbits), dtype = "int8")
  DS$ConvertToNumpyArray(fp, arr)
  as.integer(py_to_r(arr))
}

smiles_to_feats <- function(m) {
  c(
    MolWt       = Descs$MolWt(m),                          # <-- use Descs
    MolLogP     = Descs$MolLogP(m),
    TPSA        = Descs$TPSA(m),
    NumHBA      = Descs$NumHAcceptors(m),
    NumHBD      = Descs$NumHDonors(m),
    NumRotBonds = Descs$NumRotatableBonds(m),
    setNames(fp_to_bits(m), paste0("FP_", seq_len(nbits)))
  )
}

feat_list <- lapply(mols[ok], smiles_to_feats)
X <- as.data.frame(do.call(rbind, feat_list))
names(X)[1:6] <- c("MolWt","MolLogP","TPSA","NumHBA","NumHBD","NumRotBonds")

df_neuro_rdkit <- df_neuro[ok, ] |>
  dplyr::bind_cols(X) |>
  dplyr::mutate(dplyr::across(where(is.character), as.factor))

view(df_neuro_rdkit)


bad_idx <- which(vapply(mols, is.null, logical(1)))
bad_tbl <- data.frame(
  row = bad_idx,
  SMILES_raw  = as.character(df_neuro$SMILES[bad_idx]),
  SMILES_clean= clean_smiles(df_neuro$SMILES[bad_idx]),
  has_gt      = grepl(">", df_neuro$SMILES[bad_idx], fixed=TRUE),   # reactions
  has_star    = grepl("\\*", df_neuro$SMILES[bad_idx]),              # SMARTS wildcard
  has_qmark   = grepl("\\?", df_neuro$SMILES[bad_idx]),
  has_dollar  = grepl("\\$", df_neuro$SMILES[bad_idx]),
  has_space   = grepl("\\s", clean_smiles(df_neuro$SMILES[bad_idx])),
  has_pipe    = grepl("\\|", df_neuro$SMILES[bad_idx]),
  nchar       = nchar(df_neuro$SMILES[bad_idx])
)
head(bad_tbl, 20)
write.csv(bad_tbl, "bad_smiles_sample.csv", row.names=FALSE)

# histogram of key descriptors
df_neuro_rdkit %>%
  pivot_longer(cols=c(MolWt,MolLogP,TPSA,NumHBA,NumHBD,NumRotBonds)) %>%
  ggplot(aes(value)) + geom_histogram(bins=40) + facet_wrap(~name, scales="free")

# PCA of descriptors
prcomp(df_neuro_rdkit[, c("MolWt","MolLogP","TPSA","NumHBA","NumHBD","NumRotBonds")], scale.=TRUE) |> biplot()

