library(here)
library(tidyverse)
library(rcdk)
library(rJava)
library(fingerprint)
source(here("R", "feature_engineering.R"))
# all_drugs <- read.csv(here("drug_data", "literallyAllDrugs.csv"),
#                       check.names = FALSE)

all_drugs<-new_dat

all_drugs <- all_drugs %>%
  mutate(
    phase       = str_to_lower(str_squish(Phase)),
    disease_area = str_to_lower(str_squish(`Disease Area`)),
    indication   = str_squish(Indication)
  ) %>% select(-Phase, -`Disease Area`, -Indication)

neuro_launched <- all_drugs %>%
  filter(
    str_detect(coalesce(disease_area, ""), fixed("neurology/psychiatry", ignore_case = TRUE)) &
      phase == "launched"
  )
nrow(neuro_launched)

non_launched_all <- all_drugs %>%
  filter(is.na(phase) | phase != "launched")
nrow(non_launched_all)

non_launched_no_ind <- non_launched_all %>%
     filter(is.na(indication) | indication== "")
nrow(non_launched_no_ind)

  # Optional: if your modeling depends on chemical structure, keep only rows with SMILES
 non_launched_with_smiles <- non_launched_all %>% filter(!is.na(SMILES) & SMILES != "")
nrow(non_launched_with_smiles) 

neuro_launched <- neuro_launched %>% filter(!is.na(SMILES) & SMILES != "" & !is.na(Target) & Target != "" & !is.na(MOA) & MOA != "")
nrow(neuro_launched)




neuro_launched <- neuro_launched %>%
  filter(
    !is.na(SMILES), 
    !is.na(Target), 
    !is.na(MOA)
  ) %>%
  mutate(first.SMILES = map_chr(SMILES, ~ str_split(.x, ",\\s*")[[1]][[1]])) %>%
  mutate(parsed.SMILES = map(first.SMILES, ~ parse.smiles(.x)[[1]])) %>%
  filter(!map_lgl(parsed.SMILES, is.null)) %>%
  mutate(
    mass = map_dbl(parsed.SMILES, get.exact.mass),
    num.atoms = map_dbl(parsed.SMILES, get.atom.count),
    # bonds = map_dbl(parsed.SMILES, get.bonds), # doesn't work
    fp = map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")),
    xlogp = map_dbl(parsed.SMILES, get.xlogp),
    tpsa = map_dbl(parsed.SMILES, get.tpsa)
  )

view(neuro_launched)

summary(neuro_launched[, c("mass", "num.atoms", "xlogp")])

# Mass, # of Atoms, and XLogP across all Disease Areas

neuro_launched %>%
  separate_rows(disease_area, sep = ",\\s*") %>%
  group_by(disease_area) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = disease_area, x = mass)) +
  geom_boxplot() 

neuro_launched %>%
  separate_rows(disease_area, sep = ",\\s*") %>%
  group_by(disease_area) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = disease_area, x = num.atoms)) +
  geom_boxplot() 

neuro_launched %>%
  separate_rows(disease_area, sep = ",\\s*") %>%
  group_by(disease_area) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = disease_area, x = xlogp)) +
  geom_boxplot() 

# I mean I guess they're different for each disease area

# Mass, # of Atoms, and XLogP across neuro/psych drugs


neuro_launched %>%
  separate_rows(indication, sep = ",\\s*") %>%
  group_by(indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = indication, x = mass)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Exact Mass")

neuro_launched %>%
  separate_rows(indication, sep = ",\\s*") %>%
  group_by(indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = indication, x = num.atoms)) +
  geom_boxplot() +
  theme_bw() +
  xlab("# of Atoms")

neuro_launched %>%
  separate_rows(indication, sep = ",\\s*") %>%
  group_by(indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = indication, x = xlogp)) +
  geom_boxplot() +
  theme_bw() +
  xlab("XLogP")

# Both summary tables
test.smiles.p = parse.smiles(str_split(neuro_launched$SMILES[[1]], ", ")[[1]][1])
test.smiles.p

moa_counts_neuro <- neuro_launched %>% dplyr::count(MOA, sort = TRUE)
moa_top_neuro <- moa_counts_neuro %>% dplyr::slice_head(n = 20) %>% pull(MOA)
moa_top_neuro 
moa_counts_neuro

target_counts_neuro <- neuro_launched %>% dplyr::count(Target, sort = TRUE)
target_top_neuro <- target_counts_neuro %>% dplyr::slice_head(n = 20) %>% pull(Target)
target_top_neuro 
target_counts_neuro

