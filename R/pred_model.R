library(tidyverse)
library(tidyverse)
library(broom)
source(here("R", "clean_data.R"))

# -------------------------------
# Packages
# -------------------------------
library(tidyverse)
library(janitor)      # for make_clean_names
library(fingerprint)
library(rcdk)
library(purrr)
library(stringr)

# -------------------------------
# 0) Utilities
# -------------------------------

# Convert list of fingerprints -> 0/1 data.frame with stable "fp_*" names
fp_to_df <- function(fp_list) {
  if (length(fp_list) == 0) return(tibble())
  mat <- fingerprint::fp.to.matrix(fp_list)
  mat <- as.matrix(mat)
  colnames(mat) <- paste0("fp_", seq_len(ncol(mat)))
  as.data.frame(mat)
}

# Build binary columns for a set of tokens drawn from a multi-valued string column
# Example: col = "MOA" or "Target", tokens = c("SSRI", "MAO inhibitor", ...)
make_multi_hot <- function(df, col, tokens, prefix) {
  # separate each row on ", " and normalize spaces/case for matching
  long <- df %>%
    mutate(.row_id = row_number()) %>%
    separate_rows({{ col }}, sep = ",\\s*") %>%
    mutate(.tok = str_squish({{ col }}),
           .tok_lc = tolower(.tok))
  
  tokens_lc <- tolower(tokens)
  
  # keep only tokens we care about
  long_filt <- long %>% filter(.tok_lc %in% tokens_lc)
  
  # pivot to wide binary indicators
  wide <- long_filt %>%
    mutate(
      .tok_label = factor(.tok_lc, levels = tokens_lc, labels = tokens) # restore canonical labels
    ) %>%
    transmute(.row_id, key = paste0(prefix, janitor::make_clean_names(.tok_label)), val = 1L) %>%
    distinct() %>%
    pivot_wider(names_from = key, values_from = val, values_fill = 0L)
  
  # join back to original order
  out <- tibble(.row_id = seq_len(nrow(df))) %>%
    left_join(wide, by = ".row_id") %>%
    select(-.row_id) %>%
    mutate(across(everything(), ~replace_na(.x, 0L)))
  
  if (ncol(out) == 0) tibble() else out
}

# Align a new feature matrix to a reference (add missing cols as 0; drop extras)
align_to_ref <- function(X_new, X_ref) {
  missing_cols <- setdiff(colnames(X_ref), colnames(X_new))
  extra_cols   <- setdiff(colnames(X_new), colnames(X_ref))
  if (length(missing_cols)) X_new[, missing_cols] <- 0
  X_new <- X_new[, colnames(X_ref), drop = FALSE]
  attr(X_new, "extra_cols_dropped") <- extra_cols
  X_new
}

# -------------------------------
# 1) Select top-20 MOAs and Targets from neuro_launched
# -------------------------------
# NOTE: We use case-insensitive counting but return the most common canonical strings as seen.

# Top 20 MOAs
top20_moa <- neuro_launched %>%
  separate_rows(MOA, sep = ",\\s*") %>%
  mutate(MOA = str_squish(MOA)) %>%
  filter(MOA != "") %>%
  mutate(MOA_lc = tolower(MOA)) %>%
  group_by(MOA_lc, MOA) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  distinct(MOA_lc, .keep_all = TRUE) %>%  # one canonical per lowercase
  slice_head(n = 20) %>%
  pull(MOA)

# Top 20 Target proteins
top20_target <- neuro_launched %>%
  separate_rows(Target, sep = ",\\s*") %>%
  mutate(Target = str_squish(Target)) %>%
  filter(Target != "") %>%
  mutate(Target_lc = tolower(Target)) %>%
  group_by(Target_lc, Target) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  distinct(Target_lc, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  pull(Target)


# -------------------------------
# 2) Feature builders (TRAIN & SCORE)
# -------------------------------

# Build features for any dataframe that already has:
#   - mass, num.atoms, xlogp, tpsa (numeric)
#   - MOA (comma-separated string)
#   - Target (comma-separated string)
#   - parsed.SMILES (list of IAtomContainer)
# This will add:
#   - binary columns for top-20 MOAs (prefix "moa_")
#   - binary columns for top-20 Targets (prefix "targ_")
#   - fingerprint bits (fp_*)
#   - physchem numerics
build_features_full <- function(df, top_moa, top_targ, fp_type = "extended") {
  # Physchem
  base_feats <- df %>% 
    mutate(across(c(mass, num.atoms, xlogp, tpsa), ~replace_na(., 0))) %>%
    select(mass, num.atoms, xlogp, tpsa)
  
  # MOA / Target one-hots (top-N only)
  moa_bin   <- make_multi_hot(df, MOA,    tokens = top_moa,  prefix = "moa_")
  targ_bin  <- make_multi_hot(df, Target, tokens = top_targ, prefix = "targ_")
  
  # Fingerprints: compute if not present, else reuse
  fps <- if ("fp" %in% names(df) && !all(map_lgl(df$fp, is.null))) {
    df$fp
  } else {
    # compute fingerprints from parsed.SMILES (assumes present & non-null)
    map(df$parsed.SMILES, ~ get.fingerprint(.x, type = fp_type))
  }
  fp_df <- fp_to_df(fps)
  
  # Bind everything
  out <- bind_cols(base_feats, moa_bin, targ_bin, fp_df) %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  out
}

# -------------------------------
# 3) TRAIN features (neuro_launched)
# -------------------------------
X_train_full <- build_features_full(
  df = neuro_launched,
  top_moa = top20_moa,
  top_targ = top20_target,
  fp_type = "extended"   # could be "maccs" if you want MACCS keys instead
)

# Example outcome (you already set this in earlier code)
y_train_full <- if_else(
  str_detect(tolower(coalesce(neuro_launched$indication, "")), "\\bdepression\\b"),
  1L, 0L
)

# -------------------------------
# 4) SCORE features (non-launched)
# -------------------------------
# Ensure non_proc has same descriptors & parsed.SMILES as you built earlier
# If you haven't built them yet for the non-launched set, here is a safe snippet:
non_proc <- non_launched_with_smiles %>%
  mutate(
    first.SMILES  = map_chr(SMILES, ~ str_split(.x, ",\\s*")[[1]][[1]]),
    parsed.SMILES = map(first.SMILES, ~ parse.smiles(.x)[[1]])
  ) %>%
  filter(!map_lgl(parsed.SMILES, is.null)) %>%
  mutate(
    mass      = map_dbl(parsed.SMILES, get.exact.mass),
    num.atoms = map_dbl(parsed.SMILES, get.atom.count),
    fp        = map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")),
    xlogp     = map_dbl(parsed.SMILES, get.xlogp),
    tpsa      = map_dbl(parsed.SMILES, get.tpsa)
  )

X_score_full_raw <- build_features_full(
  df = non_proc,
  top_moa = top20_moa,
  top_targ = top20_target,
  fp_type = "extended"
)

# Align to training columns (very important!)
X_score_full <- align_to_ref(X_score_full_raw, X_train_full)

# Now X_train_full and X_score_full have the same columns,
# with the top-20 MOA/Target binaries plus all fingerprint bits.
# You can fit GLM/RF/etc. using X_train_full and score X_score_full.

# ===============================
# 5) MODELS + EVAL + SCORING
# ===============================
library(caret)
library(pROC)
library(randomForest)
library(broom)

# matrices
X_full <- X_train_full
y_all  <- y_train_full

# interpretable GLM design (drop fingerprints so p-values make sense)
glm_keep <- !grepl("^fp_", colnames(X_full))
X_glm    <- X_full[, glm_keep, drop = FALSE]

# ---- split launched -> train/valid/test
set.seed(123)
idx_tr <- createDataPartition(y_all, p = 0.60, list = FALSE)

Xg_tr <- X_glm[idx_tr, , drop = FALSE]; yg_tr <- y_all[idx_tr]
Xg_tm <- X_glm[-idx_tr, , drop = FALSE]; yg_tm <- y_all[-idx_tr]

Xf_tr <- X_full[idx_tr, , drop = FALSE]; yf_tr <- y_all[idx_tr]
Xf_tm <- X_full[-idx_tr, , drop = FALSE]; yf_tm <- y_all[-idx_tr]

idx_va <- createDataPartition(yg_tm, p = 0.50, list = FALSE)
Xg_va <- Xg_tm[idx_va, , drop = FALSE]; yg_va <- yg_tm[idx_va]
Xg_te <- Xg_tm[-idx_va, , drop = FALSE]; yg_te <- yg_tm[-idx_va]

Xf_va <- Xf_tm[idx_va, , drop = FALSE]; yf_va <- yf_tm[idx_va]
Xf_te <- Xf_tm[-idx_va, , drop = FALSE]; yf_te <- yf_tm[-idx_va]

# optional hygiene for GLM (keeps interpretability)
nzv <- nearZeroVar(Xg_tr)
if (length(nzv)) {
  Xg_tr <- Xg_tr[, -nzv, drop = FALSE]
  Xg_va <- Xg_va[, colnames(Xg_tr), drop = FALSE]
  Xg_te <- Xg_te[, colnames(Xg_tr), drop = FALSE]
}

# ---- GLM (binomial) on interpretable features
glm_df_tr <- cbind(Xg_tr, depression_flag = yg_tr)
glm_fit   <- glm(depression_flag ~ ., data = glm_df_tr, family = binomial())

# Inspect p-values in console:
print(summary(glm_fit))
glm_tidy <- broom::tidy(glm_fit) %>% arrange(p.value)
print(head(glm_tidy, 25))

glm_va_prob <- as.numeric(predict(glm_fit, newdata = Xg_va, type = "response"))
glm_te_prob <- as.numeric(predict(glm_fit, newdata = Xg_te, type = "response"))

# ---- Random Forest on full (incl. fingerprints)
set.seed(123)
rf_fit <- randomForest(
  x = Xf_tr, y = as.factor(yf_tr),
  ntree = 500,
  mtry = max(1, floor(sqrt(ncol(Xf_tr)))),
  importance = TRUE
)
rf_va_prob <- predict(rf_fit, newdata = Xf_va, type = "prob")[, "1"]
rf_te_prob <- predict(rf_fit, newdata = Xf_te, type = "prob")[, "1"]

# ---- thresholds (Youden’s J) on validation
best_thr <- function(y, p) {
  if (length(unique(p)) < 2) return(0.5)
  r <- roc(y, p, quiet = TRUE)
  as.numeric(coords(r, x = "best", best.method = "youden", ret = "threshold"))
}
thr_glm <- best_thr(yg_va, glm_va_prob)
thr_rf  <- best_thr(yf_va,  rf_va_prob)

acc_at <- function(p, y, thr) mean((p >= thr) == y)

# ---- report validation/test metrics
glm_va_auc <- if (length(unique(glm_va_prob)) > 1) as.numeric(auc(yg_va, glm_va_prob)) else NA_real_
rf_va_auc  <- if (length(unique(rf_va_prob))  > 1) as.numeric(auc(yf_va,  rf_va_prob))  else NA_real_

cat(sprintf("GLM  - Valid AUC: %s | Thr: %.3f | Acc@Thr: %.3f\n",
            ifelse(is.na(glm_va_auc),"NA",sprintf('%.3f',glm_va_auc)),
            thr_glm, acc_at(glm_va_prob, yg_va, thr_glm)))
cat(sprintf("RF   - Valid AUC: %s | Thr: %.3f | Acc@Thr: %.3f\n",
            ifelse(is.na(rf_va_auc),"NA",sprintf('%.3f',rf_va_auc)),
            thr_rf, acc_at(rf_va_prob,  yf_va,  thr_rf)))

glm_te_auc <- if (length(unique(glm_te_prob)) > 1) as.numeric(auc(yg_te, glm_te_prob)) else NA_real_
rf_te_auc  <- if (length(unique(rf_te_prob))  > 1) as.numeric(auc(yf_te,  rf_te_prob))  else NA_real_

glm_te_pred <- as.integer(glm_te_prob >= thr_glm)
rf_te_pred  <- as.integer(rf_te_prob  >= thr_rf)

cat(sprintf("GLM  - Test AUC: %s | Acc@Thr: %.3f\n",
            ifelse(is.na(glm_te_auc),"NA",sprintf('%.3f',glm_te_auc)),
            mean(glm_te_pred == yg_te)))
cat(sprintf("RF   - Test AUC: %s | Acc@Thr: %.3f\n",
            ifelse(is.na(rf_te_auc),"NA",sprintf('%.3f',rf_te_auc)),
            mean(rf_te_pred == yf_te)))

# ---- score NON-LAUNCHED
# GLM expects only the interpretable cols (no fp_). Align to GLM train columns:
Xg_score_raw <- X_score_full[, colnames(Xg_tr), drop = FALSE]  # same columns as GLM train
glm_score_prob <- as.numeric(predict(glm_fit, newdata = Xg_score_raw, type = "response"))
glm_score_pred <- as.integer(glm_score_prob >= thr_glm)

# RF uses full feature space aligned earlier:
rf_score_prob <- as.numeric(predict(rf_fit, newdata = X_score_full, type = "prob")[, "1"])
rf_score_pred <- as.integer(rf_score_prob >= thr_rf)

# attach to non_proc and save
scored <- non_proc %>%
  mutate(glm_prob = glm_score_prob,
         rf_prob  = rf_score_prob,
         avg_prob = (glm_prob + rf_prob)/2,
         glm_pred = glm_score_pred,
         rf_pred  = rf_score_pred) %>%
  arrange(desc(avg_prob))
scored_export <- dplyr::select(scored, where(~ !is.list(.)))
write.csv(scored_export, "C:/Users/Isaac/Downloads/STAT 482/non_launched_neuro_predictions3.csv", row.names = FALSE)
neuro_export <- dplyr::select(neuro_launched, where(~ !is.list(.)))
write.csv(neuro_export, "C:/Users/Isaac/Downloads/STAT 482/neuro_launched2.csv", row.names = FALSE)

cat("Saved: non_launched_neuro_predictions.csv\n")

# optional: RF top importance
rf_imp <- importance(rf_fit)
rf_imp_top <- as.data.frame(rf_imp) %>%
  rownames_to_column("feature") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 100)
print(rf_imp_top)



# -------------------------------
# Confusion Matrices (Test Set)
# -------------------------------

# Convert to factors with same positive class
glm_te_pred_factor <- factor(glm_te_pred, levels = c(0, 1))
rf_te_pred_factor  <- factor(rf_te_pred,  levels = c(0, 1))

glm_true_factor <- factor(yg_te, levels = c(0, 1))
rf_true_factor  <- factor(yf_te, levels = c(0, 1))

# Confusion matrices
cat("\n--- GLM Confusion Matrix (Test Set) ---\n")
cm_glm <- confusionMatrix(glm_te_pred_factor, glm_true_factor, positive = "1")
print(cm_glm)

cat("\n--- Random Forest Confusion Matrix (Test Set) ---\n")
cm_rf <- confusionMatrix(rf_te_pred_factor, rf_true_factor, positive = "1")
print(cm_rf)



glm_va_pred <- as.integer(glm_va_prob >= thr_glm)
rf_va_pred  <- as.integer(rf_va_prob  >= thr_rf)

cat("\n--- GLM Confusion Matrix (Validation Set) ---\n")
print(confusionMatrix(factor(glm_va_pred, levels=c(0,1)),
                      factor(yg_va, levels=c(0,1)),
                      positive="1"))

cat("\n--- Random Forest Confusion Matrix (Validation Set) ---\n")
print(confusionMatrix(factor(rf_va_pred, levels=c(0,1)),
                      factor(yf_va, levels=c(0,1)),
                      positive="1"))



library(tidyverse)
library(randomForest)
library(caret)

?randomForest

# Run dataPrep_1sAnd0s.R for all.drugs 

# randomForest doesn't seem to like the `word1 word2` notation,
# so we'll transform things to `word1.word2`
# names(all.drugs) = make.names(names(all.drugs))
# 
# neuro.psych = all.drugs %>%
#   filter(str_detect(Disease.Area, "neurology/psychiatry"),
#          Phase == "Launched")

#### Train/test splitting

set.seed(1)

n = nrow(neuro_launched)
random.numbers = sample(1:n) # from 1 to length of the neuro.psych df

train.index = random.numbers[1:floor(n*.70)]
test.index = random.numbers[(floor(n*.70)+1):n]

training.drugs = neuro_launched[train.index,]
testing.drugs = neuro_launched[test.index,]

# See if there's a good spread of depression drugs

length(which(training.drugs$depression == 1))
length(which(testing.drugs$depression == 1))

# There is a good amount of depression indication drugs in each group!

#### Model fitting (example for now)

top20_moa <- training.drugs %>%
  separate_rows(MOA, sep = ",\\s*") %>%
  mutate(MOA = str_squish(MOA)) %>%
  filter(MOA != "") %>%
  mutate(MOA_lc = tolower(MOA)) %>%
  group_by(MOA_lc, MOA) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  distinct(MOA_lc, .keep_all = TRUE) %>%  # one canonical per lowercase
  slice_head(n = 20) %>%
  pull(MOA)

# Top 20 Target proteins
top20_target <- training.drugs %>%
  separate_rows(Target, sep = ",\\s*") %>%
  mutate(Target = str_squish(Target)) %>%
  filter(Target != "") %>%
  mutate(Target_lc = tolower(Target)) %>%
  group_by(Target_lc, Target) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  distinct(Target_lc, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  pull(Target)

top20_moa
top20_target

# Build predictors automatically from top20_moa and top20_target
# MOA columns (already binary flags in your dataset)
moa_cols <- make_clean_names(tolower(top20_moa))        # "Serotonin Reuptake Inhibitor" -> "serotonin_reuptake_inhibitor"
moa_cols <- intersect(moa_cols, make_clean_names(names(training.drugs)))  # keep existing ones
moa_cols <- names(training.drugs)[make_clean_names(names(training.drugs)) %in% moa_cols]

# Target family columns (SLC, HTR, HRH, DRD, etc.)
# columns that are actual target-family flags in your data (no ".count")
avail_targ_cols <- grep("^[A-Z]{3,4}$", names(training.drugs), value = TRUE)

# roots from your top target strings (e.g., "HTR2A" -> "HTR2A" -> roots "HTR2A")
roots <- gsub("[^A-Za-z].*$", "", top20_target)   # keep leading letters
roots <- toupper(roots)

# try 4-letter root first (e.g., "GABRA1" -> "GABR"), else 3-letter (e.g., "CHRM1" -> "CHR")
cand4 <- substr(roots, 1, 4)
cand3 <- substr(roots, 1, 3)

mapped <- ifelse(cand4 %in% avail_targ_cols, cand4,
                 ifelse(cand3 %in% avail_targ_cols, cand3, NA))

targ_cols <- unique(na.omit(mapped))
targ_col_count <-paste0(targ_cols, ".count")

# Core molecular descriptors
core_feats <- intersect(c("xlogp", "tpsa", "num.atoms", "mass"), names(training.drugs))

# Final predictor list present in both train & test
predictors <- unique(c(core_feats, targ_cols, moa_cols,targ_col_count))
# add to your predictors list (whatever you already had)
predictors <- unique(c(predictors, chosen_fp))
predictors <- predictors[predictors %in% intersect(names(training.drugs), names(testing.drugs))]

# Optional: check what got picked up
print(predictors)
length(predictors)

# Clean up names safely (handles spaces, capitalization)
library(janitor)

# keep your original fp_to_df(); add this tiny wrapper
fp_to_df_fixed <- function(fp_list, ref_cols = NULL) {
  df <- fp_to_df(fp_list)  # your OG function
  if (is.null(ref_cols)) return(df)  # on TRAIN: learn columns
  align_to_ref(df, as.data.frame(matrix(0, nrow=0, ncol=length(ref_cols),
                                        dimnames=list(NULL, ref_cols))))
}

# TRAIN: learn the stable fp column names (this defines the reference)
fp_train <- fp_to_df_fixed(neuro_launched$fp)
fp_cols  <- colnames(fp_train)         # <- save these

# (optional) bind into neuro_launched so you can inspect in a data frame
neuro_launched_fp <- dplyr::bind_cols(neuro_launched, fp_train)

# example manual picks (edit these)
#chosen_fp <- c("fp_12", "fp_57", "fp_384", "fp_777", "fp_1012")

# guard against typos
#chosen_fp <- intersect(chosen_fp, fp_cols)
fp_ref_cols <- grep("^fp_", colnames(X_train_full), value = TRUE)

# Pick a few FP bits from big-model importance
# (edit head(...) or supply your own vector)
chosen_fp <- rf_imp_top$feature[grepl("^fp_", rf_imp_top$feature)]
chosen_fp <- intersect(chosen_fp, fp_ref_cols)
chosen_fp <- head(chosen_fp, 10)   # keep 4–5

# Non-fingerprint predictors you already use
nonfp_preds <- setdiff(predictors, fp_ref_cols)

# --- Build aligned FP matrices for TRAIN/TEST of small models -----------------
# Ensure fp list-columns exist; compute if missing
if (!"fp" %in% names(training.drugs) || all(purrr::map_lgl(training.drugs$fp, is.null))) {
  training.drugs <- training.drugs %>%
    dplyr::mutate(fp = purrr::map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")))
}
if (!"fp" %in% names(testing.drugs) || all(purrr::map_lgl(testing.drugs$fp, is.null))) {
  testing.drugs <- testing.drugs %>%
    dplyr::mutate(fp = purrr::map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")))
}

# Raw fp frames (your original fp_to_df)
fp_tr_small_raw <- fp_to_df(training.drugs$fp)
fp_te_small_raw <- fp_to_df(testing.drugs$fp)

# Align to big model's FP space so names/indices match
fp_ref_dummy <- as.data.frame(matrix(0, nrow = 0, ncol = length(fp_ref_cols),
                                     dimnames = list(NULL, fp_ref_cols)))
fp_tr_small <- align_to_ref(fp_tr_small_raw, fp_ref_dummy)
fp_te_small <- align_to_ref(fp_te_small_raw, fp_ref_dummy)

# --- Design matrices for the small RFs (non-fp + your chosen fp bits) --------
X_tr_small <- dplyr::bind_cols(
  training.drugs[, nonfp_preds, drop = FALSE],
  fp_tr_small[, chosen_fp, drop = FALSE]
)
X_te_small <- dplyr::bind_cols(
  testing.drugs[,  nonfp_preds, drop = FALSE],
  fp_te_small[, chosen_fp, drop = FALSE]
)







# predictors = c(
#   "SLC",
#    "HTR",
#   "HRH",
#    "CHR",
#   "ADR",
#   "serotonin.reuptake.inhibitor",
#   "norepinephrine.reuptake.inhibitor",
#   "monoamine.oxidase.inhibitor",
#    "sodium.channel.blocker"                              
#   ,"adrenergic.receptor.agonist"                         
#   , "adrenergic.receptor.antagonist"                      
#   , "dopamine.receptor.agonist"                           
#   , "histamine.receptor.antagonist"                       
#   , "monoamine.oxidase.inhibitor"                         
#   , "glutamate.receptor.antagonist"                       
#                       
#   , "acetylcholinesterase.inhibitor"                      
#   , "glucocorticoid.receptor.agonist"                     
#   , "phosphodiesterase.inhibitor"                         
#   , "norepinephrine.reuptake.inhibitor",                 
#   "DRD",
#   "xlogp",
#   "tpsa",
#   "num.atoms"
# )

## ---- auto-pick top-N engineered features (no dataset changes)




model.depression = randomForest(
  x = X_tr_small, y = as.factor(training.drugs$depression),
  xtest = X_te_small, ytest = as.factor(testing.drugs$depression),
  importance = TRUE,
  proximity = TRUE,
  keep.forest = TRUE
)

model.parkinsons = randomForest(
  x = X_tr_small, y = as.factor(training.drugs$`Parkinson's Disease`),
  xtest = X_te_small, ytest = as.factor(testing.drugs$`Parkinson's Disease`),
  importance = TRUE,
  proximity = TRUE,
  keep.forest= TRUE
)


model.schizo = randomForest(
  x = X_tr_small, y = as.factor(training.drugs$schizophrenia),
  xtest = X_te_small, ytest = as.factor(testing.drugs$schizophrenia),
  importance = TRUE,
  proximity = TRUE,
  keep.forest = TRUE
)

#### Look at model output & summaries

summary(model.depression)
model.depression$importance

summary(model.parkinsons)
model.parkinsons$importance

summary(model.schizo)
model.schizo$importance

# For reference:
# Sensitivity = how good the model predicts the actually positive cases
# Specificity = how good the model predicts the actually negative cases

model.depression$importance
depression.test = model.depression$test
confusionMatrix(depression.test$predicted, as.factor(testing.drugs$depression), positive = "1")


confusionMatrix(model.depression$predicted, model.depression$y, positive = "1")

test.object.d = model.depression$test
votes.d = test.object$votes

threshold.d = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.d[,2] > threshold.d)

testing.drugs[which(votes.d[,2] > threshold.d), "indication"]
testing.drugs[which(testing.drugs$depression == 1), "indication"]


model.schizo$importance
schizo.test = model.schizo$test
confusionMatrix(schizo.test$predicted, as.factor(testing.drugs$schizophrenia), positive = "1")



confusionMatrix(model.schizo$predicted, model.schizo$y, positive = "1")

test.object.s = model.schizo$test
votes.s = test.object.s$votes

threshold.s = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.s[,2] > threshold.s)

testing.drugs[which(votes.s[,2] > threshold.s), "indication"]
testing.drugs[which(testing.drugs$schizophrenia== 1), "indication"]



confusionMatrix(model.parkinsons$predicted, model.parkinsons$y, positive = "1")

model.parkinsons$importance
parkinsons.test = model.parkinsons$test
confusionMatrix(parkinsons.test$predicted, as.factor(testing.drugs$`Parkinson's Disease`), positive = "1")

test.object.p = model.parkinsons$test
votes.p = test.object.p$votes

threshold.p = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.p[,2] > threshold.p)

testing.drugs[which(votes.p[,2] > threshold.p), "indication"]
testing.drugs[which(testing.drugs$`Parkinson's Disease`== 1), "indication"]





# -------------------------------
# Predict on preclinical (non-launched) drugs
# -------------------------------

# Use the same predictors as your models
fp_pre_raw <- fp_to_df(non_proc$fp)                  # your existing fp_to_df()
fp_pre     <- align_to_ref(fp_pre_raw, fp_ref_dummy) # force same fp_* names/order as training

# --- attach FP cols to non_proc (no MOA/Target bin building here) ---
non_proc_scoring <- dplyr::bind_cols(
  non_proc,
  fp_pre[, fp_ref_cols, drop = FALSE]
)

# NOTE: 'predictors' should already include whatever you want:
#   - your existing protein/bin/count columns that are ALREADY in the dataset
#   - the chosen fp_* bits you picked (subset of fp_ref_cols), e.g. chosen_fp
# Make sure 'predictors' only contains columns present in both training.drugs and testing.drugs
# as you already did earlier.

# --- build the exact matrix for your small models ---
X_pre_small <- non_proc_scoring[, predictors, drop = FALSE]

# (recommended) final safety: align to the small-model train design to avoid mismatches
X_pre_small <- align_to_ref(X_pre_small, X_tr_small)

X_pre_small <- X_pre_small[, colnames(X_tr_small), drop = FALSE]

# Predict probabilities with the small models (use X_pre_small, not X_pre)
pred_depr_prob <- predict(model.depression,  newdata = X_pre_small, type = "prob")[, "1"]
pred_park_prob <- predict(model.parkinsons, newdata = X_pre_small, type = "prob")[, "1"]
pred_scz_prob  <- predict(model.schizo,     newdata = X_pre_small, type = "prob")[, "1"]

# Apply threshold (0.4 as in your test block)
thr <- 0.4
pred_depr_flag <- as.integer(pred_depr_prob >= thr)
pred_park_flag <- as.integer(pred_park_prob >= thr)
pred_scz_flag  <- as.integer(pred_scz_prob  >= thr)

# Attach results to your non-launched data
scored_preclinical <- non_launched_with_smiles %>%
  mutate(
    rf_depression_prob = pred_depr_prob,
    rf_parkinsons_prob = pred_park_prob,
    rf_schizophrenia_prob = pred_scz_prob,
    rf_depression_pred = pred_depr_flag,
    rf_parkinsons_pred = pred_park_flag,
    rf_schizophrenia_pred = pred_scz_flag
  ) %>%
  arrange(desc(rf_depression_prob))

# Export
preclinical_scored <- dplyr::select(scored_preclinical, where(~ !is.list(.)))
write.csv(preclinical_scored,
          "C:/Users/Isaac/Downloads/STAT 482/preclinical_predictions_rf3.csv")

cat("Saved: preclinical_predictions_rf3.csv\n")


scored_preclinical_top <- scored_preclinical %>%
  mutate(
    top_indication = c("Depression", "Parkinsons", "Schizophrenia")[max.col(
      cbind(rf_depression_prob, rf_parkinsons_prob, rf_schizophrenia_prob), ties.method = "first")]
  )
view(scored_preclinical_top)

varImpPlot(model.depression, main="Feature Importance – Depression RF")
ggplot(scored_preclinical, aes(x = rf_depression_prob, y = rf_schizophrenia_prob)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(title = "Overlap: Depression vs Schizophrenia Models",
       x = "Depression probability", y = "Schizophrenia probability")


# -------------------------------
# Summary Stats for Model Predictions
# -------------------------------


# Basic descriptive stats for each model
summary_stats <- scored_preclinical %>%
  summarise(
    n_compounds = n(),
    depression_mean = mean(rf_depression_prob, na.rm = TRUE),
    depression_sd   = sd(rf_depression_prob, na.rm = TRUE),
    depression_median = median(rf_depression_prob, na.rm = TRUE),
    parkinsons_mean = mean(rf_parkinsons_prob, na.rm = TRUE),
    parkinsons_sd   = sd(rf_parkinsons_prob, na.rm = TRUE),
    parkinsons_median = median(rf_parkinsons_prob, na.rm = TRUE),
    schizo_mean = mean(rf_schizophrenia_prob, na.rm = TRUE),
    schizo_sd   = sd(rf_schizophrenia_prob, na.rm = TRUE),
    schizo_median = median(rf_schizophrenia_prob, na.rm = TRUE)
  )

print(summary_stats)

# -------------------------------
# How many "high probability" hits per model
# -------------------------------

thr <- 0.6  # adjust threshold if you want
high_hits <- scored_preclinical %>%
  summarise(
    n_depression_high  = sum(rf_depression_prob  >= thr, na.rm = TRUE),
    n_parkinsons_high  = sum(rf_parkinsons_prob  >= thr, na.rm = TRUE),
    n_schizo_high      = sum(rf_schizophrenia_prob >= thr, na.rm = TRUE)
  )

print(high_hits)

# -------------------------------
# Compare model outputs
# -------------------------------

# Correlation between models (how similar their predictions are)
model_correlation <- scored_preclinical %>%
  select(rf_depression_prob, rf_parkinsons_prob, rf_schizophrenia_prob) %>%
  cor(use = "pairwise.complete.obs")

print(model_correlation)

# -------------------------------
# Top predicted compounds per model
# -------------------------------

top_depression <- scored_preclinical %>%
  arrange(desc(rf_depression_prob)) %>%
  select(DrugName = Name, rf_depression_prob) %>%
  head(10)

top_parkinsons <- scored_preclinical %>%
  arrange(desc(rf_parkinsons_prob)) %>%
  select(DrugName = Name, rf_parkinsons_prob) %>%
  head(10)

top_schizo <- scored_preclinical %>%
  arrange(desc(rf_schizophrenia_prob)) %>%
  select(DrugName =Name, rf_schizophrenia_prob) %>%
  head(10)

cat("\n--- Top 10 Predicted Depression-like Compounds ---\n")
print(top_depression)

cat("\n--- Top 10 Predicted Parkinson’s-like Compounds ---\n")
print(top_parkinsons)

cat("\n--- Top 10 Predicted Schizophrenia-like Compounds ---\n")
print(top_schizo)

##data split, how u decide on tree, what results (confusion matrix)

# After you fit once, you can prune & refit:
prune_by_importance <- function(imp_df, keep_n = 30, min_gini = .3) {
  imp_df %>%
    tibble::rownames_to_column("feature") %>%
    arrange(desc(MeanDecreaseGini)) %>%
    filter(MeanDecreaseGini >= min_gini) %>%
    slice_head(n = keep_n) %>%
    pull(feature)
}
impdep <- importance(model.depression) %>% as.data.frame()
impdep<- prune_by_importance(impdep)
impdep

imppark <- importance(model.parkinsons) %>% as.data.frame()
imppark<- prune_by_importance(imppark)
imppark

imps <- importance(model.schizo) %>% as.data.frame()
imps<- prune_by_importance(imps)
imps

# ===== Refit the three small RFs using pruned variables =====

# 0) safety: intersect pruned lists with available columns
keep_dep  <- intersect(colnames(X_tr_small), impdep)
keep_park <- intersect(colnames(X_tr_small), imppark)
keep_scz  <- intersect(colnames(X_tr_small), imps)

# 1) build pruned design matrices (train/test) — same structure as before
Xd_tr <- X_tr_small[, keep_dep,  drop = FALSE]
Xd_te <- X_te_small[, keep_dep,  drop = FALSE]

Xp_tr <- X_tr_small[, keep_park, drop = FALSE]
Xp_te <- X_te_small[, keep_park, drop = FALSE]

Xs_tr <- X_tr_small[, keep_scz,  drop = FALSE]
Xs_te <- X_te_small[, keep_scz,  drop = FALSE]



# 3) REFIT (same objects, same structure, just with pruned X_* matrices)
set.seed(1)
model.depression2 = randomForest(
  x = Xd_tr, y = as.factor(training.drugs$depression),
  xtest = Xd_te, ytest = as.factor(testing.drugs$depression),
  importance = TRUE,
  proximity = TRUE,
  keep.forest = TRUE
)

model.parkinsons2 = randomForest(
  x = Xp_tr, y = as.factor(training.drugs$`Parkinson's Disease`),
  xtest = Xp_te, ytest = as.factor(testing.drugs$`Parkinson's Disease`),
  importance = TRUE,
  proximity = TRUE,
  keep.forest = TRUE
)

model.schizo2 = randomForest(
  x = Xs_tr, y = as.factor(training.drugs$schizophrenia),
  xtest = Xs_te, ytest = as.factor(testing.drugs$schizophrenia),
  importance = TRUE,
  proximity = TRUE,
  keep.forest = TRUE
)



library(caret)
library(pROC)

# helper: pick Youden J threshold (guard for degenerate probs)
best_thr <- function(y, p) {
  if (length(unique(p)) < 2) return(0.5)
  r <- roc(y, p, quiet = TRUE)
  as.numeric(coords(r, x = "best", best.method = "youden", ret = "threshold"))
}

eval_rf <- function(model, Xtr, ytr, Xte, yte, label = "") {
  ytr <- factor(ytr, levels = c(0,1))
  yte <- factor(yte, levels = c(0,1))
  
  # ---- probabilities
  p_tr <- predict(model, newdata = Xtr, type = "prob")[, "1"]
  p_te <- predict(model, newdata = Xte, type = "prob")[, "1"]
  
  # ---- AUCs
  auc_tr <- if (length(unique(p_tr)) > 1) as.numeric(auc(ytr, p_tr)) else NA_real_
  auc_te <- if (length(unique(p_te)) > 1) as.numeric(auc(yte, p_te)) else NA_real_
  
  # ---- thresholds
  thr_te <- best_thr(yte, p_te)   # choose on TEST (you can swap to validation if you prefer)
  thr50  <- 0.5
  
  # ---- predictions @ thresholds
  pred_tr_50 <- factor(as.integer(p_tr >= thr50), levels = c(0,1))
  pred_te_50 <- factor(as.integer(p_te >= thr50), levels = c(0,1))
  
  pred_tr_bt <- factor(as.integer(p_tr >= thr_te), levels = c(0,1))
  pred_te_bt <- factor(as.integer(p_te >= thr_te), levels = c(0,1))
  
  # ---- confusion matrices
  cm_tr_50 <- confusionMatrix(pred_tr_50, ytr, positive = "1")
  cm_te_50 <- confusionMatrix(pred_te_50, yte, positive = "1")
  
  cm_tr_bt <- confusionMatrix(pred_tr_bt, ytr, positive = "1")
  cm_te_bt <- confusionMatrix(pred_te_bt, yte, positive = "1")
  
  # ---- OOB (last row of err.rate)
  oob <- tryCatch({
    er <- model$err.rate
    if (is.null(er)) NA_real_ else er[nrow(er), "OOB"]
  }, error = function(e) NA_real_)
  
  cat("\n==============================\n")
  cat(sprintf("Model: %s\n", label))
  cat("------------------------------\n")
  cat(sprintf("OOB error (rf): %s\n", ifelse(is.na(oob), "NA", sprintf("%.3f", oob))))
  cat(sprintf("AUC  Train: %s | Test: %s\n",
              ifelse(is.na(auc_tr), "NA", sprintf("%.3f", auc_tr)),
              ifelse(is.na(auc_te), "NA", sprintf("%.3f", auc_te))))
  cat(sprintf("Test-derived best threshold (Youden J): %.3f\n\n", thr_te))
  
  cat("--- Confusion Matrix @ 0.5 (TRAIN) ---\n"); print(cm_tr_50$table)
  cat(sprintf("Accuracy: %.3f  Sens: %.3f  Spec: %.3f\n\n",
              cm_tr_50$overall["Accuracy"],
              cm_tr_50$byClass["Sensitivity"],
              cm_tr_50$byClass["Specificity"]))
  
  cat("--- Confusion Matrix @ 0.5 (TEST) ---\n"); print(cm_te_50$table)
  cat(sprintf("Accuracy: %.3f  Sens: %.3f  Spec: %.3f\n\n",
              cm_te_50$overall["Accuracy"],
              cm_te_50$byClass["Sensitivity"],
              cm_te_50$byClass["Specificity"]))
  
  cat("--- Confusion Matrix @ BestThr (TRAIN) ---\n"); print(cm_tr_bt$table)
  cat(sprintf("Accuracy: %.3f  Sens: %.3f  Spec: %.3f\n\n",
              cm_tr_bt$overall["Accuracy"],
              cm_tr_bt$byClass["Sensitivity"],
              cm_tr_bt$byClass["Specificity"]))
  
  cat("--- Confusion Matrix @ BestThr (TEST) ---\n"); print(cm_te_bt$table)
  cat(sprintf("Accuracy: %.3f  Sens: %.3f  Spec: %.3f\n",
              cm_te_bt$overall["Accuracy"],
              cm_te_bt$byClass["Sensitivity"],
              cm_te_bt$byClass["Specificity"]))
  
  invisible(list(
    auc_train = auc_tr, auc_test = auc_te,
    thr_test = thr_te,
    cm_train_05 = cm_tr_50, cm_test_05 = cm_te_50,
    cm_train_bt = cm_tr_bt, cm_test_bt = cm_te_bt
  ))
}

# ==== run evals on your three NEW models ====
res_dep  <- eval_rf(model.depression2,  Xd_tr, training.drugs$depression,            Xd_te, testing.drugs$depression,            label = "Depression")
res_park <- eval_rf(model.parkinsons2,  Xp_tr, training.drugs$`Parkinson's Disease`,  Xp_te, testing.drugs$`Parkinson's Disease`,  label = "Parkinson's")
res_scz  <- eval_rf(model.schizo2,      Xs_tr, training.drugs$schizophrenia,          Xs_te, testing.drugs$schizophrenia,          label = "Schizophrenia")



model.depression2$importance
depression2.test = model.depression2$test
confusionMatrix(depression2.test$predicted, as.factor(testing.drugs$depression), positive = "1")


confusionMatrix(model.depression2$predicted, model.depression2$y, positive = "1")

test.object.d2 = model.depression2$test
votes.d2 = test.object.d2$votes

threshold.d2 = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.d2[,2] > threshold.d2)

testing.drugs[which(votes.d2[,2] > threshold.d2), "indication"]
testing.drugs[which(testing.drugs$depression == 1), "indication"]


model.schizo2$importance
schizo.test2 = model.schizo2$test
confusionMatrix(schizo.test2$predicted, as.factor(testing.drugs$schizophrenia), positive = "1")



confusionMatrix(model.schizo2$predicted, model.schizo2$y, positive = "1")

test.object.s2 = model.schizo2$test
votes.s2 = test.object.s2$votes

threshold.s2 = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.s2[,2] > threshold.s2)

testing.drugs[which(votes.s2[,2] > threshold.s2), "indication"]
testing.drugs[which(testing.drugs$schizophrenia== 1), "indication"]



confusionMatrix(model.parkinsons2$predicted, model.parkinsons2$y, positive = "1")

model.parkinsons2$importance
parkinsons.test2 = model.parkinsons2$test
confusionMatrix(parkinsons.test2$predicted, as.factor(testing.drugs$`Parkinson's Disease`), positive = "1")

test.object.p2 = model.parkinsons2$test
votes.p2 = test.object.p2$votes

threshold.p2 = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes.p2[,2] > threshold.p2)

testing.drugs[which(votes.p2[,2] > threshold.p2), "indication"]
testing.drugs[which(testing.drugs$`Parkinson's Disease`== 1), "indication"]


