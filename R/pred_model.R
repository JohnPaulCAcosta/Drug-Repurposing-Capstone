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

# Clean up names safely (handles spaces, capitalization)
library(janitor)

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


# Core molecular descriptors
core_feats <- intersect(c("xlogp", "tpsa", "num.atoms", "mass"), names(training.drugs))

# Final predictor list present in both train & test
predictors <- unique(c(core_feats, targ_cols, moa_cols))
predictors <- predictors[predictors %in% intersect(names(training.drugs), names(testing.drugs))]

# Optional: check what got picked up
print(predictors)
length(predictors)


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
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$depression),
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$depression),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
)

model.parkinsons = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$`Parkinson's Disease`),
  data = training.drugs,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$`Parkinson's Disease`),
  importance = TRUE,
  proximity = TRUE,
  keep.forest= TRUE
)


model.schizo = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$schizophrenia),
  data = training.drugs,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$schizophrenia),
  importance = TRUE,
  proximity = TRUE,
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
vote.ds = test.object$votes

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

which(votes.p[,2] > threshold.s)

testing.drugs[which(votes.p[,2] > threshold.p), "indication"]
testing.drugs[which(testing.drugs$`Parkinson's Disease`== 1), "indication"]





# -------------------------------
# Predict on preclinical (non-launched) drugs
# -------------------------------

# Use the same predictors as your models
X_pre <- non_proc[, predictors, drop = FALSE]

# Predict probabilities
pred_depr_prob <- predict(model.depression,  newdata = X_pre, type = "prob")[, "1"]
pred_park_prob <- predict(model.parkinsons, newdata = X_pre, type = "prob")[, "1"]
pred_scz_prob  <- predict(model.schizo,     newdata = X_pre, type = "prob")[, "1"]

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