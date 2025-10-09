library(tidyverse)
library(tidyverse)
library(broom)

# # ------------------------------------------------
# # 1) Prepare outcome
# # ------------------------------------------------
 neuro_launched <-neuro_launched %>%select(-parsed.SMILES)
 neuro_labeled <- neuro_launched %>%
  mutate(
    depression_flag = if_else(
      str_detect(tolower(coalesce(indication, "")), "\\bdepression\\b"), 1L, 0L
    )
  )%>%select(-fp)

# # ------------------------------------------------
# # 2) Build feature matrix (your helper already)
# # ------------------------------------------------
# X_neuro <- neuro_labeled
# y_neuro <- neuro_labeled$depression_flag
# 
# moa_counts_neuro <- neuro_labeled %>% dplyr::count(MOA, sort = TRUE)
# moa_top_neuro <- moa_counts_neuro %>% dplyr::slice_head(n = 20) %>% pull(MOA)
# moa_top_neuro 
# moa_counts_neuro
# 
# target_counts_neuro <- neuro_labeled %>% dplyr::count(Target, sort = TRUE)
# target_top_neuro <- target_counts_neuro %>% dplyr::slice_head(n = 20) %>% pull(Target)
# target_top_neuro 
# target_counts_neuro
# 
# 
# # ---- pick top-20 MOA and Target (UNGROUPED: split by commas first)
# moa_top_neuro <- neuro_labeled %>%
#   separate_rows(MOA, sep = ",\\s*") %>%
#   mutate(MOA = str_squish(MOA)) %>%
#   filter(MOA != "") %>%
#   dplyr::count(MOA, sort = TRUE) %>%
#   slice_head(n = 10) %>%
#   pull(MOA)
# 
# target_top_neuro <- neuro_labeled %>%
#   separate_rows(Target, sep = ",\\s*") %>%
#   mutate(Target = str_squish(Target)) %>%
#   filter(Target != "") %>%
#   dplyr::count(Target, sort = TRUE) %>%
#   slice_head(n = 10) %>%
#   pull(Target)
# 
# 
# # ---- helper: case-insensitive multi-hot with canonical names ----
# make_multi_hot_ci <- function(df, col, tokens, prefix) {
#   if (length(tokens) == 0) return(tibble())
#   canon <- tibble(tok = tokens, tok_lc = tolower(tokens))
#   
#   long <- df %>%
#     mutate(.row_id = row_number()) %>%
#     tidyr::separate_rows({{ col }}, sep = ",\\s*") %>%
#     mutate(.tok = str_squish({{ col }}),
#            .tok_lc = tolower(.tok)) %>%
#     filter(.tok_lc != "")
#   
#   long_map <- long %>%
#     inner_join(canon, by = c(".tok_lc" = "tok_lc")) %>%
#     transmute(.row_id,
#               key = paste0(prefix, janitor::make_clean_names(tok)),
#               val = 1L) %>%
#     distinct()
#   
#   out <- tibble(.row_id = seq_len(nrow(df))) %>%
#     left_join(long_map, by = ".row_id") %>%
#     select(-.row_id) %>%
#     tidyr::pivot_wider(
#       names_from  = key,
#       values_from = val,
#       # IMPORTANT: use a *named list* for older tidyr versions
#       values_fill = list(val = 0L)
#     )
#   
#   if (ncol(out) == 0) tibble() else dplyr::mutate(out, dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0L)))
# }
# 
# # Re-run these (note: it's "moa_bin", not "oa_bin")
# moa_bin  <- make_multi_hot_ci(neuro_labeled, MOA,    tokens = moa_top_neuro,    prefix = "moa_")
# targ_bin <- make_multi_hot_ci(neuro_labeled, Target, tokens = target_top_neuro, prefix = "targ_")
# 
# 
# # ---- outcome on launched ----
# neuro_labeled <- neuro_launched %>%
#   mutate(depression = as.integer(str_detect(tolower(coalesce(indication, "")), "\\bdepression\\b")))
# 
# # ---- ensure fingerprints exist ----
# neuro_labeled <- neuro_labeled %>%
#   mutate(first.SMILES = if (!"first.SMILES" %in% names(.)) purrr::map_chr(SMILES, ~ str_split(.x, ",\\s*")[[1]][[1]]) else first.SMILES)
# 
# if (!("parsed.SMILES" %in% names(neuro_labeled)) || any(purrr::map_lgl(neuro_labeled$parsed.SMILES, is.null))) {
#   neuro_labeled <- neuro_labeled %>%
#     mutate(parsed.SMILES = purrr::map(first.SMILES, ~ rcdk::parse.smiles(.x)[[1]])) %>%
#     filter(!purrr::map_lgl(parsed.SMILES, is.null))
# }
# if (!("fp" %in% names(neuro_labeled)) || all(purrr::map_lgl(neuro_labeled$fp, is.null))) {
#   neuro_labeled <- neuro_labeled %>%
#     mutate(fp = purrr::map(parsed.SMILES, ~ fingerprint::get.fingerprint(.x, type = "extended")))
# }
# 
# # ---- BUILD FEATURES: include physchem + MOA/Target + ALL fp_* ----
# physchem <- neuro_labeled %>%
#   mutate(across(c(mass, num.atoms, xlogp, tpsa), ~tidyr::replace_na(.x, 0))) %>%
#   select(mass, num.atoms, xlogp, tpsa)
# 
# moa_bin  <- make_multi_hot_ci(neuro_labeled, MOA,    tokens = moa_top_neuro,    prefix = "moa_")
# targ_bin <- make_multi_hot_ci(neuro_labeled, Target, tokens = target_top_neuro, prefix = "targ_")
# fp_df    <- fp_to_df(neuro_labeled$fp)
# 
# X_glm_full <- bind_cols(physchem, moa_bin, targ_bin, fp_df)
# 
# # ---- ONLY trim near-zero variance on fingerprint bits (keep physchem & one-hots) ----
# fp_cols <- grepl("^fp_", colnames(X_glm_full))
# nzv_idx <- if (any(fp_cols)) caret::nearZeroVar(X_glm_full[, fp_cols, drop = FALSE]) else integer(0)
# if (length(nzv_idx)) {
#   # map nzv indices back to full colnames for fp subset
#   drop_names <- colnames(X_glm_full)[fp_cols][nzv_idx]
#   X_glm_full <- X_glm_full[, setdiff(colnames(X_glm_full), drop_names), drop = FALSE]
# }
# 
# # ---- sanity check: you should now SEE these in names() ----
# # print(colnames(X_glm_full)[1:50])  # optional
# # stopifnot(all(c("mass","tpsa","xlogp","num.atoms") %in% colnames(X_glm_full)))
# 
# # ---- fit the GLM with ALL requested features ----
# glm_df <- bind_cols(X_glm_full, depression = neuro_labeled$depression)
# glm_fit_all <- glm(depression ~ ., data = glm_df, family = binomial())
# 
# summary(glm_fit_all)  # should now show mass/xlogp/tpsa + moa_* / targ_* + many fp_*
# 
# 
# 
# model.example = glm(
#   data = neuro_launched,
#   family = "binomial",
#   formula = depression ~.)
# 
# #### Look at model output & summaries
# 
# summary(model.example)
# 
# 
# 









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

write.csv(scored, "non_launched_neuro_predictions.csv", row.names = FALSE)
cat("Saved: non_launched_neuro_predictions.csv\n")

# optional: RF top importance
rf_imp <- importance(rf_fit)
rf_imp_top <- as.data.frame(rf_imp) %>%
  rownames_to_column("feature") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 25)
print(rf_imp_top)



library(tidyverse)
library(randomForest)
library(caret)

?randomForest

# Run dataPrep_1sAnd0s.R for all.drugs 

# randomForest doesn't seem to like the `word1 word2` notation,
# so we'll transform things to `word1.word2`
names(all.drugs) = make.names(names(all.drugs))

neuro.psych = all.drugs %>%
  filter(str_detect(Disease.Area, "neurology/psychiatry"),
         Phase == "Launched")

#### Train/test splitting

set.seed(0)
names(neuro_launched) = make.names(names(neuro_launched))
n = nrow(neuro_launched)
random.numbers = sample(1:n) # from 1 to length of the neuro.psych df

train.index = random.numbers[1:floor(n*.80)]
test.index = random.numbers[(floor(n*.80)+1):n]

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
predictors = c(
  "SLC",
   "HTR",
  "HRH",
   "CHR",
  "ADR",
  "serotonin.reuptake.inhibitor",
  "norepinephrine.reuptake.inhibitor",
  "monoamine.oxidase.inhibitor",
   "sodium.channel.blocker"                              
  ,"adrenergic.receptor.agonist"                         
  , "adrenergic.receptor.antagonist"                      
  , "dopamine.receptor.agonist"                           
  , "histamine.receptor.antagonist"                       
  , "monoamine.oxidase.inhibitor"                         
  , "glutamate.receptor.antagonist"                       
                      
  , "acetylcholinesterase.inhibitor"                      
  , "glucocorticoid.receptor.agonist"                     
  , "phosphodiesterase.inhibitor"                         
  , "norepinephrine.reuptake.inhibitor",                 
  "DRD",
  "xlogp",
  "tpsa",
  "num.atoms"
)

model.example = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$depression),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$depression)
)

model.example2 = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$Parkinson.s.Disease),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$Parkinson.s.Disease)
)


model.example3 = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$schizophrenia),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$schizophrenia)
)

#### Look at model output & summaries

summary(model.example)
model.example$importance

# For reference:
# Sensitivity = how good the model predicts the actually positive cases
# Specificity = how good the model predicts the actually negative cases

confusionMatrix(model.example$predicted, model.example$y, positive = "1")

test.object = model.example$test
votes = test.object$votes

threshold = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes[,2] > threshold)

testing.drugs[which(votes[,2] > threshold), "indication"]
testing.drugs[which(testing.drugs$schizophrenia == 1), "indication"]

