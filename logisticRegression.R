library(tidyverse)
library(pROC)
library(caret)

?glm

# Run dataPrep_1sAnd0s.R for all.drugs 
names(all.drugs) = make.names(names(all.drugs))

neuro.psych = all.drugs %>%
  filter(str_detect(`Disease.Area`, "neurology/psychiatry"),
         Phase == "Launched")

depression.data = neuro.psych %>%
  filter(str_detect(Indication, "depression"))
inverse.depression.data = neuro.psych %>%
  filter(!str_detect(Indication, "depression"))

parkinsons.data = neuro.psych %>%
  filter(str_detect(Indication, "Parkinson's Disease"))
inverse.parkinsons.data = neuro.psych %>%
  filter(!str_detect(Indication, "Parkinson's Disease"))

schizophrenia.data = neuro.psych %>%
  filter(str_detect(Indication, "schizophrenia"))
inverse.schizophrenia.data = neuro.psych %>%
  filter(!str_detect(Indication, "schizophrenia"))

#### Train/test splitting

train.percentage = 0.70

set.seed(1)

n.true = nrow(depression.data)
random.numbers.true = sample(1:n.true)
n.false = nrow(inverse.depression.data)
random.numbers.false = sample(1:n.false)

depression.train = rbind(
  depression.data[random.numbers.true[1:floor(n.true*train.percentage)],],
  inverse.depression.data[random.numbers.false[1:floor(n.false*train.percentage)],]
)

depression.test = rbind(
  depression.data[random.numbers.true[(floor(n.true*train.percentage)+1):n.true],],
  inverse.depression.data[random.numbers.false[(floor(n.false*train.percentage)+1):n.false],]
)

n.true = nrow(parkinsons.data)
random.numbers.true = sample(1:n.true)
n.false = nrow(inverse.parkinsons.data)
random.numbers.false = sample(1:n.false)

parkinsons.train = rbind(
  parkinsons.data[random.numbers.true[1:floor(n.true*train.percentage)],],
  inverse.parkinsons.data[random.numbers.false[1:floor(n.false*train.percentage)],]
)

parkinsons.test = rbind(
  parkinsons.data[random.numbers.true[(floor(n.true*train.percentage)+1):n.true],],
  inverse.parkinsons.data[random.numbers.false[(floor(n.false*train.percentage)+1):n.false],]
)

n.true = nrow(schizophrenia.data)
random.numbers.true = sample(1:n.true)
n.false = nrow(inverse.schizophrenia.data)
random.numbers.false = sample(1:n.false)

schizophrenia.train = rbind(
  schizophrenia.data[random.numbers.true[1:floor(n.true*train.percentage)],],
  inverse.schizophrenia.data[random.numbers.false[1:floor(n.false*train.percentage)],]
)

schizophrenia.test = rbind(
  schizophrenia.data[random.numbers.true[(floor(n.true*train.percentage)+1):n.true],],
  inverse.schizophrenia.data[random.numbers.false[(floor(n.false*train.percentage)+1):n.false],]
)

#### Model fitting starting with features found in RF

set.seed(1)

predictors.depression = c(
  # "SLC",
  "SLC.count",
  # "HTR",
  "HTR.count",
  # "HRH",
  # "HRH.count",
  # "CHR",
  # "CHR.count",
  # "ADR",
  # "ADR.count",
  # "serotonin.reuptake.inhibitor",
  "serotonin.receptor.agonist",
  # "serotonin.receptor.antagonist",
  # "norepinephrine.reuptake.inhibitor",
  # "monoamine.oxidase.inhibitor",
  # "T.type.calcium.channel.blocker",
  # "serotonin.receptor.antagonist",
  "xlogp",
  "tpsa" ,
  "num.atoms",
  # f.numbers # ,
  "fp_19", # presence of a 7 membered ring
  # "fp_58",
  "fp_82", # presence of a methylene group bridging together a two heavy atoms, one of which is bonded to hydrogen
  "fp_104", # presence of two heavys atom bonded to hydrogen connected by methylene and a heavy atom 
  "fp_144", # presence of three non-aromatic atoms chained together
  "fp_154", # presence of a carbonyl group
  "fp_156" # presence of a nitrogen atom bonded to three heavy atoms
)

predictors.parkinsons = c(
  # "HTR",
  # "HTR.count",
  # "DRD",
  # "DRD.count",
  # "ADR",
  # "ADR.count",
  "dopamine.receptor.agonist",
  "xlogp",
  "tpsa",
  "num.atoms",
  # f.numbers # ,
  "fp_17" #, # presence of a carbon atom triple bonded to another carbon atom
  # "fp_54",
  # "fp_91",
  # "fp_96" #,
  # "fp_155",
  # "fp_133"
)

predictors.schizophrenia = c(
  "HTR",
  "HTR.count",
  "HRH",
  "HRH.count",
  "DRD",
  "DRD.count",
  # "CHR",
  # "CHR.count",
  # "ADR",
  # "ADR.count",
  "dopamine.receptor.antagonist",
  "serotonin.receptor.antagonist",
  "xlogp",
  "tpsa",
  "num.atoms",
  "mass",
  # f.numbers # ,
  # "fp_19",
  # "fp_41",
  # "fp_57",
  # "fp_66",
  # "fp_74",
  # "fp_75",
  # "fp_82",
  # "fp_93",
  # "fp_95",
  # "fp_98",
  # "fp_104",
  # "fp_105",
  # "fp_111",
  # "fp_113",
  # "fp_120",
  # "fp_121",
  # "fp_128",
  # "fp_144",
  # "fp_147",
  "fp_151"
  # "fp_154",
  # "fp_156"
)

## Model fitting & ROC curves

# Depression

depression.cutoff = 0.7 # Empirically, I found that this was the best for testing,
# though from 0.5 and beyond the model did well

depression.full = glm(
  data = depression.train,
  family = "binomial",
  formula = as.formula(paste0("depression ~", "`",paste0(predictors.depression, collapse = "`+`"), "`"))
)

test.predictions = predict(
  object = depression.full,
  newdata = depression.test,
  type = "response"
)

# Training Confusion Matrix
confusionMatrix(as.factor(as.numeric(depression.full$fitted.values > depression.cutoff)), depression.train$depression, positive = "1")
# Testing Confusion Matrix
confusionMatrix(as.factor(as.numeric(test.predictions > depression.cutoff)), depression.test$depression, positive = "1")

summary(depression.full)

# Subset predictors for which p-value < .1 in fuller model to get an example
# of a subset of predictors
# 
# depression.sub = glm(
#   data = depression.train,
#   family = "binomial",
#   formula = as.formula(paste0("depression ~", "`",paste0(c(
#     "SLC.count",
#     "tpsa",
#     "fp_19",
#     "fp_82",
#     "fp_104"
#   ), collapse = "`+`"), "`"))
# )
# 
# test.predictions = predict(
#   object = depression.sub,
#   newdata = depression.test,
#   type = "response"
# )
# 
# roc.curve = roc(depression.test$depression, test.predictions)
# plot(roc.curve, col = "maroon", main = "ROC Curve", print.auc = TRUE)
# 
# # Training Confusion Matrix
# confusionMatrix(as.factor(as.numeric(depression.sub$fitted.values > depression.cutoff)), depression.train$depression, positive = "1")
# # Testing Confusion Matrix
# confusionMatrix(as.factor(as.numeric(test.predictions > depression.cutoff)), depression.test$depression, positive = "1")
# 
# summary(depression.sub)
# 
# roc.curve = roc(depression.test$depression, test.predictions)
# plot(roc.curve, col = "maroon", main = "ROC Curve",
#      print.auc = TRUE)


# Parkinson's Disease

parkinsons.cutoff = 0.5

parkinsons.full = glm(
  data = parkinsons.train,
  family = "binomial",
  formula = as.formula(paste0("Parkinson.s.Disease ~", "`",paste0(predictors.parkinsons, collapse = "`+`"), "`"))
)

test.predictions = predict(
  object = parkinsons.full,
  newdata = parkinsons.test
)

roc.curve = roc(parkinsons.test$depression, test.predictions)
plot(roc.curve, col = "maroon", main = "ROC Curve", 
     print.auc = TRUE)

# Training Confusion Matrix
confusionMatrix(as.factor(as.numeric(parkinsons.full$fitted.values > parkinsons.cutoff)), parkinsons.train$Parkinson.s.Disease, positive = "1")
# Testing Confusion Matrix
confusionMatrix(as.factor(as.numeric(test.predictions > parkinsons.cutoff)), parkinsons.test$Parkinson.s.Disease, positive = "1")

summary(parkinsons.full)

# parkinsons.sub = glm(
#   data = parkinsons.train,
#   family = "binomial",
#   formula = as.formula(paste0("Parkinson.s.Disease ~", "`",paste0(c(
#     "dopamine.receptor.agonist",
#     "tpsa",
#     "num.atoms"
#   ), collapse = "`+`"), "`"))
# )
# 
# test.predictions = predict(
#   object = parkinsons.sub,
#   newdata = parkinsons.test
# )
# 
# roc.curve = roc(parkinsons.test$depression, test.predictions)
# plot(roc.curve, col = "maroon", main = "ROC Curve", 
#      print.auc = TRUE)
# 
# # Training Confusion Matrix
# confusionMatrix(as.factor(as.numeric(parkinsons.sub$fitted.values > parkinsons.cutoff)), parkinsons.train$Parkinson.s.Disease, positive = "1")
# # Testing Confusion Matrix
# confusionMatrix(as.factor(as.numeric(test.predictions > parkinsons.cutoff)), parkinsons.test$Parkinson.s.Disease, positive = "1")
# 
# summary(parkinsons.sub)

# Schizophrenia

schizophrenia.cutoff = 0.6

schizophrenia.full = glm(
  data = schizophrenia.train,
  family = "binomial",
  formula = as.formula(paste0("schizophrenia ~", "`",paste0(predictors.schizophrenia, collapse = "`+`"), "`"))
)

test.predictions = predict(
  object = schizophrenia.full,
  newdata = schizophrenia.test
)

roc.curve = roc(schizophrenia.test$depression, test.predictions)
plot(roc.curve, col = "maroon", main = "ROC Curve",
     print.auc = TRUE)

# Training Confusion Matrix
confusionMatrix(as.factor(as.numeric(schizophrenia.full$fitted.values > schizophrenia.cutoff)), schizophrenia.train$schizophrenia, positive = "1")
# Testing Confusion Matrix
confusionMatrix(as.factor(as.numeric(test.predictions > schizophrenia.cutoff)), schizophrenia.test$schizophrenia, positive = "1")

summary(schizophrenia.full)

schizophrenia.sub = glm(
  data = schizophrenia.train,
  family = "binomial",
  formula = as.formula(paste0("schizophrenia ~", "`",paste0(c(
    "HTR",
    "dopamine.receptor.antagonist",
    "serotonin.receptor.antagonist",
    "num.atoms",
    "mass"
  ), collapse = "`+`"), "`"))
)

test.predictions = predict(
  object = schizophrenia.sub,
  newdata = schizophrenia.test
)

roc.curve = roc(schizophrenia.test$depression, test.predictions)
plot(roc.curve, col = "maroon", main = "ROC Curve", 
     print.auc = TRUE)

# Training Confusion Matrix
confusionMatrix(as.factor(as.numeric(schizophrenia.sub$fitted.values > schizophrenia.cutoff)), schizophrenia.train$schizophrenia, positive = "1")
# Testing Confusion Matrix
confusionMatrix(as.factor(as.numeric(test.predictions > schizophrenia.cutoff)), schizophrenia.test$schizophrenia, positive = "1")

summary(schizophrenia.sub)

