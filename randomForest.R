library(tidyverse)
library(randomForest)
library(caret)
library(rpart)
library(rpart.plot)

# Run dataPrep_1sAnd0s.R for all.drugs 

# randomForest doesn't seem to like the `word1 word2` notation,
# so we'll transform things to `word1.word2`
names(all.drugs) = make.names(names(all.drugs))

# I will type the important ones manually, but "full" model includes this
f.numbers = paste0("fp_", 1:166)

neuro.psych = all.drugs %>%
  filter(str_detect(Disease.Area, "neurology/psychiatry"),
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


# There is a good amount of depression indication drugs in each group!

#### Model fitting

# For reference:
# Sensitivity = how good the model predicts the actually positive cases
# Specificity = how good the model predicts the actually negative cases

## Depression

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

model.depression = randomForest(
  x = depression.train[, predictors.depression],
  y = as.factor(depression.train$depression),
  data = depression.train,
  importance = TRUE,
  proximity = TRUE,
  xtest = depression.test[, predictors.depression],
  ytest = as.factor(depression.test$depression),
  keep.forest = TRUE
)

# Peek at model output & summaries to see how good it is
model.depression$importance
confusionMatrix(model.depression$predicted, model.depression$y, positive = "1")

# Confusion matrix for the testing data
depression.model.test = model.depression$test
confusionMatrix(depression.model.test$predicted, depression.test$depression, positive = "1")

## Parkinson's Disease

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

model.parkinsons = randomForest(
  x = parkinsons.train[, predictors.parkinsons],
  y = as.factor(parkinsons.train$Parkinson.s.Disease),
  data = parkinsons.train,
  importance = TRUE,
  proximity = TRUE,
  xtest = parkinsons.test[, predictors.parkinsons],
  ytest = as.factor(parkinsons.test$Parkinson.s.Disease),
  keep.forest = TRUE
)

# Peek at model output & summaries to see how good it is
model.parkinsons$importance
confusionMatrix(model.parkinsons$predicted, model.parkinsons$y, positive = "1")

# Confusion matrix for the testing data
parkinsons.model.test = model.parkinsons$test
confusionMatrix(parkinsons.model.test$predicted, parkinsons.test$Parkinson.s.Disease, positive = "1")

## Schizophrenia

set.seed(0)

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

model.schizophrenia = randomForest(
  x = schizophrenia.train[, predictors.schizophrenia],
  y = as.factor(schizophrenia.train$depression),
  data = schizophrenia.train,
  importance = TRUE,
  proximity = TRUE,
  xtest = schizophrenia.test[, predictors.schizophrenia],
  ytest = as.factor(schizophrenia.test$depression)
)

# Peek at model output & summaries to see how good it is
model.schizophrenia$importance
confusionMatrix(model.schizophrenia$predicted, model.schizophrenia$y, positive = "1")

# Confusion matrix for the testing data
schizophrenia.model.test = model.schizophrenia$test
confusionMatrix(schizophrenia.model.test$predicted, schizophrenia.test$schizophrenia, positive = "1")

### Rename fingerprints for plots

names(depression.train)[names(depression.train) == "SLC.count"] = "# of SLC Target Proteins"
names(depression.train)[names(depression.train) == "monoamine.oxidase.inhibitor"] = "Monoamine Oxidase Inhibitor"
names(depression.train)[names(depression.train) == "fp_19"] = "7M Ring"
names(depression.train)[names(depression.train) == "fp_154"] = "Presence of Carbonyl Group"
names(depression.train)[names(depression.train) == "fp_156"] = "Presence of Nitrogen Connected to 2 Other Atoms"

predictors.depression[which(predictors.depression == "SLC.count")] = "# of SLC Target Proteins"
predictors.depression[which(predictors.depression == "monoamine.oxidase.inhibitor")] = "Monoamine Oxidase Inhibitor"
predictors.depression[which(predictors.depression == "fp_19")] = "7M Ring"
predictors.depression[which(predictors.depression == "fp_154")] = "Presence of Carbonyl Group"
predictors.depression[which(predictors.depression == "fp_156")] = "Presence of Nitrogen Connected to 2 Other Atoms"

names(parkinsons.train)[names(parkinsons.train) == "SLC.count"] = "# of SLC Target Proteins"
names(parkinsons.train)[names(parkinsons.train) == "dopamine.receptor.agonist"] = "Dopamine Receptor Agonist"

names(neuro.psych)[names(neuro.psych) == "SLC.count"] = "# of SLC Target Proteins"
names(neuro.psych)[names(neuro.psych) == "dopamine.receptor.agonist"] = "Dopamine Receptor Agonist"

predictors.parkinsons[which(predictors.parkinsons == "SLC.count")] = "# of SLC Target Proteins"
predictors.parkinsons[which(predictors.parkinsons == "dopamine.receptor.agonist")] = "Dopamine Receptor Agonist"

#### Single trees using rpart

set.seed(0)

# Depression
depression.tree = rpart(
  formula = as.formula(paste0("depression ~", "`",paste0(predictors.depression, collapse = "`+`"), "`")),
                        data = depression.train[, c("depression", predictors.depression)],
  method = "class"
)

rpart.plot(depression.tree)
title(main = "Depression Tree")

# Parkinson's Disease
parkinsons.tree = rpart(
  formula = as.formula(paste0("Parkinson.s.Disease ~", "`",paste0(predictors.parkinsons, collapse = "`+`"), "`")),
  data = neuro.psych[, c("Parkinson.s.Disease", predictors.parkinsons)],
  method = "class"
)

rpart.plot(parkinsons.tree)
title(main = "Parkinson's Disease Tree")

# Schizophrenia
schizophrenia.tree = rpart(
  formula = as.formula(paste0("Parkinson.s.Disease ~", "`",paste0(predictors.schizophrenia, collapse = "`+`"), "`")),
  data = schizophrenia.train[, c("Parkinson.s.Disease", predictors.schizophrenia)],
  method = "class"
)

rpart.plot(schizophrenia.tree)
title(main = "Schizophrenia Tree")
