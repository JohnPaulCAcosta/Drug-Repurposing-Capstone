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

n = nrow(neuro.psych)
random.numbers = sample(1:n) # from 1 to length of the neuro.psych df

train.index = random.numbers[1:floor(n*.80)]
test.index = random.numbers[(floor(n*.80)+1):n]

training.drugs = neuro.psych[train.index,]
testing.drugs = neuro.psych[test.index,]

# See if there's a good spread of depression drugs

length(which(training.drugs$depression == 1))
length(which(testing.drugs$depression == 1))

# There is a good amount of depression indication drugs in each group!

#### Model fitting

# For reference:
# Sensitivity = how good the model predicts the actually positive cases
# Specificity = how good the model predicts the actually negative cases

## Depression

predictors.depression = c(
  "SLC",
  # "HTR",
  # "HRH",
  # "CHR",
  # "ADR",
  # "serotonin.reuptake.inhibitor",
  "norepinephrine.reuptake.inhibitor",
  "monoamine.oxidase.inhibitor",
  # "T.type.calcium.channel.blocker",
  # "serotonin.receptor.antagonist",
  "xlogp",
  "tpsa",
  "num.atoms"
)

model.depression = randomForest(
  x = training.drugs[, predictors.depression],
  y = as.factor(training.drugs$depression),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors.depression],
  ytest = as.factor(testing.drugs$depression)
)

# Peek at model output & summaries to see how good it is

model.depression$importance

confusionMatrix(model.depression$predicted, model.depression$y, positive = "1")

## Parkinson's Disease

predictors.parkinsons = c(
  # "HTR",
  # "DRD",
  # "ADR",
  "dopamine.receptor.agonist",
  "xlogp",
  "tpsa",
  "num.atoms"
)

model.parkinsons = randomForest(
  x = training.drugs[, predictors.parkinsons],
  y = as.factor(training.drugs$Parkinson.s.Disease),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors.parkinsons],
  ytest = as.factor(testing.drugs$Parkinson.s.Disease)
)

# Peek at model output & summaries to see how good it is

model.parkinsons$importance

confusionMatrix(model.parkinsons$predicted, model.parkinsons$y, positive = "1")

## Schizophrenia

predictors.schizophrenia = c(
  # "HTR",
  # "DRD",
  # "ADR",
  # "dopamine.receptor.antagonist",
  "xlogp",
  "tpsa",
  "num.atoms"
)

model.schizophrenia = randomForest(
  x = training.drugs[, predictors.schizophrenia],
  y = as.factor(training.drugs$depression),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors.schizophrenia],
  ytest = as.factor(testing.drugs$depression)
)

# Peek at model output & summaries to see how good it is

model.schizophrenia$importance

confusionMatrix(model.schizophrenia$predicted, model.schizophrenia$y, positive = "1")

#### What are the results of our models?

test.object = model.depression$test
votes = test.object$votes

threshold = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression, somehow I think

which(votes[,2] > threshold)

testing.drugs[which(votes[,2] > threshold), "Indication"]
testing.drugs[which(testing.drugs$depression == 1), "Indication"]

