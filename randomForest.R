library(tidyverse)
library(randomForest)

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

which(training.drugs$depression == 1)
which(testing.drugs$depression == 1)

# There is a good amount of depression indication drugs in each group!

#### Model fitting (example for now)

predictors = c(
  "SLC",
  "HTR",
  "CHR",
  "ADR",
  "serotonin.reuptake.inhibitor",
  "norepinephrine.reuptake.inhibitor",
  "monoamine.oxidase.inhibitor",
  "T.type.calcium.channel.blocker",
  "serotonin.receptor.antagonist",
  "xlogp",
  "tpsa",
  "num.atoms"
)

# model.example = randomForest(
#   formula = as.factor(depression) ~ 
#     SLC + HTR + CHR + ADR +
#     serotonin.reuptake.inhibitor +
#     norepinephrine.reuptake.inhibitor +
#     monoamine.oxidase.inhibitor +
#     T.type.calcium.channel.blocker +
#     serotonin.receptor.antagonist +
#     xlogp + tpsa + num.atoms,
#   data = training.drugs,
#   importance = TRUE,
#   proximity = TRUE
#   # xtest = testing.drugs,
#   # ytest = as.factor(testing.drugs$depression),
# )

model.example = randomForest(
  x = training.drugs[, predictors],
  y = as.factor(training.drugs$depression),
  data = training.drugs,
  importance = TRUE,
  proximity = TRUE,
  xtest = testing.drugs[, predictors],
  ytest = as.factor(testing.drugs$depression)
)

#### Look at model output & summaries

summary(model.example)
model.example$importance

test.object = model.example$test
votes = test.object$votes

threshold = 0.4 # this can help us determine what probability we would say is
# "good enoough" to say a drug is able to be used for depression

which(votes[,2] > threshold)

testing.drugs[which(votes[,2] > threshold), "Indication"]
testing.drugs[which(testing.drugs$depression == 1), "Indication"]

