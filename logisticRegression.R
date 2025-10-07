library(tidyverse)

?glm

# Run dataPrep_1sAnd0s.R for all.drugs 

neuro.psych = all.drugs %>%
  filter(str_detect(`Disease Area`, "neurology/psychiatry"),
         Phase == "Launched")

#### Train/test splitting

#### Model fitting (example for now)

model.example = glm(
  data = neuro.psych,
  family = "binomial",
  formula = depression ~ 
    SLC + HTR + CHR + ADR +
    `serotonin reuptake inhibitor` +
    `norepinephrine reuptake inhibitor` +
    `monoamine oxidase inhibitor` +
    `T-type calcium channel blocker` +
    `serotonin receptor antagonist` +
    xlogp + tpsa + num.atoms
)

model.example.2 = glm(
  data = neuro.psych,
  family = "binomial",
  formula = depression ~ 
    SLC + HTR + CHR + ADR +
    `norepinephrine reuptake inhibitor` +
    `monoamine oxidase inhibitor` +
    xlogp + tpsa
)

#### Look at model output & summaries

summary(model.example)
summary(model.example.2)

# summary(all)

# Compare the AICs of the different model/parameter combos we try
model.example$aic
model.example.2$aic


