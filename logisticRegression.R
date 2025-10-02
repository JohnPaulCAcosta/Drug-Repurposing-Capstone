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
    SLC + CHR + ADR +
    `serotonin reuptake inhibitor` +
    `norepinephrine reuptake inhibitor` +
    `monoamine oxidase inhibitor` +
    xlogp + tpsa + num.atoms
)

#### Look at model output & summaries

summary(model.example)

# Compare the AICs of the different model/parameter combos we try
model.example$aic


