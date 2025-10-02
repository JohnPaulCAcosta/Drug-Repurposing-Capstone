library(tidyverse)

?glm

# Run dataPrep_1sAnd0s.R for all.drugs 

neuro.psych = all.drugs %>%
  filter(str_detect(`Disease Area`, "neurology/psychiatry"),
         Phase == "Launched")

logistic.model = glm(
  data = neuro.psych,
  family = "binomial",
  formula = depression ~ 
    SLC + CHR + ADR +
    `serotonin reuptake inhibitor` +
    `norepinephrine reuptake inhibitor` +
    `monoamine oxidase inhibitor` +
    xlogp + tpsa
)

summary(logistic.model)
