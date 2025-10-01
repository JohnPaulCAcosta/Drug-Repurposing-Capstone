library(tidyverse)

?glm

new = read.csv("allDrugs1s0s.csv")

parkinsons = all.drugs %>%
  filter(`Parkinson's Disease` == 1,
         !is.na(SMILES),
         !is.na(MOA))

logistic.model = glm(
  data = parkinsons,
  family = "binomial",
  formula = `Parkinson's Disease` ~ SLC
)
