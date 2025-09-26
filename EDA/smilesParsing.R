library(tidyverse)
library(readxl)
library(rcdk)
library(fingerprint)

all.drugs = read_xls("literallyAllDrugs.xls") %>%
  filter(!is.na(SMILES)) %>%
  mutate(first.SMILES = map_chr(SMILES, ~ str_split(.x, ",\\s*")[[1]][[1]])) %>%
  mutate(parsed.SMILES = map(first.SMILES, ~ parse.smiles(.x)[[1]])) %>%
  filter(!map_lgl(parsed.SMILES, is.null)) %>%
  mutate(
    mass = map_dbl(parsed.SMILES, get.exact.mass),
    num.atoms = map_dbl(parsed.SMILES, get.atom.count),
    # bonds = map_dbl(parsed.SMILES, get.bonds), # doesn't work
    fp = map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")),
    xlogp = map_dbl(parsed.SMILES, get.xlogp).
    get.mol2formula()
  )

test = rcdk::get.fingerprint(all.drugs$parsed.SMILES[[1]])

summary(all.drugs[, c("mass", "num.atoms", "xlogp")])

# Mass, # of Atoms, and XLogP across all Disease Areas

all.drugs %>%
  separate_rows(`Disease Area`, sep = ",\\s*") %>%
  group_by(`Disease Area`) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = `Disease Area`, x = mass)) +
  geom_boxplot() 

all.drugs %>%
  separate_rows(`Disease Area`, sep = ",\\s*") %>%
  group_by(`Disease Area`) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = `Disease Area`, x = num.atoms)) +
  geom_boxplot() 

all.drugs %>%
  separate_rows(`Disease Area`, sep = ",\\s*") %>%
  group_by(`Disease Area`) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = `Disease Area`, x = xlogp)) +
  geom_boxplot() 

# Mass, # of Atoms, and XLogP across neuro/psych drugs

neuro.psych = all.drugs %>%
  filter(`Disease Area` == "neurology/psychiatry", Phase == "Launched")

summary(neuro.psych[,c("mass", "num.atoms", "xlogp")])

neuro.psych %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  group_by(Indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = Indication, x = mass)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Exact Mass")

neuro.psych %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  group_by(Indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = Indication, x = num.atoms)) +
  geom_boxplot() +
  theme_bw() +
  xlab("# of Atoms")

neuro.psych %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  group_by(Indication) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  ggplot(aes(y = Indication, x = xlogp)) +
  geom_boxplot() +
  theme_bw() +
  xlab("XLogP")
