library(tidyverse)
library(readxl)
library(rlang)
library(rcdk)

########################################################################
## Code to create counts the for dataset on the number of occurrences ##
## for each Disease Area, Target Protein, and Indication + Plotting   ##
########################################################################

all.drugs = read_xls("literallyAllDrugs.xls") %>%
  filter(
    !is.na(SMILES), 
    !is.na(Target), 
    !is.na(MOA)
  ) %>%
  mutate(first.SMILES = map_chr(SMILES, ~ str_split(.x, ",\\s*")[[1]][[1]])) %>%
  mutate(parsed.SMILES = map(first.SMILES, ~ parse.smiles(.x)[[1]])) %>%
  filter(!map_lgl(parsed.SMILES, is.null)) %>%
  mutate(
    mass = map_dbl(parsed.SMILES, get.exact.mass),
    num.atoms = map_dbl(parsed.SMILES, get.atom.count),
    # bonds = map_dbl(parsed.SMILES, get.bonds), # doesn't work
    fp = map(parsed.SMILES, ~ get.fingerprint(.x, type = "extended")),
    xlogp = map_dbl(parsed.SMILES, get.xlogp),
    tpsa = map_dbl(parsed.SMILES, get.tpsa)
  )

head(all.drugs)

## Distribution of the Targets for Launched - I do not know what this means
# Maybe it means what is seen in a patient for the treatment to be dispatched?
# Can help us choose which cancer (or health issue in general) to focus on

launched = all.drugs %>%
  filter(Phase == "Launched")

# How many are there in each disease area? (launched drugs)

launched$DAvector = strsplit(launched$`Disease Area`, ", ")

counts = list()

counter = function(vec) {
  
  for (i in vec) {
    
    if (is.null(counts[[i]])) {
      counts[[i]] <<- 0
    }
    
    counts[[i]] <<- counts[[i]] + 1
    
  }
  
}

for (vecT in launched$DAvector) {
  counter(vecT)
}

ggplot(data.frame(
  Names = names(counts),
  Counts = as.numeric(counts)
), aes(x = Counts, y = Names, fill = Names)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")

## Distribution of the Targets in Launched for Neurology/Psychiatry

nP = launched %>% # 384 launched
  filter(str_detect(`Disease Area`, "neurology/psychiatry"))

nP$TargetVector = strsplit(nP$Target, ", ")

counts = list()

counter = function(vec) {
  
  for (i in vec) {
    
    if (is.null(counts[[i]])) {
      counts[[i]] <<- 0
    }
    
    counts[[i]] <<- counts[[i]] + 1
    
  }
  
}

for (vecT in nP$TargetVector) {
  counter(vecT)
}

# Counts of the prevalent target proteins

count.filter = 7

ggplot(data.frame(
  Names = names(counts[which(counts > count.filter)]),
  Counts = as.numeric(counts[which(counts > count.filter)])
), aes(x = Counts, y = Names, fill = Names)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        legend.title = element_text("Bee")) +
  guides(fill = "none") +
  theme_bw()

# Counts by the prefixes of the Target Proteins for all launched
# neurological/psychiatric drugs

nP = launched %>% 
  filter(str_detect(`Disease Area`, "neurology/psychiatry"),
         !is.na(SMILES))

nP$TargetVector = strsplit(nP$Target, ", ")

nP$TargetPrefixVector = lapply(nP$TargetVector, function(vec) {
  unique(sapply(vec, function(x) {
    sub("^([A-Za-z]+).*", "\\1", x)
  }))
}
)

counts.pref = list()

counter.pref = function(vec) {
  
  for (i in vec) {
    
    if (is.null(counts.pref[[i]])) {
      counts.pref[[i]] <<- 0
    }
    
    counts.pref[[i]] <<- counts.pref[[i]] + 1
    
  }
  
}

for (vecT in nP$TargetPrefixVector) {
  counter.pref(vecT)
}

# Big unfiltered plot

ggplot(data.frame(
  Names = names(counts.pref),
  Counts = as.numeric(counts.pref)
), aes(x = Counts, y = Names, fill = Names)) +
  geom_bar(stat = "identity") +
  guides(fill = "none") +
  ylab("Target Protein") +
  xlab("# of Occurrences") +
  theme_bw()

# to interpret this, we can say that we are clustering target proteins
# for which the first letter grouping encompases more than
# count.filter / 370 of the launched drugs - 2 for no SMILES

count.filter = 7

# I think this is important for EDA

ggplot(data.frame(
  Names = names(counts.pref[which(counts.pref > count.filter)]),
  Counts = as.numeric(counts.pref[which(counts.pref > count.filter)])
), aes(x = Counts, y = Names)) +
  geom_bar(stat = "identity", fill = "maroon") +
  guides(fill = "none") +
  ylab("Target Protein") +
  xlab("# of Occurrences") +
  ggtitle("Distribution of Target Protein Prefixes Across Neuro/Psych Drugs") +
  theme(axis.text = element_text(size = 5)) +
  theme_bw()

########################################################################
## Code to create columns for dataset indicating whether a particular ##
## Target Protein or Indication is associated with a particular drug  ##
########################################################################

# Based on EDA, we put the popular target protein prefixes below,
# along with a comment of their occurrence count in neuro/psych drugs

popular.prefixes = c("SLC", # 50+
                     "SCN", # 25+
                     "PTGS", # 25+,
                     "MAO", # 20+
                     "HTR", # 75+
                     "HRH", # 50+
                     "GABR", # 100+
                     "DRD", # 75+,
                     "CHR", # 60+
                     "ADR") # 60+

#### Create presence and count columns for each of the popular protein groups ####

for (prefix in popular.prefixes) {
  
  all.drugs = all.drugs %>%
    rowwise() %>%
    mutate(!!prefix := as.factor(if_else(
      any(str_starts(
        str_split(Target, ", ")[[1]],
        as.character(prefix)
      )),
      1,
      0
    ))) %>%
    mutate(!!paste0(prefix, ".count") := sum(str_starts(
      str_split(Target, ", ")[[1]],
      as.character(prefix)
    ))
    ) %>%
    ungroup()
  
}

#### Create presence and count columns for each of the popular indication groups groups ####

# at this time, popular indications included are these, but
# can be recalculated at the bottom

popular.indication = c("Parkinson's Disease", "depression",
                       "pain relief", "schizophrenia")

for (ind in popular.indication) {
  
  all.drugs = all.drugs %>%
    rowwise() %>%
    mutate(!!paste0(ind) := as.factor(sum(str_starts(
      str_split(Indication, ", ")[[1]],
      as.character(ind)
    )))
    ) %>%
    ungroup()
  
}

# sum(all.drugs$`Parkinson's Disease`)

all.drugs$`Parkinson's Disease` = if_else(
  is.na(all.drugs$`Parkinson's Disease`),
  as.factor(0),
  all.drugs$`Parkinson's Disease`
)

all.drugs$`depression` = if_else(
  is.na(all.drugs$`depression`),
  as.factor(0),
  all.drugs$`depression`
)

all.drugs$`pain relief` = if_else(
  is.na(all.drugs$`pain relief`),
  as.factor(0),
  all.drugs$`pain relief`
)

all.drugs$`schizophrenia` = if_else(
  is.na(all.drugs$`schizophrenia`),
  as.factor(0),
  all.drugs$`schizophrenia`
)

# sum(all.drugs$schizophrenia) # should NOT be the same
# sum(!is.na(all.drugs$Indication)) # should be the same

#### Create presence and count columns for each of the MOAs groups ####

# at this time, popular indications included are these, but
# can be recalculated at the bottom

np.moas = unique((nP %>%
                    separate_rows(MOA, sep = ",\\s*") %>%
                    select(MOA))$MOA)

for (moa in np.moas) {
  
  all.drugs = all.drugs %>%
    rowwise() %>%
    mutate(!!paste0(moa) := as.factor(sum(str_starts(
      str_split(MOA, ", ")[[1]],
      as.character(moa)
    )))) %>%
    ungroup()
  
}


# all.drugs %>%
#   select(-first.SMILES, -parsed.SMILES) %>%
#   write.csv(file = "allDrugs1s0s.csv")

