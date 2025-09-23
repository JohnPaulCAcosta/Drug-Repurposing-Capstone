library(tidyverse)
library(readxl)

########################################################################
## Code to create counts the for dataset on the number of occurrences ##
## for each Disease Area, Target Protein, and Indication + Plotting   ##
########################################################################

all.drugs = read_xls("dataCapstone.xls") %>%
  filter(!is.na(SMILES))

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

nP = launched %>% # 368 launched
  filter(str_detect(`Disease Area`, "neurology/psychiatry"),
         !is.na(SMILES))

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

nP = launched %>% # 368 launched
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
), aes(x = Counts, y = Names, fill = Names)) +
  geom_bar(stat = "identity") +
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
    mutate(!!prefix := if_else(
      any(str_starts(
        str_split(Target, ", ")[[1]],
        as.character(prefix)
      )),
      1,
      0
    )) %>%
    mutate(!!paste0(prefix, ".count") := sum(str_starts(
      str_split(Target, ", ")[[1]],
      as.character(prefix)
    ))
    ) %>%
    ungroup()
  
}

#### Create columns for the presence of every neuro/psych indication ####

neuro.psych = all.drugs %>%
  filter(str_detect(`Disease Area`, "neurology/psychiatry"), 
         Phase == "Launched")

neuro.psych$IndicationVector = strsplit(neuro.psych$Indication, ", ")
counts.indication = list()
counter.indication = function(vec) {
  
  for (i in vec) {
    
    if (is.null(counts.indication[[i]])) {
      counts.indication[[i]] <<- 0
    }
    
    counts.indication[[i]] <<- counts.indication[[i]] + 1
    
  }
  
}
for (vecT in neuro.psych$IndicationVector) {
  counter.indication(vecT)
}

unique.ind = sort(names(counts.indication)) # 187 unique indications!

for (ind in unique.ind) {
  
  neuro.psych = neuro.psych %>%
    rowwise() %>%
    mutate(!!ind := if_else(
      any(str_starts(
        str_split(Indication, ", ")[[1]],
        as.character(ind)
      )),
      1,
      0
    )) %>%
    ungroup()
  
}

#### Analyze distribution of target proteins across various neuro/psych indications ####

counts.indication

# The most common indications are depression and schizophrenia
counts.indication[which(counts.indication == max(as.numeric(counts.indication)))]

counts.indication[which(counts.indication > 10)]

# Kind of hard to plot since there's a lot of columns of 1s and 0, so we can try the
# tidyverse "unnesting" route for more complex relationships

# But this route certainly made basic counting easy! Maybe for model fitting too?

dim(neuro.psych)
colnames(neuro.psych)

