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

names(all.drugs) = make.names(names(all.drugs))

#### MACCs fingerprints, these are more interpretable I think, looking into what
# 'extended' meant, there's some hashing and hopefully this way we can still
# get some performance boosts while actually giving a meaning

fp_to_df <- function(fp_list) {
  if (length(fp_list) == 0) return(tibble())
  mat <- fingerprint::fp.to.matrix(fp_list)
  mat <- as.matrix(mat)
  colnames(mat) <- paste0("fp_", seq_len(ncol(mat)))
  as.data.frame(mat)
}

fingerprints = map(all.drugs$parsed.SMILES, ~ get.fingerprint(.x, type = "maccs"))

all.drugs = cbind(all.drugs, fp_to_df(fingerprints))


#### MOA by Indication ####

all.drugs %>%
  separate_rows(MOA, sep = ",\\s*") %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  separate_rows(`Disease.Area`, sep = ",\\s*") %>%
  filter(`Disease.Area` == "neurology/psychiatry",
         Indication == c("depression", "schizophrenia",
                         "Parkinson's Disease")) %>%
  group_by(MOA) %>%
  filter(n() > 1) %>%
  ggplot(aes(y = MOA)) +
  geom_bar(stat = "count", fill = "maroon") +
  facet_wrap(vars(Indication), scales = "free_x") +
  ylab("Mechanism of Action (MOA)") +
  xlab("# of Occurences") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", # get rid of ggplot grey face headers
                                        color = "white"),
        strip.text = element_text(size = 14)) +
  guides(fill = "none") +
  scale_x_continuous(
    breaks = function(x) {
      rng <- floor(min(x)):ceiling(max(x))
      even <- rng[rng %% 2 == 0]
      even
    }
  )


#### Target Proteins by Indication

all.drugs %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  separate_rows(Target, sep = ",\\s*") %>%
  mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target)) %>%
  group_by(TargetPrefix, Indication) %>%
  # filter(n() >= 5) %>%
  ungroup() %>%
  # filter(Indication %in% c("depression", "Parkinson's Disease", "schizophrenia", "pain relief", "seizures")) %>%
  filter(Indication %in% c("depression", "Parkinson's Disease",
                           "schizophrenia", "pain relief", "seizures"),
         TargetPrefix %in% c(
           "SLC",
           "SCN",
           "PTGS",
           "HTR",
           "HRH",
           # "GRIN",
           "GABRA",
           "DRD",
           "CHRM",
           "CACNG",
           "CACNA",
           "ADRB",
           "ADRA"
         )) %>%
  mutate(
    Indication = factor(
      Indication,
      levels = c(
        "depression",
        "Parkinson's Disease",
        "schizophrenia",
        "seizures",
        "pain relief"
      )
    )
  ) %>%
  ggplot(aes(y = TargetPrefix)) +
  geom_bar(stat = "count", fill = "navy") +
  facet_wrap(vars(Indication), scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", # get rid of ggplot grey face headers
                                    color = "white"),
    strip.text = element_text(size = 14),
    # panel.grid.major.x = element_line(color = "white"),
    # panel.grid.minor = element_line(color = "white")
  ) +
  labs(
    x = "Count",
    y = "Target Protein Prefix"
  )

all.drugs %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  separate_rows(Target, sep = ",\\s*") %>%
  filter(Indication %in% c("depression", "Parkinson's Disease",
                           "schizophrenia", "pain relief", "seizures")) %>%
  group_by(Target, Indication) %>%
  # filter(n() >= 2) %>%
  ungroup() %>%
  mutate(
    Indication = factor(
      Indication,
      levels = c(
        "depression",
        "Parkinson's Disease",
        "schizophrenia",
        "seizures",
        "pain relief"
      )
    )
  ) %>%
  ggplot(aes(y = Target)) +
  geom_bar(stat = "count", fill = "navy") +
  facet_wrap(vars(Indication), scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", # get rid of ggplot grey face headers
                                    color = "white"),
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 8, face = "italic")
    # panel.grid.major.x = element_line(color = "white"),
    # panel.grid.minor = element_line(color = "white")
  ) +
  labs(
    x = "Count",
    y = "Target Protein Prefix"
  )

# all.drugs %>%
#   separate_rows(Indication, sep = ",\\s*") %>%
#   separate_rows(Target, sep = ",\\s*") %>%
#   mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target)) %>%
#   group_by(TargetPrefix, Indication) %>%
#   filter(n() > 5) %>%
#   ungroup() %>%
#   filter(Indication %in% c("depression", "Parkinson's Disease", "schizophrenia")) %>%
#   ggplot(aes(y = TargetPrefix)) +
#   geom_bar(
#     aes(x = after_stat(count / sum(count))),
#     stat = "count"
#   ) +
#   facet_wrap(vars(Indication), scales = "free_x") +
#   labs(x = "Percentage", y = "Target Prefix")




# all.drugs %>%
#   select(-first.SMILES, -parsed.SMILES) %>%
#   write.csv(file = "allDrugs1s0s.csv")



#### Number of unique target proteins ####

unique(all.drugs %>%
         separate_rows(Target, sep = ",\\s*") %>%
         pull(Target))

# 2019 unique target proteins
length(unique(all.drugs %>%
                separate_rows(Target, sep = ",\\s*") %>%
                pull(Target)))

# 1004 unique protein groups based on character prefix
length(unique(all.drugs %>%
                separate_rows(Target, sep = ",\\s*") %>%
                mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target)) %>%
                pull(TargetPrefix))
       )

# 142 unique target proteins for our indications
length(unique(all.drugs %>%
                separate_rows(Indication, sep = ",\\s*") %>%
                filter(Indication %in% c("depression", "Parkinson's Disease", "schizophrenia")) %>%
                separate_rows(Target, sep = ",\\s*") %>%
                pull(Target)))

# 69 unique protein groups based on character prefix for our indication
length(unique(all.drugs %>%
                separate_rows(Indication, sep = ",\\s*") %>%
                filter(Indication %in% c("depression", "Parkinson's Disease", "schizophrenia")) %>%
                separate_rows(Target, sep = ",\\s*") %>%
                mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target)) %>%
                pull(TargetPrefix))
)

# 18 unique SLC proteins

unique(all.drugs %>%
  separate_rows(Disease.Area, sep = ",\\s*") %>%
  filter(str_detect(Disease.Area, "neurology/psychiatry")) %>%
  separate_rows(Target, sep = ",\\s*") %>%
  filter(str_detect(Target, "SLC")) %>%
  pull(Target))

length(
  unique(all.drugs %>%
           separate_rows(Disease.Area, sep = ",\\s*") %>%
           filter(str_detect(Disease.Area, "neurology/psychiatry")) %>%
           separate_rows(Target, sep = ",\\s*") %>%
           filter(str_detect(Target, "SLC")) %>%
           pull(Target))
)

# 14 unique HTR proteins

unique(all.drugs %>%
         separate_rows(Disease.Area, sep = ",\\s*") %>%
         filter(str_detect(Disease.Area, "neurology/psychiatry")) %>%
         separate_rows(Target, sep = ",\\s*") %>%
         filter(str_detect(Target, "HTR")) %>%
         pull(Target))

length(
  unique(all.drugs %>%
           separate_rows(Disease.Area, sep = ",\\s*") %>%
           filter(str_detect(Disease.Area, "neurology/psychiatry")) %>%
           separate_rows(Target, sep = ",\\s*") %>%
           filter(str_detect(Target, "HTR")) %>%
           pull(Target))
)

# all.drugs %>%
#   separate_rows(Indication, sep = ",\\s*") %>%
#   separate_rows(Target, sep = ",\\s*") %>%
#   mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target)) %>%
#   group_by(TargetPrefix, Indication) %>%
#   filter(Indication %in% c("depression", "Parkinson's Disease", "schizophrenia", "pain relief", "seizures")) %>%
  

# 636 unique indications
length(unique(all.drugs %>%
                separate_rows(Indication, sep = ",\\s*") %>%
                pull(Indication)))

all.drugs %>%
  filter(str_detect(`Disease.Area`, "neurology/psychiatry"),
         Phase == "Launched") %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  group_by(Indication) %>%
  summarise(n = n()) %>%
  arrange(desc(n))



# Does target protein count look related to # of indications treated?

all.drugs %>%
  filter(Indication != "<NA>") %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  separate_rows(Target, sep = ",\\s*") %>%
  group_by(Name) %>%
  summarize(
    num_indications = n_distinct(Indication),
    num_targets = n_distinct(Target)
  ) %>%
  ggplot(aes(x = num_targets, y = num_indications)) +
  geom_point() +
  geom_bin_2d(binwidth = c(1,1)) +
  geom_smooth(method = "lm") +
  labs(
    title = "Relationship Between Number of Targets and Indications",
    x = "# of Targets",
    y = "# of Indications"
  ) +
  theme_minimal()



