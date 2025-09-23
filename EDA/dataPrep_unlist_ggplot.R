library(tidyverse)
library(readxl)
library(rcdk)

########################################################################
## Code to "unlist" Targets & Indication columns, allowing for easier ##
## usage of ggplot                                                    ##
########################################################################

all.drugs = read_xls("dataCapstone.xls") %>%
  filter(!is.na(SMILES))

# From previous EDA in 1s & 0s / counting-based R code

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

#### "Unlist" the Targets, Indications, & Disease Area ####

unlisted.drugs = all.drugs %>%
  separate_rows(Target, sep = ",\\s*") %>%
  separate_rows(Indication, sep = ",\\s*") %>%
  separate_rows(`Disease Area`, sep = ",\\s*") %>%
  mutate(TargetPrefix = sub("^([A-Za-z]+).*", "\\1", Target))

#### ggplotting, must specify two levels of Target, Indications, ####
#### or Disease Area or it may break (I'm not sure)              ####

## Depression, Schizophrenia, Parkinson's, & pain relief
## are the most prevalent neuro/psych Indications in our dataset, 
## so let's look at the effect of the target protein prefixes
## and clustering vs. no clustering

# No clustering at all
unlisted.drugs %>%
  filter(`Disease Area` == "neurology/psychiatry",
         Indication == c("depression", "schizophrenia",
                         "Parkinson's Disease", "pain relief")) %>%
  ggplot(aes(y = Target)) +
  geom_bar(stat = "count") +
  facet_wrap(vars(Indication)) +
  ylab("Target Protein") +
  xlab("# of Occurences") +
  theme_bw()

# Generic first-letters-only clustering
unlisted.drugs %>%
  filter(`Disease Area` == "neurology/psychiatry",
         Indication == c("depression", "schizophrenia",
                         "Parkinson's Disease", "pain relief")) %>%
  ggplot(aes(y = TargetPrefix)) +
  geom_bar(stat = "count",
           fill = "black") +
  facet_wrap(vars(Indication), scales = "free_x") +
  ylab("Target Protein") +
  xlab("# of Occurences") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", # get rid of ggplot grey face headers
                                        color = "white"),
        strip.text = element_text(size = 14))

# Our defined version of clustering - didn't work

# unlisted.drugs %>%
#   filter(`Disease Area` == "neurology/psychiatry",
#          Indication == c("depression", "schizophrenia",
#                          "Parkinson's Disease", "pain relief")) %>%
#   mutate(
#     TargetPrefix = if_else(
#       str_sub(TargetPrefix, 1, 3) %in% popular.prefixes,  # first 3 letters
#       str_sub(TargetPrefix, 1, 3),                        # replace with popular prefix
#       TargetPrefix                                         # else keep original
#     ) %>%
#   ggplot(aes(y = TargetPrefix)) +
#   geom_bar(stat = "count") +
#   facet_wrap(vars(Indication), scales = "free_x") +
#   ylab("Target Protein") +
#   xlab("# of Occurences") +
#   theme_bw()

