library(tidyverse)
library(readxl)

#### JOURNAL ARTICLE TITLE SCREENING ####

# load search query results
papers <- as.data.frame(SQuery_results_feb_12 <- read_excel("Input/SQuery_results_feb _12.xls"))

# Add study ID
study_id <- seq(1:dim(papers)[1])
papers$study_id <- study_id

#designate search terms 
search_terms <- c('BIO', 'ECO', 'CONSERV', 'AQUAT', 'MARINE', 'NATUR', 'FISH', 'ANIM', 'ZOOL')

#filter articles by journal | filtered from 670 to 422, 248 articles removed via journal name   
jrnl_filt_ppers <- as.data.frame(filter(papers, grepl(paste(search_terms, collapse = '|'), `Source Title`, ignore.case = T)))

#Save papers
write.csv(jrnl_filt_ppers, file = "Input/SQFilteredGrepl.csv")

#### TITLE SCREENING ####
#Load manually title screened list of papers
title <- read_csv("SysLitRevSheets/SQFilteredTitle.csv")

#Add study_id
title$study_id <- jrnl_filt_ppers$study_id

#Change NAs in manually screened journal article column to 0
title$removed_journ_title_man <-replace(title$removed_journ_title_man, is.na(title$removed_journ_title_man), 0)

#Filter by manually screened columns | 422 to 269, 153 articles removed via title screening
filtered_title <- title %>%
  filter(removed_title==0) %>%
  filter(removed_journ_title_man==0) %>% 
  dplyr::select(!c(removed_title, removed_journ_title_man))

# Save list of papers filtered and manually screened by article title and journal article 
write.csv(filtered_title, file = "Input/SQFilteredByTitle.csv")

#### ABSTRACT SCREENING ####
#Load manually screened abstracts
abstract <- read_csv("SysLitRevSheets/SQFilteredAbstract.csv")

#add study_id
abstract$study_id <- filtered_title$study_id

#Filter by manually screened columns | 269 to 129, 140 articles removed via abstract screening
filtered_abstract <- abstract %>%
  filter(screened_abstract==0) %>%
dplyr::select(!screened_abstract)

#Save list of papers for potential effect sizes 
write.csv(filtered_abstract, file = "Input/SQPotentialEffectSizePapers.CSV")

