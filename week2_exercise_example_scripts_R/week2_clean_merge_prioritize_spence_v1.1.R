# title: exercise 2 (cleaning, merging, sort/filter datasets)
# date: 2024-03-11
# author: Spencer Ericksen
# email: ssericksen@wisc.edu

# set working directory
setwd("/home/ssericksen/disk2/Dropbox/DDC/ProjectsTemp/PHMSCI-756-2024SP/week2/week2_exercise_example_scripts_R")
# system('pwd')
library(readr)
library(dplyr)

# import data, 1460=primary, 1468=secondary, 1463=counter
AID_1460 <- read_delim("./data/AID_1460_datatable_week2_exercise.csv")
AID_1468 <- read_delim("./data/AID_1468_datatable_week2_exercise.csv")
AID_1463 <- read_delim("./data/AID_1463_datatable_week2_exercise.csv")


# get cols and rows
test_1460 <- AID_1460 %>%
  select(PUBCHEM_CID, Max_Response, Fit_CurveClass, Potency) %>%
  slice(6:n())

test_1468 <- AID_1468 %>%
  select(PUBCHEM_CID, Max_Response, Fit_CurveClass, Potency) %>%
  slice(6:n())

test_1463 <- AID_1463 %>%
  select(PUBCHEM_CID, Max_Response, Fit_CurveClass, Potency) %>%
  slice(6:n())

# de-duplicate PUBCHEM_CIDs (keep first?)
test_1460 <- test_1460[!duplicated(test_1460$PUBCHEM_CID),]
test_1468 <- test_1468[!duplicated(test_1468$PUBCHEM_CID),]
test_1463 <- test_1463[!duplicated(test_1463$PUBCHEM_CID),]

# setting data types because filter below wasn't working as expected
class(test_1460$PUBCHEM_CID) = "integer"
class(test_1460$Fit_CurveClass) = "double"
class(test_1460$Potency) = "double"
class(test_1460$Max_Response) = "double"

class(test_1468$PUBCHEM_CID) = "integer"
class(test_1468$Fit_CurveClass) = "double"
class(test_1468$Potency) = "double"
class(test_1468$Max_Response) = "double"

class(test_1463$PUBCHEM_CID) = "integer"
class(test_1463$Fit_CurveClass) = "double"
class(test_1463$Potency) = "double"
class(test_1463$Max_Response) = "double"

# remove rows with missing vals
test_1460 <- drop_na( test_1460 )
test_1468 <- drop_na( test_1468 )
#test_1463 <- drop_na( test_1463 )

# apply filters, get hits in 1460 (primary)
test_1460 <- test_1460 %>%
  filter( Fit_CurveClass >= -2.2 ) %>%
  filter( Fit_CurveClass <= 1.0  ) %>%
  filter( Max_Response <= -40.0 ) %>%
  filter( Potency <= 5.0 )

# apply filters, get hits in 1468 (secondary)
test_1468 <- test_1468 %>%
  filter( Fit_CurveClass >= -2.2 ) %>%
  filter( Fit_CurveClass <= 1.0  ) %>%
  filter( Max_Response <= -25.0 ) %>%
  filter( Potency <= 50.0 )

# apply filters, get negatives in 1463 (counter)
test_1463 <- test_1463 %>%
  # remove fluoresence cpds (<= -40.0)
  filter( Max_Response >= -25.0 ) %>%
  # remove quenchers (>= 40.0)
  filter( Max_Response <= 25.0 )

# Rename the columns of each data frame, then rename PUBCHEM_SID
names(test_1460) <- paste0(names(test_1460), "_1460")
names(test_1468) <- paste0(names(test_1468), "_1468")
names(test_1463) <- paste0(names(test_1463), "_1463")

# Rename the first column of each subset to "PUBCHEM_SID"
colnames(test_1460)[1] <- "PUBCHEM_CID"
colnames(test_1468)[1] <- "PUBCHEM_CID"
colnames(test_1463)[1] <- "PUBCHEM_CID"

# merge hits from 1460 and 1468
inner12 <- inner_join( test_1460, test_1468, by='PUBCHEM_CID')
# merge hits with negatives from 1463 (counterscreen)
inner123 <- inner_join( inner12, test_1463, by='PUBCHEM_CID')

# add new columns to use for prioritization of hits
inner123 <- inner123 %>% 
  mutate( mean_curve_class = (Fit_CurveClass_1460 + Fit_CurveClass_1468)/2.0 ) %>%
  mutate( mean_max_response = (Max_Response_1460 + Max_Response_1468) / 2.0 ) 

# sort by mean_curve_class and then by mean_max_response (AID1460 + AID1468)
inner123_rank <- inner123[ with( inner123, order(desc(mean_curve_class), mean_max_response)), ]

#save as excel file
#write_xlsx(inner123, "active_HTS_hits.xlsx")
write_csv(inner123_rank, "active_HTS_hits.csv")