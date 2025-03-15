# title: HTS Hit Identification

library(tidyverse)
library(writexl)

#Read data files
ThT_1558_raw <- read.csv("AID_1558_datatable_week3_exercise.csv")
FPmP_1559_raw <- read.csv("AID_1559_datatable_week3_exercise.csv")
FPtotF_1694_raw <- read.csv("AID_1694_datatable_week3_exercise.csv")

raw_files <- list(ThT_1558_raw, FPmP_1559_raw, FPtotF_1694_raw)

# Clean data files ThT_1558
ThT_1558_analysis <- ThT_1558_raw[-c(1:5),]
ThT_1558_analysis <- select(ThT_1558_analysis, 
                        PUBCHEM_CID, PUBCHEM_EXT_DATASOURCE_SMILES, Potency, Efficacy,
                        contains("Fit", F), Max_Response)

####### figure out data types of col
str(ThT_1558_analysis)

####### fix coltypes
ThT_1558_analysis[,1] <- as.character(ThT_1558_analysis[,1])
for (x in 2:ncol(ThT_1558_analysis)) { 
  if (x != 2){
   ThT_1558_analysis[,x] <- as.double(ThT_1558_analysis[,x]) 
  }
  colnames(ThT_1558_analysis)[x] <- paste(colnames(ThT_1558_analysis)[x],"1558",sep="_") 
}

# Clean data files FPmp_1559
FPmP_1559_analysis <- FPmP_1559_raw[-c(1:5),]
FPmP_1559_analysis <- select(FPmP_1559_analysis, 
                            PUBCHEM_CID, PUBCHEM_EXT_DATASOURCE_SMILES, Potency, Efficacy,
                            contains("Fit", F), Max_Response)


####### fix coltypes
FPmP_1559_analysis[,1] <- as.character(FPmP_1559_analysis[,1])
for (x in 2:ncol(FPmP_1559_analysis)) { 
  if (x != 2){
    FPmP_1559_analysis[,x] <- as.double(FPmP_1559_analysis[,x]) 
  }
  colnames(FPmP_1559_analysis)[x] <- paste(colnames(FPmP_1559_analysis)[x],"1559",sep="_") 
}

####### figure out data types of col
str(FPmP_1559_analysis)


# Clean data files FPtotF_1694
FPtotF_1694_analysis <- FPtotF_1694_raw[-c(1:5),]
FPtotF_1694_analysis <- select(FPtotF_1694_analysis, 
                             PUBCHEM_CID, PUBCHEM_EXT_DATASOURCE_SMILES, Potency, Efficacy,
                             contains("Fit", F), Max_Response)


####### fix coltypes
FPtotF_1694_analysis[,1] <- as.character(FPtotF_1694_analysis[,1])
for (x in 2:ncol(FPtotF_1694_analysis)) { 
  if (x != 2){
    FPtotF_1694_analysis[,x] <- as.double(FPtotF_1694_analysis[,x]) 
  }
  colnames(FPtotF_1694_analysis)[x] <- paste(colnames(FPtotF_1694_analysis)[x],"1694",sep="_") 
}

####### figure out data types of col
str(FPtotF_1694_analysis)



#merge
full_data_analysis <- inner_join(ThT_1558_analysis, FPmP_1559_analysis, by = join_by(PUBCHEM_CID), multiple = "first")
full_data_analysis <- inner_join(full_data_analysis, FPtotF_1694_analysis, by = join_by(PUBCHEM_CID), multiple = "first")

#write_xlsx(full_data_analysis, "full_HTS_hits_validation.xlsx")


####################################################################################


#SwissADME Data
swissADME_raw <- read.csv("swissadme HTS Validation.csv")
str(swissADME_raw)
swissADME_raw[,1] <- as.character(swissADME_raw[,1])

#merge
full_data_analysis <- inner_join(full_data_analysis, swissADME_raw, by = join_by(PUBCHEM_CID == Molecule), multiple = "first")
str(full_data_analysis)

#filter substances
##### PAINS / Lipinski = no violations
active_samples <- filter(full_data_analysis, PAINS_alerts == 0, Lipinski_violations==0)
active_samples <- active_samples[(abs(active_samples$Fit_CurveClass_1558) <= 2),]
active_samples <- active_samples[(abs(active_samples$Fit_CurveClass_1559) <= 2),]
active_samples <- arrange(active_samples, Potency_1558)

#write_xlsx(active_samples, "active_samples_HTS_hits_validation.xlsx")

active_samples_clustering_raw <- read.csv("active_samples_HTS_clustering.csv")
str(active_samples_clustering_raw)
active_samples_clustering_raw[,1] <- as.character(active_samples_clustering_raw[,1])

active_clustered_samples <- right_join(active_samples, active_samples_clustering_raw, by = join_by(PUBCHEM_CID == ids), multiple = "first")
active_clustered_samples <- active_clustered_samples %>% 
  group_by(CLID_0.4) %>%
  arrange(Potency_1558, .by_group = TRUE)
write_xlsx(active_clustered_samples, "active_sorted.xlsx")
