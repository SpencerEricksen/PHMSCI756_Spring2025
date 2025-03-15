# PHMSCI 756 MINI PROJ ####
# PULL OUT "Fit_Curve_Class", "Potency" AND USE TO FIND 40 HITS

#load in libraries
library("readr")

# import data
AID_1558 <- read_delim("AID_1558_datatable_week3_exercise.csv")
AID_1559 <- read_delim("AID_1559_datatable_week3_exercise.csv")
AID_1694 <- read_delim("AID_1694_datatable_week3_exercise.csv")

# select relevant columns from each data frame
subset_1558 <- AID_1558[6:nrow(AID_1558), c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "Fit_CurveClass", "Potency")]

subset_1559 <- AID_1559[6:nrow(AID_1559), c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "Fit_CurveClass", "Potency")]

subset_1694 <- AID_1694[6:nrow(AID_1694), c("PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "Fit_CurveClass", "Potency")]

# Rename the columns of each data frame, then rename PUBCHEM_SID & SMILES
names(subset_1558) <- paste0(names(subset_1558), "_1558")
names(subset_1559) <- paste0(names(subset_1559), "_1559")
names(subset_1694) <- paste0(names(subset_1694), "_1694")

# Rename the first column of each subset to "PUBCHEM_SID"
colnames(subset_1558)[1] <- "PUBCHEM_SID"
colnames(subset_1559)[1] <- "PUBCHEM_SID"
colnames(subset_1694)[1] <- "PUBCHEM_SID"

# Rename the second column of each subset to "PUBCHEM_CID"
colnames(subset_1558)[2] <- "PUBCHEM_CID"
colnames(subset_1559)[2] <- "PUBCHEM_CID"
colnames(subset_1694)[2] <- "PUBCHEM_CID"

# Rename the third column of each subset to "PUBCHEM_EXT_DATASOURCE_SMILES"
colnames(subset_1558)[3] <- "PUBCHEM_EXT_DATASOURCE_SMILES"
colnames(subset_1559)[3] <- "PUBCHEM_EXT_DATASOURCE_SMILES"
colnames(subset_1694)[3] <- "PUBCHEM_EXT_DATASOURCE_SMILES"

# Put all data frames in a list
data_frames <- list(subset_1558, subset_1559, subset_1694)

# Use Reduce to merge all data frames
merged_data <- Reduce(function(x, y) merge(x, y, by = "PUBCHEM_SID", all = TRUE), data_frames)

# Remove extra SMILES
merged_data$PUBCHEM_EXT_DATASOURCE_SMILES.x <- NULL
merged_data$PUBCHEM_EXT_DATASOURCE_SMILES.y <- NULL

# Remove extra PUBCHEM_CID
merged_data$PUBCHEM_CID.x <- NULL
merged_data$PUBCHEM_CID.y <- NULL

swissADME <- read_delim("swissADME.csv")
swissADME <- swissADME[, c("Lipinski #violations", "PAINS #alerts", "BBB permeant")]

# Attach swissADME to merged data
merged_data <- cbind(merged_data, swissADME)

# Remove row with Lipinski & PAINS violations
merged_data <- subset(merged_data, `Lipinski #violations` == 0 & `PAINS #alerts` == 0)

num_rows <- nrow(subset(merged_data, `Potency_1694` > 5))
num_rows

# Filter rows where Potency_1558 < 80
merged_data <- subset(merged_data, `Potency_1558` < 80)

# Filter rows where Potency_1559 < 80
merged_data <- subset(merged_data, `Potency_1559` < 80)

# Sort by average Fit_CurveClass of primary assays, assign to new column Fit_CurveClass_avg

merged_data$Fit_CurveClass_1558 <- as.numeric(as.character(merged_data$Fit_CurveClass_1558))
merged_data$Fit_CurveClass_1559 <- as.numeric(as.character(merged_data$Fit_CurveClass_1559))
merged_data$Fit_curveclass_avg <- rowMeans(merged_data[, c("Fit_CurveClass_1558", "Fit_CurveClass_1559")])

# Sort by Fit_CurveClass_avg, choose 40 lowest 

merged_data <- merged_data[order(merged_data$Fit_curveclass_avg), ]
merged_data <- merged_data[1:40, ]

# write SMILES to csv
write.csv(merged_data[, c("PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_CID", "PUBCHEM_SID")], file = "smiles40.csv", row.names = FALSE)

