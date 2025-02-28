
# load the tidyverse libraries 
source('load_tidyverse_packs.r')
library("reshape2")

# read in CSV of the raw dose-response data
df1 <- read.csv('module4_week1_exercise_15cpds_dose-response_data.csv')

# remove the first 5 lines and keep only necessary columns
df2 <- df1 %>% 
    slice( 6:n() ) %>%
    select( PUBCHEM_CID, PUBCHEM_EXT_DATASOURCE_SMILES, contains("Activity.at"))

# can we melt this data into a tidy format like this...
# PUBCHEM_CID, SMILES, CONC, Activity
# 5310645, CC1=CC2=C(C=CC=C2N=N1)N, 0.049, -2.95
# 5310645, CC1=CC2=C(C=CC=C2N=N1)N, 0.111, -3.1
# 5310645, CC1=CC2=C(C=CC=C2N=N1)N, 0.223, -2.75
# 5310645, CC1=CC2=C(C=CC=C2N=N1)N, 0.389, -3.15
# ...

# fix column names
x <- colnames(df2)
x <- gsub( "Activity.at.", "", x )
x <- gsub( ".uM", "", x)
x[2] <- "SMILES"
colnames(df2) <- x

# isolate smiles
df_smiles <- df2 %>%
    select( c("PUBCHEM_CID","SMILES") )

# isolate dose-response (remove SMILES column)
df_dr <- df2 %>%
    select( -c("SMILES") )

# melt df_dr
tidy_df2 <- df_dr %>%
    melt( id.vars="PUBCHEM_CID", variable.name='conc_uM', value.name='activity' )

# now add back in the SMILES with a left join
tidy_df3 <- tidy_df2 %>%
    left_join( df_smiles, by="PUBCHEM_CID")

# fix types
tidy_df3$conc_uM <- as.numeric( as.character( tidy_df3$conc_uM))
tidy_df3$activity <- as.numeric( tidy_df3$activity)

# remove rows missing activity values (NA)
#tidy_df3 <- tidy_df3 %>%
#    filter( !is.na("activity") )
tidy_df4 <- tidy_df3 %>%
    drop_na( activity )

# flip sign on assay responses (negative) so now desired activity is postive-valued
tidy_df4['activity'] <- tidy_df4['activity'] * -1.000

# write to CSV file (without index numbers)
write.csv( tidy_df4, "tidy_15cpds_dose-response_data.csv", row.names=FALSE )
