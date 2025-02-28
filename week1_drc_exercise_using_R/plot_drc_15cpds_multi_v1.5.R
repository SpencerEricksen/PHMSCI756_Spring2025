
library('drc')
library('ggplot2')
library('broom')
#source('load_tidyverse_packs.r')
library('tidyr')
library('dplyr')

# load the CSV (tidy format) 
df_drc15 <- read.csv("tidy_15cpds_dose-response_data.csv")


# get max_obs values
max_act <- tapply( df_drc15$activity, df_drc15$PUBCHEM_CID, max )
df_max <- data.frame(max_act)
df_max <- cbind( newColName=rownames(df_max), df_max )
rownames(df_max) <- NULL
colnames(df_max) <- c('PUBCHEM_CID','max_activity')
df_max$PUBCHEM_CID <- as.integer(df_max$PUBCHEM_CID)

# fit 4-pt logistic regression to FP data
# fit 4PL to FP data set (15 different samples)

drc15.m1 <- drm(
    data = df_drc15,
    formula = activity~conc_uM,
    curveid = PUBCHEM_CID, 
    fct=LL.4( names=c("Hill_slope","min","max","EC50" )) 
)

# for 4PL, the parameters b, c, d, e are estimated
# f(x) = c + (d-c)/{1 + exp[ b(log(x) - log(e)) ] }
# c=min
# d=max
# b=Hill slope
# e=potency (EC50, IC50)

# dump fit params to CSVs
df_params <- tidy(drc15.m1)
colnames(df_params) = c("param","PUBCHEM_CID","coeff","stderr","t-val","p-val")
df_params$PUBCHEM_CID <- as.integer(df_params$PUBCHEM_CID)

# add in the max values
df_params <- df_params %>%
    left_join( df_max, by="PUBCHEM_CID" )

# output in Tidy format
write.csv( df_params, "drc15_4PL_fit_params_tidy.csv", row.names=FALSE )

# output in wide format
#df_params_wide <- pivot_wider( df_params, names_from='param', values_from=c('coeff','stderr','t-val','p-val','max_activity') )
df_params_wide <- pivot_wider( df_params, names_from='param', values_from=c('coeff','stderr','t-val','p-val') )
write.csv( df_params_wide, "drc15_4PL_fit_params_wide.csv", row.names=FALSE)


# open file for graphics output
png('drc15_4PL_fit_curves.png', units="in", width=10, height=10, res=600)

# loop through all 15 cpds

# set up panel image
par(mfrow=c(4,4))
cid_list <- unique( df_drc15$PUBCHEM_CID )
for( cid in cid_list ) {
    plot( drc15.m1, level=c(cid), main=cid )
}
# close png file
dev.off()


