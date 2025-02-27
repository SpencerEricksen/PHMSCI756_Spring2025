
library('drc')
library('ggplot2')

# DRC package includes many sample dose-response data 
# sets. One such set is called "spinach"

# or just load the CSV I created for it
# write.csv( spinach, "spinach.csv", row.names=FALSE)
 
df_spinach <- read.csv("spinach.csv")

# fit 4-pt logistic regression to ryegrass data (triplicate)
# fit 4PL to spinach data set (5 different samples)

spinach.m1 <- drm(SLOPE~DOSE, CURVE, data=df_spinach, fct=LL.4() )
#spinach.m1 <- drm(SLOPE~DOSE, CURVE, data=spinach, fct=LL.4() )
fit_params <- summary( spinach.m1 )
df_params <- data.frame( fit_params[3] )
write.csv( df_params, "spinach_drc_4PL_fit_params.csv" )

#sink(file='spinach_drc_4PL_fit_params.txt')
#summary(spinach.m1)
#sink(file=NULL)

# open file for graphics output
png('spinach_drc_4PL_fit_curves.png', units="in", width=5, height=8, res=600)

# set up panel image
par(mfrow=c(3,2))
# plot each sample
plot(spinach.m1, level=c(1), main="c1" )
#text( x=10, y=1.5, 
plot(spinach.m1, level=c(2), main="c2" )
plot(spinach.m1, level=c(3), main="c3" )
plot(spinach.m1, level=c(4), main="c4" )
plot(spinach.m1, level=c(5), main="c5" )

# close png file
dev.off()

