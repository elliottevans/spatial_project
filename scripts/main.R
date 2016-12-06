######################################################################
#This is the main script to perform the spatial analysis. A more detailed
#account of the mich_grad_dat data frame is included in the prep.R file.
######################################################################

setwd("~/spatial_project")
#Load the libraries, create main dataframe mich_grad_dat
source('scripts/prep.R')

############################################
#Create a map of graduation rts per county.
###########################################
michigan_map <- map("county", 'michigan',plot=FALSE, fill = TRUE)
#michigan_map$names<-michigan_map$names[-42]

IDs<-michigan_map$names

michigan_sp <- map2SpatialPolygons(michigan_map, IDs = IDs,
                            proj4string = CRS("+proj=longlat +ellps=WGS84"))

plot(michigan_sp,col="white",axes=T)

plotvar <- mich_grad_dat$graduation_rt
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(mich_grad_dat$graduation_rt,  probs = c(0,20,40,60,80,100)/100))
class
### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)

plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[60%, 81)","[81, 84)","[84, 87)","[87, 89)","[89, 98%]")
legend("bottomleft",legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
title(main="Michigan High School\n Graduation Rates")
############################################
#Create a map of graduation rts per county.
###########################################


############################################
#Create the adjacency matrix, calculate Moran's I
###########################################
mich.nb <- poly2nb(michigan_sp)
mich.coord <- coordinates(michigan_sp)

## Here we get the adjacency weights, together with the information on the number of weights for each county
mich.weights <- nb2WB(mich.nb)
adj.mich <- mich.weights$adj
weights.mich <- mich.weights$weights
num.mich <- mich.weights$num
mich.listw <- nb2listw(mich.nb)

## Here we are computing Moran's I to see if the data is spatially correlated
#P-VALUE IS SMALL, MORAN's I IS POSITIVE SO THERE IS EVIDENCE OF STRONG SPATIAL ASSOCIATION
moran.test(mich_grad_dat$graduation_rt, mich.listw, rank=TRUE)
############################################
#Create the adjacency matrix, calculate Moran's I
###########################################


############################################
#Fit a model to the large scale trend for graduation_rt,
#then use minimizing AIC as criteria for determining good
#predictors of graduating 
###########################################

fit<-lm(graduation_rt~percent_male+percent_female+percent_econ_disadvan+percent_hispanic+percent_white+percent_black+percent_asian,data=mich_grad_dat)

#This gives us percent_econ_disadvan and percent_black as the best covariates
step(fit)

#center the variables so that the intercept is the mean graduation rt
per_econ_dis_cent<-mich_grad_dat$percent_econ_disadvan-mean(mich_grad_dat$percent_econ_disadvan)
per_black_cent<-mich_grad_dat$percent_black-mean(mich_grad_dat$percent_black)

fit<-lm(mich_grad_dat$graduation_rt~per_econ_dis_cent+per_black_cent)

#Moran's I on residuals show's there is still spatial variation that is unexplained by the covariates
moran.test(fit$residuals, mich.listw, rank=TRUE)

############################################
#Fit a model to the large scale trend for graduation_rt,
#then use minimizing AIC as criteria for determining good
#predictors of graduating 
###########################################





