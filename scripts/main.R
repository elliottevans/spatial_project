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

makepic('mich_grad_rts',width=4,height=4)
plot(michigan_sp,col="white",axes=F)

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
box(col = 'black')
dev.off()
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


############################################
#Create a map of % black per county.
###########################################
makepic('percent_black',width=4,height=4)
plot(michigan_sp,col="white",axes=F)

plotvar <- mich_grad_dat$percent_black
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(mich_grad_dat$percent_black,  probs = c(0,20,40,60,80,100)/100))
class
### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)

plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[0.0%, 0.2)","[0.2, 2.8)","[2.8, 4.4)","[4.4, 9.1)","[9.1, 12.5%]")
legend("bottomleft",legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
title(main="Michigan High Schools\n Percent African-American")
box(col = 'black')
dev.off()
############################################
#Create a map of % black per county.
###########################################


############################################
#Create a map of % economically disadvantaged per county.
###########################################
makepic('percent_econ_disadvan',width=4,height=4)
plot(michigan_sp,col="white",axes=F)

plotvar <- mich_grad_dat$percent_econ_disadvan
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(mich_grad_dat$percent_econ_disadvan,  probs = c(0,20,40,60,80,100)/100))
class
### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)

plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[19.8%, 34.6)","[34.6, 41.4)","[41.4, 46.2)","[46.2, 50.9)","[50.9, 100%]")
legend("bottomleft",legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
title(main="Michigan High Schools\n Percent Economically Disadvantaged")
box(col = 'black')
dev.off()
############################################
#Create a map of % economically disadvantaged per county.
###########################################


##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with an improper CAR prior
##############################################
### Now we fill the adjacency matrix W using the following trick 
W <- matrix(0,nrow(mich_grad_dat),nrow(mich_grad_dat))

rep.mich <- rep(1:nrow(mich_grad_dat),num.mich)

for(i in 1:nrow(mich_grad_dat)){
	W[i,adj.mich[rep.mich==i]] <- rep(1,num.mich[i])	
}

formula <- mich_grad_dat$graduation_rt~per_econ_dis_cent+per_black_cent
model.car <- S.CARleroux(formula=formula,
                         W=W,
                         family="gaussian",
                         fix.rho=TRUE,
                         rho=1,
                         burnin=70000, 
                         n.sample=200000,
                         thin=1, 
                         prior.mean.beta=NULL, 
                         prior.var.beta=NULL, 
                         prior.nu2=NULL,
                         prior.tau2=NULL, 
                         verbose=TRUE)

samples.eta <- model.car$samples$phi
model.car$summary.results
2*pnorm(abs(model.car$summary.results[,7]), mean = 0, sd = 1, lower.tail = FALSE)

stargazer(round(model.car$summary.results[,-c(4,5)],3))

#Trace plot for beta hats on intercept, per_econ_dis_cent, and per_black_cent
savepdf('trace_plots',width=17,height=16)
plot(model.car$samples$beta,main="Beta Coefficients")
dev.off()

#Trace plots for nu2 and tau2
plot(model.car$samples$nu2,main='Nu2')
plot(model.car$samples$tau2,main='Tau2')

makepic('med_spat_rand',width=4,height=4)
#Posterior median of spatial effects
post.median.eta <- as.numeric(apply(samples.eta,2,median))
plotvar <- post.median.eta
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(post.median.eta,  probs = c(0,20,40,60,80,100)/100))
class
colcode <- findColours(class,plotclr)

plot(michigan_sp,border="black",axes=F)
title(main="Median of Spatial Random Effects")
plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[-0.02, -0.003)","[-0.003, 0.0007)","[0.0007, 0.004)","[0.004, 0.01)","[0.01, 2.07]")
legend("bottomleft",legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
box(col = 'black')
dev.off()

###Posterior SD of spatial effects
makepic('sd_spat_rand',width=4,height=4)
post.sd.eta <- as.numeric(apply(samples.eta,2,sd))
plotvar <- post.sd.eta
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(post.sd.eta,  probs = c(0,20,40,60,80,100)/100))
class
colcode <- findColours(class,plotclr)

plot(michigan_sp,border="black",axes=F)
title(main="SD of Spatial Random Effects")
plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[0.32, 0.37)","[0.37, 0.43)","[0.43, 0.54)","[0.54, 1.00)","[1.00, 1.63]")
legend("bottomleft",legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
box(col = 'black')
dev.off()

#Get credible intervals for each state
lower<-c()
upper<-c()
for(i in 1:dim(samples.eta)[2]){
  lower<-append(lower,quantile(samples.eta[,i],  probs = c(.025,.975))[1])
  upper<-append(upper,quantile(samples.eta[,i],  probs = c(.025,.975))[2])
}

cbind(lower,upper)
#all intervals contain 0,
# no significantly positive/negative spatial random effects

##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with an improper CAR prior
##############################################



##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with a PROPER CAR model
##############################################

mich.w.proper.car <- nb2listw(mich.nb,style="B")

model.proper.car <- spautolm(formula=formula,
                             listw=mich.w.proper.car,
                             family="CAR")
summary(model.proper.car)

##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with a PROPER CAR model
##############################################



##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with SAR model
##############################################

model.sar <- spautolm(formula=formula,listw=mich.listw,family="SAR")
summary(model.sar)

##############################################
#Here, we fit the Bayesian hierarchical spatial
#linear model with SAR model
##############################################




###############################################
#Disease Modeling Approach
###############################################

### The function S.CARbym fits spatial model to the Poisson data with two sets of random effects, the spatial
### and non spatial ones.
### Specifically, the model that will be fit is a Poisson model for the observed counts
### with mean (parameter) equal to E*relative risk where E represents the expected counts.
### In the S.CARbym function, we specify a formula that describes how the log of the
### mean of the Poisson distribution for the observed counts depends on covariates.
### Taking the log of E*relative risk, we have log(E)+log(relative risk).
### So the formula command expresses how the log(relative risk) depends on covariates and adds as offset
### the log(E), or the log of the expected counts.

### Here, since we assume that the log relative risk depends on the percentage of population 
### involved in Agriculture, Fisheries and Forestries, we have:

num_graduating<-as.integer(mich_grad_dat$num_graduating)
num_cohort<-mich_grad_dat$cohort_count

formula2 <- num_graduating~per_econ_dis_cent + per_black_cent + offset(log(num_cohort))

model.poi <- S.CARbym(formula=formula2, family="poisson", W=W, burnin=30000, n.sample=170000)

model.poi$summary.results

savepdf('trace_plots_poi',width=16,height=15)
#Trace plots
plot(model.poi$samples$beta)
dev.off()
###############################################
#Disease Modeling Approach
###############################################



###############################################
#Plot the propensity to graduate
###############################################
#Make a cloropleth of the estimated propensity to graduate

### Estimated relative risks / PROPENSITY TO GRADUATE
### The relative risks are given by exp(X*beta+random effects) and are stored in the fitted values
 
###  This gets the component of the log relative risk that is given by X*beta.
###  The result of this is a matrix with as many rows as MCMC samples we have for the beta coefficients and as many columns
###  as number of areal units
RR.covariate <- model.poi$samples$beta%*%t(matrix(cbind(rep(1,nrow(mich_grad_dat)),per_econ_dis_cent,per_black_cent),nrow(mich_grad_dat),3))

###  This gets the component of the log relative risk that is given by the random effects
###  The result of this is a matrix with as many rows as MCMC samples we have for the spatial random effects (phi) and as many columns
###  as number of areal units

RR.random <- model.poi$samples$psi
#+model.scotland$residuals[,3]
RR.all <- exp(RR.covariate+RR.random)

post.median.RR <- apply(RR.all,2,median)
plotvar <- post.median.RR
nclr<-5

plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(post.median.RR,  probs = c(0,20,40,60,80,100)/100))
class

colcode <- findColours(class,plotclr)

makepic('propensity',width=4,height=4)
plot(michigan_sp,border="black",axes=F)
title(main="Estimated Propensity to\n Graduate High School")
plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[0.74, 0.83)","[0.83, 0.85)","[0.85, 0.87)","[0.87, 0.88)","[0.88, 0.91]")
legend('bottomleft',legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
box(col = 'black')
dev.off()
###############################################
#Plot the propensity to graduate
###############################################




###############################################
#The posterior median of the spatial random effects
###############################################
post.median.psi <- as.numeric(apply(model.poi$samples$psi,2,median))
plotvar <- post.median.psi
nclr<-5

plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=quantile(post.median.psi,  probs = c(0,20,40,60,80,100)/100))
class

colcode <- findColours(class,plotclr)

makepic('post_med_poi',width=4,height=4)
plot(michigan_sp,border="black",axes=F)
title(main="Posterior Median of\n Spatial Random Effects")
plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("[-0.70, -0.02)","[-0.02, -0.001)","[-0.001, 0.01)","[0.01, 0.02)","[0.02, 0.05]")
legend('bottomleft',legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
box(col = 'black')
dev.off()

#the lower bound of their 95% CI
post.lower.psi <- as.numeric(apply(model.poi$samples$psi,2,quantile,probs=c(2.5)/100))


#the upper bound of their 95% CI
post.upper.psi <- as.numeric(apply(model.poi$samples$psi,2,quantile,probs=c(97.5)/100))

# counties with significantly positive spatial random effects
nrow(mich_grad_dat[(post.lower.psi>0),])

# counties with significantly negative spatial random effects
nrow(mich_grad_dat[post.upper.psi<0,])

plotvar<-c()
for(i in 1:nrow(mich_grad_dat)){
  #sig. positive spatial random effects
  if(post.lower.psi[i]>0){
    plotvar<-append(plotvar,post.lower.psi[i])
  }
  #sig. negative spatial random effects
  else if(post.upper.psi[i]<0){
    plotvar<-append(plotvar,post.lower.psi[i])
  }
  else{
    plotvar<-append(plotvar,0)
    #print(as.character(election_dat$State[i]))

  }
}

nclr<-2

plotclr <- c("red","grey80")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(plotvar), 0, max(plotvar)), intervalClosure='left')


colcode <- findColours(class,plotclr)

makepic('sig_effects',width=4.2,height=4)
plot(michigan_sp,border="black",axes=F)
title(main="Significant Spatial Random Effects")
plot(michigan_sp,col=colcode,add=T)

leg.txt<-c("Sig. Negative","Not Significant")
legend('bottomleft',legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
box(col = 'black')
dev.off()
###############################################
#The posterior median of the spatial random effects
###############################################

