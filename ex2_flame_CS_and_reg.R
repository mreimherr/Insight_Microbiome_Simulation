# In this example we carry out FLAME on abundances 
# as well as a penalized regression for the FB ratio and diet


# To install FLAME one needs RCPP properly running
# install.packages("RFunctions/flm_0.1.tar.gz",rep=NULL)
# This is based on simulated data which was constructed to mimic
# the structure of the original data

# Load libraries and Data
rm(list = ls())
library(flm)
library(fda)
library(fdapace)
load("flame_CS_sim_data.RData")

# Extract needed objects
T<-flame_CS_Data$T # observation times
Y<-flame_CS_Data$Y # growth index values
N<-dim(Y)[1]; M<-dim(Y)[2] # sample size and obs per curve
X<-flame_CS_Data$X # Scaled abundance values

# First PACE to get data filled in on grid
Lt<-list();Ly<-list()
for(i in 1:N){Lt[[i]] = T[i,];Ly[[i]]=Y[i,]}
FPCA_fit<-FPCA(Ly,Lt,optns = list(userBwCov=5,userBwMu=5))
Y_imp<-fitted(FPCA_fit)
Y_grid<-FPCA_fit$workGrid
matplot(Y_grid,t(Y_imp),type="l")


# Nest run FLAME on centered data
Y_imp_c<-scale(Y_imp,scale=FALSE)
Y_list<-list(time_domain=Y_grid,data=Y_imp_c)
set.seed(1234)
FLAME_fit<-FLAME(Y=Y_list,X=X,type_kernel="sobolev", ratio_lambda=0.01)


# Examine results
# True predictors are 12 35 41 54 61
# Should miss 12th, but capture others.
FLAME_fit$predictors
betas<-FLAME_fit$beta$data[FLAME_fit$predictors,]
matplot(Y_grid,t(betas),type="l")




# Now we turn to penalized regression
# Nothing should come out significant as 
# the data was simulated for FLAME
# and is independent of the diet and FB ratio
Xreg<-flame_CS_Data$Xreg
Xreg_spts<-flame_CS_Data$Xreg_spt
Y.f<-Data2fd(Y_grid,t(Y_imp_c))
Y.f<-Y.f[Xreg_spts]
source("RFunctions/FS_penreg_v2.R")

myfit<-suppressWarnings(FS_penreg(Y.f,Xreg,q=6))
myfit$pvals

# However, we can add an artificial signal 
# to see the effect
beta_fb<-matrix(0,nrow=4,ncol=length(Y_grid))
beta_fb[4,]<-rep(0.02,times=length(Y_grid))
#
signal<-Data2fd(Y_grid,t(Xreg%*%beta_fb))
Y_sig<-Y.f+signal

myfit2<-suppressWarnings(FS_penreg(Y_sig,Xreg,q=6))
myfit2$pvals



