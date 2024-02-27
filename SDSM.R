#######################################################################################################
#                                                                                                     #
#       This is the code for Spatial Data Subset Model (SDSM).                                        #
#       This code can run as it is with the dataset simulated using the DataSimulator.R code.         #
#       This code generates delta from SRS without replacement (refer to line 125).                   #
#       The choices for subsample sizes can be changed in line 46.                                    #
#       The choice for the sampling design can be changed in line 125.                                #
#       To run the code with a different dataset, modify the code accordingly as mentioned below:     #
#       (i)   Line 40: change the dataset to be loaded.                                               #
#                      Make sure that "Nn" stores the total size of the data (observed + missing),    #
#                      "A" stores the set of index of prediction/missing locations/data,              #
#                      "locations" stores all the locations (observed + missing),                     #
#                      "X" stores the predictors,                                                     #
#                      "y_latent" stores the latent process,                                          #
#                      and "y" stores the data.                                                       #
#       (ii)  Line 46: change the choices for subsample sizes (length of mult_n should be >=1).       #
#       (iii) Line 112: change the support of phi (if needed).                                        #
#       (iv)  Line 125: change the sampling design (if needed).                                       #
#                                                                                                     #
#######################################################################################################

#clear the environment
rm(list=ls())

#load libraries
library(MASS)
library(pscl)
library(mvtnorm)
library(coda)
library(Metrics)
library(scoringRules)
library(scoringutils)
library(LaplacesDemon)
library(bayestestR)
library(parallel)
library(foreach)
library(doParallel)

#load data
load("SimulatedData.RData")

#number of missing locations
m <- length(A)

#set multiple subsample sizes (n)
mult_n <- c(10,20,50,100,200,500,1000,2000)

#number of predictors
p <- dim(X)[2]

#size of MCMC
G <- 5000

#burn-in
g0 <- 1001

#index of training set (non-missing loactions)
index_training <- c(1:Nn)[-A]

#data frame to store results
df_result <- data.frame(n=mult_n,mae=rep(0,length(mult_n)),
                        rmse=rep(0,length(mult_n)),crps=rep(0,length(mult_n)),
                        int=rep(0,length(mult_n)),cvg=rep(0,length(mult_n)),
                        model_time=rep(0,length(mult_n)),pred_time=rep(0,length(mult_n)))

#3D matrix to store mean and variance of the posterior samples of latent y
what_mat <- array(0,dim=c(m,2,length(mult_n))) #what_mat[,1,] stores mean and what_mat[,2,] stores variance

#function for updating phi
func_phi <- function(ph){
  Hphi_val <-  sigma2[,i]*cov_func(ph,DistMat[(m+1):(m+n),(m+1):(m+n)],n)
  return(dmvnorm(as.vector(nu_delta), matrix(0,n,1), Hphi_val, log=FALSE))
}

#function to sample phi
sample_phi <- function(prob_vals){
  if(sum(prob_vals)>0)
  {
    if(sum(prob_vals)!=Inf)
    {
      prob_vals <- prob_vals/sum(prob_vals)
      return(sample(phis, 1, prob = prob_vals))
    }
    else
      return(sample(phis[which(prob_vals==Inf)],1))
  }
  else
    return(sample(phis,1))
}

#set number of clusters to perform prediction in parallel
cl <- makeCluster(detectCores()-1)

#loop over multiple n
for(k in 1:length(mult_n))
{
  #start clusters
  registerDoParallel(cl)
  
  #subsample size
  n <- mult_n[k]
  
  #matrix to store predictions
  w_pred <- matrix(0,nrow=m,ncol=G)
  
  #initialization
  tau2 <- matrix(rigamma(1,1,1),nrow=1,ncol=G)
  sigma2 <- matrix(rigamma(1,1,1),nrow=1,ncol=G)
  sigma2_beta <- matrix(rigamma(1,1,1),nrow=1,ncol=G)
  beta <- matrix(rnorm(p,0,sqrt(sigma2_beta[,1])),nrow=p,ncol=G)
  phi <- matrix(1,nrow=1,ncol=G)
  phis <- seq(0.001,1,0.001) #support of discrete uniform for phi
  
  #track CPU time
  time_model <- 0
  time_pred <- 0
  
  #Gibbs Sampler
  for(i in 2:G)
  {
    #start time for model time
    st <- Sys.time()
    
    #sample delta of size n, using SRS without replacement (different sampling strategies can be used here)
    delta_i <- sample(index_training,n)
    
    #Calculating H(phi) which is an (m+n)*(m+n) matrix
    DistMat <- as.matrix(dist(locations[c(A,delta_i),], method="euclidean"))
    Hphi <- exp(-phi[,i-1]*DistMat)
    Hphi_delta <- Hphi[(m+1):(m+n),(m+1):(m+n)]
    Hphi_delta_inv <- solve(Hphi_delta)
    Hphi_m <- Hphi[1:m,(m+1):(m+n)] #H_m(phi) matrix
    
    #update nu
    nuCovar <- solve((1/tau2[,i-1])*diag(n) + (1/sigma2[,i-1])*Hphi_delta_inv)
    nuMean <- (1/tau2[,i-1])*(nuCovar%*%(y[delta_i,]-X[delta_i,]%*%beta[,i-1]))
    nu_delta <- as.matrix(mvrnorm(1,nuMean,nuCovar,tol=1e-3))
    
    #update beta
    betaCovar <- solve((1/tau2[,i-1])*t(X[delta_i,])%*%X[delta_i,] + (1/sigma2_beta[,i-1])*diag(p))
    betaMean <- (1/tau2[,i-1])*(betaCovar%*%t(X[delta_i,])%*%(y[delta_i,]-nu_delta))
    beta[,i] <- as.matrix(mvrnorm(1,betaMean,betaCovar,tol=1e-3))
    
    #update tau2
    a_tau2 <- 1+(n/2)
    b_tau2 <- 1 + 0.5*(t(y[delta_i,]-(X[delta_i,]%*%beta[,i])-nu_delta)%*%(y[delta_i,]-(X[delta_i,]%*%beta[,i])-nu_delta))
    tau2[,i] <- rigamma(1,a_tau2,b_tau2)
    
    #update sigma2
    a_sigma2 <- 1+(n/2)
    b_sigma2 <- 1 + 0.5*(t(nu_delta)%*%Hphi_delta_inv%*%nu_delta)
    sigma2[,i] <- rigamma(1,a_sigma2,b_sigma2)
    
    #update sigma2_beta
    a_sigma2_beta <- 1+(p/2)
    b_sigma2_beta <- 1 + 0.5*(t(beta[,i]%*%beta[,i]))
    sigma2_beta[,i] <- rigamma(1,a_sigma2_beta,b_sigma2_beta)
    
    #update phi
    prob_vals <- sapply(phis,func_phi)
    phi[,i] <- sample_phi(prob_vals)
    
    #track model time
    time_model <- time_model + as.numeric(Sys.time()-st,unit="secs")
    
	  #prediction step
    if(i>=g0)
    {
      #start time for prediction time
      st1 <- Sys.time()
      
      #exporting variables to cluster for parallel operation
      clusterExport(cl,varlist = c("m","n","i","Hphi_m","Hphi_delta_inv","nu_delta","sigma2"))
      
  	  #updating nu_A
      nu_A <- parSapply(cl,1:m,func <- function(j){
      
        #calculate z=h_{mi.}'*H_{delta}(phi)
  		  z <- matrix(Hphi_m[j,],1,n)%*%Hphi_delta_inv
  		
  		  #mean=z*nu_delta, variance is sigma2*(h_{Aii}-z*h_{mi.}), where h_{Aii}=1
        return(rnorm(1,as.numeric(z%*%nu_delta),sqrt(sigma2[,i-1]*(1-as.numeric(z%*%as.matrix(Hphi_m[j,]))))))
      })
	  
      #calculate predicted y
      w_pred[,i] <- X[A,]%*%beta[,i] + as.matrix(nu_A)
      
      #track prediction time
      time_pred <- time_pred + as.numeric(Sys.time()-st1,unit="secs")
    }
    
    #print iteration number to keep track
    print(paste("Iteration:",i))
  }
  
  #mean of posterior samples
  what_mat[,1,k] <- apply(w_pred[,g0:G],1,mean)
  
  #variance of posterior samples
  what_mat[,2,k] <- apply(w_pred[,g0:G],1,var)
  
  #calculate MAE
  maeval <- mae(y_latent[A,],what_mat[,1,k])
  
  #calculate RMSE
  rmseval <- rmse(y_latent[A,],what_mat[,1,k])
  
  #calculate CRPS
  sd_pred <- apply(w_pred[,g0:G],1,sd) #sd of predictions
  crpsval <- mean(crps.numeric(y_latent[A,],family="normal",mean=what_mat[,1,k],sd=sd_pred))
  
  #calculate INT
  lval <- qnorm(0.05/2,mean=what_mat[,1,k],sd=sd_pred) #lower limit of 95% credible interval
  uval <- qnorm(1-0.05/2,mean=what_mat[,1,k],sd=sd_pred) #upper limit of 95% credible interval
  intval <- mean(interval_score(y_latent[A,],lower=lval,upper=uval,interval_range = rep(95,length(A)),weigh=F))
  
  #calculate CVG
  cvgval <- mean((lval <= y_latent[A,]) & (y_latent[A,] <= uval))
  
  #save result
  df_result[k,] <- c(n,maeval,rmseval,crpsval,intval,cvgval,time_model,time_pred)
  
  #print n to keep track
  print(paste("Subsample Size:",n))
}

#stop cluster
stopCluster(cl)

#save environment
save.image("SDSM_output.RData")
