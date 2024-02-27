#######################################################################################################
#                                                                                                     #
#       This is the code to generate simulated data for Spatial Data Subset Model (SDSM).             #
#       Note that "Nn" stores the total size of observed and missing data added,                      #
#       and "m" stores the size of only prediction/missing data.                                      #
#       To generate data with different size, modify the code accordingly as mentioned below:         #
#       (i)   Line 18: change the total size of the data (observed + missing).                        #
#       (ii)  Line 19: change the size of prediction/missing data.                                    #
#       (iii) Line 22: change the grid size in the direction of x-axis.                               #
#       (iv)  Line 23: change the grid size in the direction of y-axis.                               #
#                                                                                                     #
#######################################################################################################

#clear the environment
rm(list=ls())

#size of the data (20M observed locations and 1000 missing locations)
Nn <- 20001000
m <- 1000

#simulate grid of locations
xnum <- 6667
ynum <- 3000
locations <- matrix(0,Nn,2)
xcoords <- seq(0,1,by=1/(xnum-1))
ycoords <- seq(0,1,by=1/(ynum-1))
for(i in 1:length(xcoords))
  locations[(((i-1)*ynum)+1):(i*ynum),1] <- xcoords[i]
locations[,2] <- ycoords

#spectral simulation setup to generate nu
K <- 1000
phi <- 3
omega2 <- matrix(rcauchy(2*K,0,1/phi), 2, K)
kappa <- runif(K,-pi,pi)

#function to generate nu
func_nu <- function(loc){
  return(sqrt(2/K)*sum(cos(loc%*%omega2 + kappa)))
}

#simulate nu
nutilde.K <- as.matrix(apply(locations,1,func_nu))

#simulate predictor X=(x1,x2)
X <- matrix(rnorm(2*Nn),Nn,2)

#true beta
beta_true <- as.matrix(c(2,3))

#simulate latent y
y_latent <- as.matrix(X%*%beta_true + nutilde.K)

#simulate observed y
sigma2_true <- 0.1
y <- rnorm(Nn,y_latent,sigma2_true)

#set index of missing locations
A <- sample(c(1:Nn),m)

#enforce missing locations into the data
y[A] <- NA

#final data frame
sim.data <- as.data.frame(cbind(Lon=locations[,1],Lat=locations[,2],Temp=y))

#save environment
save.image("SimulatedData.RData")
