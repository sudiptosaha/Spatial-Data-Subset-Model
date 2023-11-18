#clear the environment
rm(list=ls())

#size of N (20M observed locations and 1000 missing locations)
N <- 20001000

#simulate grid of locations
xnum <- 6667
ynum <- 3000
locations = matrix(0,N,2)
xcoords <- seq(0,1,by=1/(xnum-1))
ycoords <- seq(0,1,by=1/(ynum-1))
for(i in 1:length(xcoords))
  locations[(((i-1)*n)+1):(i*n),1] <- xcoords[i]
locations[,2] <- ycoords

#spectral simulation to generate nu
K = 1000
phi = 3
omega2 = matrix(rcauchy(2*K,0,1/phi), 2, K)
kappa = runif(K,-pi,pi)
nutilde.K=matrix(0,N,1)
for (i in 1:N)
{
  nutilde.K[i] = sqrt(2/K)*sum(cos(locations[i,]%*%omega2 + kappa))
  if(i%%100000==0)
    print(i)
}

#simulate predictor X=(x1,x2)
X <- matrix(rnorm(2*N),N,2)

#true beta
beta_true <- as.matrix(c(2,3))

#simulate latent y
y_latent <- as.matrix(X%*%beta_true + nutilde.K)

#simulate observed y
sigma2_true <- 0.1
y <- rnorm(1,y_latent,sigma2_true)

#set index of missing locations
A <- sample(c(1:N),1000) #number of missing locations = 1000

#final data frame
df <- as.data.frame(cbind(xLoc=locations[,1],yLoc=locations[,2],yLatent=y_latent[,1],yObs=y[,1]))
save.image("data_20M.RData")
