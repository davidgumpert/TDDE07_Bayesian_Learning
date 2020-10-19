# Exercise 1 - Normal model, mixture of normal model with semi-conjugate prior.
data <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab3 - MCMC using Gibbs sampling and Metropolis-Hastings/rainfall.txt", header = FALSE)
rainfallData <- data$V1
n <- length(rainfallData)

# (a) - Gibbs sampler from the joint posterior of Normal model
# Setting parameter values for the full conditional posteriors
# for mu and sigma2
dataMean <- mean(rainfallData)
dataSD <- sd(rainfallData)
hist(rainfallData, breaks = 50, freq = FALSE)
dataSamples <- c()

# Setting arbitrary starting values for mu and sigmaSq
mu <- 1
sigmaSq <- 1

# values representing our prior beliefs
v0 <- 10
sigmaSq0 <- 10
tauSq0 <- 1
  
# defining the gibbs sampling loop
for (i in 1:1000) {
  # Calculate tauSqN, muN and vN
  tauSqN <- 1/((n/sigmaSq) + (1/tauSq0))
  w <- (n/sigmaSq)/((n/sigmaSq) + (1/tauSq0))
  muN <- w * dataMean + (1-w) * mu
  vN <- v0 + n
  
  # every other sample is of mu and the other of sigmaSq
  if(i%%2 == 0) {
    mu <- rnorm(1, muN, sqrt(tauSqN))
  } else {
    sumSqDiffMu <- sum((rainfallData - mu)^2)
    variancechiSq <- ((v0*sigmaSq0 + sumSqDiffMu)/(n+v0))
    sigmaSq <- (vN*variancechiSq)/rchisq(1, vN)
  }
  
  sample <- rnorm(1, mu, sqrt(sigmaSq))
  sampleRow <- c(sample, mu, sigmaSq)
  dataSamples <- rbind(dataSamples, sampleRow)
}

dataNames <- c("Sampled Rainfall", "mu", "sigmaSq")
colnames(dataSamples) <- dataNames

# Plotting to see how values converge towards the true posterior
# mean and sigmaSq
par(mfrow = c(1,2))
plot(dataSamples[1:1000,2], dataSamples[1:1000,3], type="o",
     col="darkgreen", main="Mu and sigma without taking burn-in values into 
     account", xlab="mu", ylab="sigmaSq")
plot(dataSamples[50:1000,2], dataSamples[50:1000,3], type="o",
     col="darkgreen", main="Mu and sigma taking burn-in values into 
     account", xlab="mu", ylab="sigmaSq")
plot(dataSamples[1:1000,1], type="l", col="darkblue", main="Sampled Rainfall")


# (b) - Gibbs sampler from the joint posterior of Mixture of Normal models
# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com

##########    BEGIN USER INPUT #################
# Data options
rawData <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab3 - MCMC using Gibbs sampling and Metropolis-Hastings/rainfall.txt", header = FALSE)
rawData$Rainfall <- rawData$V1
rawData$V1 = c()

x <- as.matrix(rawData['Rainfall'])

# Model options
# Choosing 2 components according to assignement
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
# Choose 10 as our prior mean
muPrior <- rep(10,nComp) # Prior mean of mu
# Choosing 20 as our prior std
tau2Prior <- rep(20,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 500 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.001 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
#ylim <- c(0,2*max(hist(x)$density))

piVector <- c()
sigma2Vector <- c()
muVector <- c()
mixedDensVector <- c()
for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    # hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      # lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    # lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    # legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
    #       col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  sigma2Vector <- rbind(sigma2Vector, sigma2)
  muVector <- rbind(muVector, mu)
  piVector <- rbind(piVector, pi[1,])
  mixedDensVector <- rbind(mixedDensVector, mixDensMean)
}

par(mfrow = c(1,1))
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

#########################    Helper functions    ##############################################

par(mfrow = c(1,2))
plot(muVector[,1],sqrt(sigma2Vector[,1]), type="o", col="darkblue", main="First component with burnin",
     xlab="Mu", ylab="sigma squared")
plot(muVector[100:500,1],sqrt(sigma2Vector[100:500,1]), type="o", col="darkblue", 
     main="First component without burnin", xlab="Mu", ylab="Sigma squared")

plot(muVector[,2],sqrt(sigma2Vector[,2]), type="o", col="darkred", main="Second component with burnin",
     xlab="Mu", ylab="Sigma squared")
plot(muVector[100:500,2],sqrt(sigma2Vector[100:500,2]), type="o", col="darkred", 
     main="Second component without burnin", xlab="Mu", ylab="Sigma squared")

########################### part c ############################################################
# Calculating posterior mean and std for normal model in assignement a
muhata = mean(dataSamples[,2])
sigma2hata = mean(dataSamples[,3])

par(mfrow = c(1,1))
hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = muhata, sd = sqrt(sigma2hata)) ,type = "l", lwd = 2, col = "blue")
#legend("topright", c("Original Data", "Normal model", "Mixed Normal Model"), col= c("black", "darkred", "darkblue"), lty=1, lwd=2)



# EXERCISE 2 - Metropolis Random Walk for Poisson regression
rm(list = ls())
data2 <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab3 - MCMC using Gibbs sampling and Metropolis-Hastings/eBayNumberOfBidderData.txt", header = TRUE)
data2 <- data2[,-2]

# (a)
glmFit <- glm(nBids ~ ., data = data2, family = "poisson")
summary(glmFit)

X = as.matrix(data2[,2:9])
y = data2$nBids


# (b)
library("mvtnorm")
#Defining the posterior function for the poisson model
postPoisson <- function(betaVect, y, X, mu, Sigma) {
  nPara <- length(betaVect)
  linPred <- X%*%betaVect
  
  # Calculating log likelihood for
  logLik <- sum(y*linPred - exp(linPred))
  
  # Likelihood is not finite, stear the optimizer away from here!
  if (abs(logLik) == Inf) logLik = -20000
  
  # Evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), sigma = Sigma, log = TRUE)
  
  return(logLik + logPrior)
} 

data <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab3 - MCMC using Gibbs sampling and Metropolis-Hastings/eBayNumberOfBidderData.txt", header = TRUE)

# Defining values as input for optim
chooseCov <- c(2:10)

X <- as.matrix(data[,1:10])
nPara <- dim(X)[2];
y <- as.vector(data[,1])
covNames <- names(data)[1:length(names(data))];
X <- X[,chooseCov]; # Here we pick out the chosen covariates.
covNames <- covNames[chooseCov];
mu <- as.vector(rep(0,nPara))
Sigma <- 100*solve(as.matrix(t(X))%*%as.matrix(X))
initVal <- as.vector(rep(0,dim(X)[2])); 

OptimResults<-optim(initVal,postPoisson,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing results
postMode <- OptimResults$par
postCov <- -solve(OptimResults$hessian)
names(postMode) <- covNames

print("posterior mode: ")
print(postMode)

print("Posteror covariance matrix: ")
print(postCov)


# (c)
# Defining the log posterior poisson function to be as input into the MH random
# walk fucntion
logPostPoisson <- function(theta, y, X, mu, Sigma) {
  nPara <- length(theta)
  linPred <- X%*%theta
  
  # Calculating log likelihood for
  logLik <- sum(y*linPred - exp(linPred))
  
  # Likelihood is not finite, stear the optimizer away from here!
  if (abs(logLik) == Inf) logLik = -20000
  
  # Evaluating the prior
  logPrior <- dmvnorm(theta, matrix(0,nPara,1), sigma = Sigma, log = TRUE)
  
  return(logLik + logPrior)
} 

MHRandomWalk <- function(postCov, c, logPostFunct, ...) {
  theta <- as.vector(rep(0, dim(X)[2]))
  thetaVector <- c()
  for(i in 1:5000){
    newTheta <- as.vector(rmvnorm(1, theta, c*postCov))
    logPost <- logPostFunct(theta, ...)
    logPostNew <- logPostFunct(newTheta, ...)
    
    # Calculating acceptance level
    alpha <- min(1, exp(logPost - logPostNew))
    alphaCompare <- runif(1,0,1)
    if(alphaCompare > alpha) {
      theta <- newTheta
    }
    thetaVector <- rbind(thetaVector, theta)
  }
  return(thetaVector)
}

c <- 1
res <- MHRandomWalk(postCov, c, logPostPoisson, y, X, mu, Sigma)

# Plot how all coefficients converge (or walk)
par(mfrow = c(1,3))
plot(res[,1], type="o", col="blue", cex=0.2, main="How Constant converge", 
     xlab="iteration", ylab="constant value")    
plot(res[,2], type="o", col="blue", cex=0.2, main="How PowerSeller converge", 
     xlab="iteration", ylab="PowerSeller value")  
plot(res[,3], type="o", col="blue", cex=0.2, main="How VerifyID converge", 
     xlab="iteration", ylab="VerifyID value") 
plot(res[,4], type="o", col="blue", cex=0.2, main="How Sealed converge", 
     xlab="iteration", ylab="Sealed value") 
plot(res[,5], type="o", col="blue", cex=0.2, main="How Minblem converge", 
     xlab="iteration", ylab="Minblem value") 
plot(res[,6], type="o", col="blue", cex=0.2, main="How MajBlem converge", 
     xlab="iteration", ylab="MajBlem value") 
plot(res[,7], type="o", col="blue", cex=0.2, main="How LargNeg converge", 
     xlab="iteration", ylab="LargeNeg value") 
plot(res[,8], type="o", col="blue", cex=0.2, main="How LogBook converge", 
     xlab="iteration", ylab="LogBook value") 
plot(res[,9], type="o", col="blue", cex=0.2, main="How MinBidShare converge", 
     xlab="iteration", ylab="MinBidShare value") 


# (d)
x <- c(1,1,1,1,0,0,0,1,0.5)
resBurn <- res[4001:5000,]
meanBetas = c()

for(i in 1:ncol(resBurn)) {
  meanBetas <- append(meanBetas, mean(resBurn[,i]))
}

poissonRegression <- exp(meanBetas%*%x)[1]
simulation <- rpois(50000, poissonRegression)
par(mfrow=c(1,1))
hist(simulation)

pNoBidders <- sum(simulation==0)/length(simulation)
pNoBidders
