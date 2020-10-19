# Libraries and settings
library(mvtnorm)
library(matlib)
set.seed(12345)

# --------------------------------------------------------------------- #
# (a)
# Import the data
data <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab2 - Polynomial regression and classification with logistic regression/temp.txt", header = TRUE)
head(data)

# Function to generate scaledInvChiSq draws
scaledInvChiSqDraw <-function(NDraws, sigmaSquare, DoF){
  postDraws <- matrix(0,NDraws,1)
  postDraws[,1]<-((DoF)*sigmaSquare)/rchisq(NDraws, DoF)
  return(postDraws)
}

# Function to calculate the temperature as a polynomial regression
calcTemp = function(beta, sigmaSquare){
  temp <- c()
  for (t in data$time) {
    # Draw an error from the normal distribution with mean = 0
    error <- rnorm(1 , 0, sqrt(sigmaSquare))
    # Predict new temperature with quadratic regression
    newtemp <-beta[1]+beta[2]*t+beta[3]*t^2 + error
    # Add prediction of temp in time t to vector of all temp predictions
    temp <- c(temp,newtemp)
  }
  return(temp)
}

# Initial parameters
mu0 <- c(-10, 100, -100)
omega0 <- 0.01*diag(3)
v0 <- 4
sigmaSquared0 <- 1
NDraws <- 10

# Simulate 10 sigma draws over a sequens of 366 days wit hthe initial priors
NDraws <- 10
par(mfrow <- c(1,2))
tempVector <- c()
plot(tempVector, ylim = c(-20,40), xlim = c(0,370),
     main = "Prior from original parameter values", xlab = "Day", ylab = "Temperature")

for (i in 1:NDraws) {
  sigmaDraw <- scaledInvChiSqDraw(1, sigmaSquared0, v0 )
  betadraw <- rmvnorm(1, mean = mu0, sigma = sigmaDraw[1]*solve(omega0))
  temps <- calcTemp(betadraw[1,], sigmaSquared0)
  tempVector <- rbind(tempVector,temps)
  lines(temps)
}

# Set new priors to get better prediction curves
# From the initial priors we see that the variance is quite large.
# To lower this we increase omega0. ALso we make the curvature more steep by
# changing  mu0.
omega0 <- 1*diag(3)
mu0 <- matrix(c(-10,130,-130))
tempVector <- c()
tempinfoVector <- c()
plot(tempVector, ylim = c(-20,40), xlim = c(0,370), main = "Prior from modified parameter values", xlab="Day", ylab="Temperature")

for (i in 1:NDraws) {
  sigmaDraw <- scaledInvChiSqDraw(1, sigmaSquared0, v0 )
  betadraw <- rmvnorm(1, mean = mu0, sigma = sigmaDraw[1]*solve(omega0))
  temps <- calcTemp(betadraw[1,], sigmaSquared0)
  tempVector <- rbind(tempVector,temps)
  lines(temps)
}

# --------------------------------------------------------------------- #

# (b)
# Sampling now from the posterior distribution
# Creating the feature matrix
X <- c()
for(time in data$time) {
  xNew <- c(1,time,time^2)
  X <- rbind(X,xNew)
}

# Defining the new parameters
n <- length(data$temp)
XX <- t(X) %*% X
betaHat <- solve(XX) %*% t(X) %*% data$temp
mun <- solve(XX + omega0) %*% (XX %*% betaHat + omega0 %*% mu0)
omegan <- XX + omega0
vn <- v0 + n
sigmaSquaren <- (v0*sigmaSquared0 + (t(data$temp) %*% data$temp + t(mu0) %*% omega0 %*% mu0 - t(mun) %*% omegan %*% mun))/vn

# Draw beta and sigmasquare from joint posterior distribution and plot
tempVector <- c()
betaDraws <- c()
sigmaDraws <- c()
par(mfrow = c(1,1))
plot(c(), ylim=c(-20,40), col = "blue", xlim = c(0,370), main="Draws from the posterior", xlab="Day", ylab="Temperature")
legend(x = 1, y = 35, c("Temperature data", "Curve over the posterior median"), col = c("blue", "black"), pch = 16)
points(data$temp, col ="blue")

NDraws <- 100
for (i in 1:NDraws) {
  sigmaDraw <- scaledInvChiSqDraw(1, sigmaSquaren, vn)
  sigmaDraws = c(sigmaDraws, sigmaDraw)
  
  betadraw <- rmvnorm(1, mean = mun, sigma = sigmaDraw[1]*solve(omegan))
  betaDraws = rbind(betaDraws, betadraw)
  
  temps <- calcTemp(betadraw[1,], sigmaDraw)
  tempVector <- rbind(tempVector,temps)
}

# Calculating the median betas to be able to calculate the median posterior
medianBeta0 <- median(betaDraws[,1])
medianBeta1 <- median(betaDraws[,2])
medianBeta2 <- median(betaDraws[,3])

# Computing the median posterior of the regression function
medianPosteriorTemp <- medianBeta0 + medianBeta1*data$time + medianBeta2*data$time^2
lines(medianPosteriorTemp)

# Computing the posterior credible interval
lowerLimitVector <- c()
upperLimitVector <- c()
for (i in 1:ncol(tempVector)) {
  lowerLimit <- quantile(tempVector[,i], probs = c(0.025))
  lowerLimitVector <- c(lowerLimitVector, lowerLimit[[1]])
  
  upperLimit <- quantile(tempVector[,i], probs = c(0.975))
  upperLimitVector <- c(upperLimitVector, upperLimit[[1]])
}

lines(lowerLimitVector, col = "red")
lines(upperLimitVector, col = "red")


# Plot histogram of sigma and beta
hist(sigmaDraws, main = "Marginal posterior for sigma", breaks = 20)
hist(betaDraws[,1], main = "Marginal posterior for Beta 0", breaks = 20)
hist(betaDraws[,2], main = "Marginal posterior for Beta 1", breaks = 20)
hist(betaDraws[,3], main = "Marginal posterior for Beta 2", breaks = 10)

# --------------------------------------------------------------------- #

# (c)
# posterior distribution of xTilde
# find this by taking the derivative and setting = 0
xTildeVector <- c()
for (row in 1:nrow(betaDraws)) {
  betaRow <- betaDraws[row,]
  xTilde <- (-betaRow[2]/(2*betaRow[3])*366)
  xTildeVector <- c(xTildeVector,xTilde)
}
hist(xTildeVector, breaks = 10)


# ----------------------------- EXERCISE 2 ---------------------------- #
# (a)
# Importing and prepping data
data <- read.table("/Users/davidgumpert/Desktop/TDDE07/Labs/Lab2 - Polynomial regression and classification with logistic regression/WomenWork.txt", header = TRUE)
head(data)
yLabels <- as.vector(data[,1])
xMatrix <- as.matrix(data[,2:ncol(data)]) #including the bias constant vector
covNames <- names(data)[2:ncol(data)]

# Function to optimize on with optim, i.e. the regression coefficients
LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  #Computes the number of parameters
  nPara <- length(betaVect);
  
  # Makes a linear prediction on the y-value by multiplying the feature matrix
  # with the Beta vector
  linPred <- X%*%betaVect;
  
  # Calculating the log-likelihood                                    
  logLik <- sum( linPred*y - log(1 + exp(linPred)));
  
  # Likelihood is not finite, stear the optimizer away from here!
  if (abs(logLik) == Inf) logLik = -20000; 
  
  # Draw the prior Beta
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # Add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
}

# Variable definition
nParam <- ncol(data) - 1
tauSq <- 10^2
tauSqI <- tauSq * diag(nParam)
intitialBetaVal <-  as.vector(rnorm(dim(xMatrix)[2]))
mu <- as.vector(rep(0,nParam))

# Calling the optim function to optimize my Beta values given the initial 
# Beta vector (can be initialized to anything), the function to optimize on,
# the labels of the data, the training data, mu  and sigma
OptimResults<-optim(intitialBetaVal,
                    LogPostLogistic,
                    gr=NULL,
                    yLabels,
                    xMatrix,
                    mu,
                    tauSqI,
                    method=c("BFGS"),
                    control=list(fnscale=-1),
                    hessian=TRUE)

postMode <- OptimResults$par
postCovariance <- -solve(OptimResults$hessian)
names(postMode) <- covNames
approxPostStd <- sqrt(diag(postCovariance))
names(approxPostStd) <- covNames
multiVarNormDraws <- rmvnorm(1000, postMode, postCovariance)


# investigating the variable NSmallChild
nDraws <- 1000
NSmallChildDraws <- rnorm(nDraws, postMode["NSmallChild"],
                          approxPostStd["NSmallChild"])

# Compute the approximate credible interval of NSmallChild
lowerLimit <- quantile(NSmallChildDraws, 0.025)
upperLimit <- quantile(NSmallChildDraws, 0.975)

# Plotting the distribution
hist(NSmallChildDraws, probability = TRUE)
lines(density(NSmallChildDraws))
abline(v = c(lowerLimit, upperLimit), col = "red")
# The credible interval does not include 0 and has a mean of -1.4 it is believed
# to be an important determinant whether a woman works or not

# --------------------------------------------------------------------- #
# (b)
library(Rlab)
# Function that simulates from the predictive distribution of the response 
# variable in a logistic regression.
predDist <- function(featureVect, postMode, covariance, nDraws) {
  simulatedBetas <- rmvnorm(nDraws, postMode, covariance)
  # calculate approcimate probability
  expValue <- exp(featureVect %*% t(simulatedBetas))
  prob <- expValue/(1 + expValue)
  return(prob)
}

featureVect <- c(1, 10, 8, 10, 1, 40, 1, 1)
predDistVect <- predDist(featureVect, postMode, postCovariance, 1000)
hist(predDistVect, xlab = "Probability to work", main = "Histogram of the probability distribution")
# --------------------------------------------------------------------- #
# (c)
binomial <- 0
for(i in 1:10){
  pred <- ifelse(predDist(featureVect, postMode, postCovariance, 1000)>0.5, 1, 0)
  binomial <- binomial + pred
}
barplot(table(binomial), main="Histogram for the number out of 10 women working")

## CORRECT WAY
pred <- mean(predDistVect)
bin <- rbinom(1000, 10, pred)
hist(bin, xlim=c(0,10), main="Correct way") 

