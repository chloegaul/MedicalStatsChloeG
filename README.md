# MedicalStatsChloeG

#We first want to define a function in R that will give us the minimum sample size n that will give us the desired type I and type II rates. This will give us an idea of what values to choose for the number of patients we recruit.  

#We adapt the type II error formula (for two sample, continuous variable trial) to apply to our binary, one sample case, and rearrange for n to find a formula that gives us a minimum sample size n. 

#For alpha = type 1 error, beta = type 2 error, theta0 = theta under null hypothesis, theta1 = theta under alternate hypothesis, we define a function:

samplesize <- function(alpha, beta, theta0, theta1){
  
  nvarestimate <- (theta0*(1 - theta0))  #n x (variance of estimate), to get numerator of variance.
  
  deltasquared <- (theta1 - theta0)^2    #delta = theta1 -theta0
  
  numerator <- nvarestimate*(qnorm(1-alpha) - qnorm(beta))^2 #numerator of expression for n 
  
  n <- numerator/deltasquared  #Calculate n 
  
  return(n)
  
}

#Now we input the values we have been given, to find the minimum sample size we need.

samplesize(0.05, 0.2, 0.5, 0.7)

#This tells us we need a sample size of at least 38.64098, which we will round up to 39 participants.

#We want to find the type I and type II error for various values of lambda and gamma in our threshold function. This will guide us on the decision of which lambda and gamma to use.

#We can find the type I and type II error through monte carlo methods: Finding the proportion of a set of simulated trials where the null was rejected,

#To find the type I error, we will find the proportion of trials where the null is rejected when theta = theta0. 

#To find the type II error, we find the proportion of trials where the null was rejected, when theta = theta1, and minus this from 1. 


lambda <- seq(0.1, 1, 0.1) #Lambda is between 0 and 1. 
gamma<- seq(0.1 , 1, 0.1) #Gamma needs to be larger than 0

LambGamCombo <- expand.grid(Lambda = lambda, Gamma = gamma) #Find all possible combinations of those values. 

#In our function we want to use a combination of lambda and gamma expressed as x = (lambda, gamma). This needs to be an argument of the function.

ProbofRejectedNull <- function(x, n1, n2, theta){
  
  lambda <- x[1]
  gamma <- x[2]
  
  M_Stage1 <- 10^5 #We want 10^5 simulated trials, in order to have a large sample for Monte Carlo.
  
  y1 <- rbinom(M_Stage1, n1, theta)  #Generate a 'number of responses from first stage of trial' variable for each M in M_Stage1
  
  a1 <- 0.5 + y1 #Work out parameters for posterior, 0.5 has been used for a0
  b1 <- 0.5 + n1 - y1 #Work out parameters for posterior, 0.5 has been used for b0
  
  probfut1 <- pbeta(0.5, a1, b1) #Find probability of futility using the posterior. 
  
  threshold1 <- 1 - lambda*(n1/n2)^gamma
  
  #The trials that pass this stage will then go to the second stage 
  
  M_Stage2 <- sum(probfut1 <= threshold1)
  
  #We now have M_Stage2 trials that are being assessed again. 
  
  y2 <- rbinom(M_Stage2, n2-n1, theta) #Generate a 'number of responses from second stage of trial' variable for each M in M_Stage2
  
  a2 <- 0.5 + y1 + y2 #Updating the posterior 
  b2 <- 0.5 + n2 - y2 - y1 #Updating the posterior  
  
  probfut2 <- pbeta(0.5, a2, b2) #New probability of futility based on new posterior. 
  
  threshold2 <- 1 - lambda*(n2/n2)^gamma #New threshold to compare probability of futility for second stage of trial to.
  
  rejectednulls <- sum(probfut1 <= threshold1 & probfut2 <= threshold2) #The Null is rejected if both thresholds are not crossed.
   
  return(rejectednulls/M_Stage1) #Finding proportion of trials that rejected the null
}

#So we now have a function that will give us the propability of rejecting the null. 

LambGamCombo

TypeII <- apply(LambGamCombo, 1, ProbofRejectedNull, n1 = 30, n2 = 60, theta =0.5)
TypeI <- 1 -apply(LambGamCombo, 1, ProbofRejectedNull, n1 =30, n2 =60, theta =0.7)

TypeIAndTypeII <- cbind(TypeI, TypeII)
TypeIAndTypeII

ChosenCombos <- which(TypeIAndTypeII[ ,1] <= 0.5 & TypeIAndTypeII[ ,2] <= 0.2)
ChosenCombos

ChosenLambdaGamma <- LambGamCombo[ChosenCombos, ]
ChosenLambdaGamma

#So now we have a list of gamma and lambda that will give us the desired operating characteristics.

#We now calculate the expected sample size with each one of these errors, and find the smallest expected sample size. 

ExpectedSampleSize <- function(lambda, gamma, n1, n2) {
  
  M <- 10^5 #M simulations
  
  N <- rep(NA, M) #Creating a vector in which we can store the expected sample sizes.
  
  for (i in 1:M) {
  
    theta <- rbeta(1, 0.5, 0.5) #Generate a value of theta from its prior, defined as a0 = 0.5 and b0 = 0.7
    y1 <- rbinom(1, n1, theta) #Find the number of responses from the trial, using its distribution.
    
    a1 <- 0.5 + y1 #Calculate posterior parameters
    b1 <- 0.5 + n1 - y1 #Calculate posterior paraeters
    
    probfut1_SampleSize <- pbeta(0.5, a1, b1)
    
    threshold1_SampleSize <- lambda * (n1 / n2)^gamma
    
    if (probfut1_SampleSize > threshold1_SampleSize) {
      N[i] <- n1
    } else {
      N[i] <- n2
    }
  }
  
  return(mean(N))
  }

ChosenLambdaGamma

ExpectedSampleSize(lambda = 0.8, gamma = 0.1, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.1, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.8, gamma = 0.2, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.2, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.8, gamma = 0.3, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.3, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.8, gamma = 0.4, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.4, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.5, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.6, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.7, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.8, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 0.9, n1 = 30, n2 =60)
ExpectedSampleSize(lambda = 0.9, gamma = 1, n1 = 30, n2 =60)

#For Evaluation, we want to look into what happens if we do this with two interim analyses rather than just one. 
#We want to adjust the code to allow for 2 interim analyses, so 3 sets of patients we are testing. 
#We then want to see what happens to elements such as: Possible lambdas and gammas to use and minimum expected sample size. 
#If it requires less participants to achieve the same power, we may conclude that using two interim analysis could be more efficient.












