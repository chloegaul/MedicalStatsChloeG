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

#We do this through monte carlo methods: Finding the proportion of a set of simulated trials where the null was rejected.

#To find the type I error, we will find the proportion of trials where the null is rejected when theta = theta0. 

#To find the type II error, we find the proportion of trials where the null was rejected, when theta = theta1, and minus this from 1. 




