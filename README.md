# MedicalStatsChloeG

#Define a function in R that outputs the minimum sample size n required for the desired, given type I and type II rates, for the purposes of guiding our decision of n1 and n2.  

#Adapt the type II error formula (for two sample, continuous variable trial) to apply to our binary, one sample case, and rearrange for n.

samplesize <- function(alpha, beta, theta0, theta1){ #Arguments of function: type I error, type II error, null theta and alternate theta
  
  nvarestimate <- (theta0*(1 - theta0))  #n x (variance of estimate), to get numerator of variance.
  
  deltasquared <- (theta1 - theta0)^2    #delta = theta1 -theta0
  
  numerator <- nvarestimate*(qnorm(1-alpha) - qnorm(beta))^2 #numerator of expression for n 
  
  n <- numerator/deltasquared  #Calculate n 
  
  return(n)
  
}

#Input given values to find the minimum n2 (total sample size) we need to use.

samplesize(0.05, 0.2, 0.5, 0.7)

#So this suggests we need at least n = 38.64098, which we will round up to 39 participants. We choose n1 = 25, n1 = 50. 

#Find type I and type II error for various values of lambda and gamma, for the purposes of guiding choices of lambda and gamma. 

#We use monte carlo methods: Finding the proportion of a set of simulated trials where the null was rejected. 

lambda <- seq(0.1, 1, 0.05) #Lambda is between 0 and 1, create set of 10 values for lambda in this range
gamma<- seq(0.1 , 1, 0.05) #Gamma needs to be larger than 0, create set of 10 values for gamma in this range 

LambGamCombo <- expand.grid(Lambda = lambda, Gamma = gamma) #Find all possible combinations of those values. 

nrow(LambGamCombo)

ProbofRejectedNull <- function(x, n1, n2, theta){ #combinations inputted into the function as x = (lambda, gamma)
  
  lambda <- x[1] 
  gamma <- x[2] 
  
  M_Stage1 <- 10^5 #We want 10^5 simulated trials, in order to have a large sample for Monte Carlo.
  
  y1 <- rbinom(M_Stage1, n1, theta)  #Generate a 'number of responses from first stage of trial' variable for each M in M_Stage1
  
  a1 <- 0.5 + y1 #First parameter for posterior theta distribution, 0.5 has been used for a0
  b1 <- 0.5 + n1 - y1 #Second parameter for posterior theta distribution, 0.5 has been used for b0
  
  probfut1 <- pbeta(0.5, a1, b1) #Find probability of futility using the posterior. 
  
  threshold1 <- 1 - lambda*(n1/n2)^gamma #Use lambda and gamma to calculate threshold for probfut1
  
  M_Stage2 <- sum(probfut1 <= threshold1) #The trials that pass this stage will then go to the second stage to be assessed again.
  
  y2 <- rbinom(M_Stage2, n2-n1, theta) #Generate a 'number of responses from second stage of trial' variable for each M in M_Stage2
  
  a2 <- 0.5 + y2 #Updating the posterior 
  b2 <- 0.5 + (n2 - (y1 + y2)) #Updating the posterior  
  
  probfut2 <- pbeta(0.5, a2, b2) #New probability of futility based on new posterior. 
  
  threshold2 <- 1 - lambda*(n2/n2)^gamma #New threshold to compare to probability of futility for second stage of trial.
  
  rejectednulls <- sum(probfut1 <= threshold1 & probfut2 <= threshold2) #The Null is rejected if both thresholds are not crossed.
   
  return(rejectednulls/M_Stage1) #Finding proportion of trials that rejected the null
}


#To find the type I error, we will find the proportion of trials where the null is rejected when theta = theta0. 

#To find the type II error, we find the proportion of trials where the null was rejected, when theta = theta1, and minus this from 1.

TypeI <- apply(LambGamCombo, 1, ProbofRejectedNull, n1 = 25, n2 = 50, theta =0.5) #type I error is proportion of trials where null is rejected when theta = theta0.

RejectedUnderTheta1 <- apply(LambGamCombo, 1, ProbofRejectedNull, n1 =25, n2 =50, theta =0.7)

TypeII <- 1 - RejectedUnderTheta1 #Type II error is proportion of trials where the null was rejected, when theta = theta1, subtracted from 1.
TypeII

TypeIAndTypeII <- cbind(TypeI, TypeII) #Form matrix of type 1 and type 2 errors for combinations of lambda and gamma. 

ChosenCombos <- which(TypeIAndTypeII[ ,1] <= 0.05 & TypeIAndTypeII[ ,2] <= 0.2) #Choose rows that satisfy constraints on type I and II errors

ChosenLambdaGamma <- LambGamCombo[ChosenCombos, ] #Extract combinations that gives us desired error rates.

nrow(ChosenLambdaGamma) #Calculate how many combinations of lambda and gamma give us desirable error rates. 

ExpectedSampleSize <- function(lambda, gamma, n1, n2) { #Define a function that calculates expected sample size.
  
  M <- 10^5 #M simulations
  
  N <- rep(NA, M) #Creating a vector in which we can store the expected sample sizes.
  
  for (i in 1:M) {
  
    theta <- rbeta(1, 0.5, 0.5) #Generate a value of theta from its prior, defined as a0 = 0.5 and b0 = 0.5
    
    y1 <- rbinom(1, n1, theta) #Find the number of responses from the trial, using its distribution.
    
    a1 <- 0.5 + y1 #Calculate posterior parameters
    b1 <- 0.5 + n1 - y1 #Calculate posterior paraeters
    
    probfut1_SampleSize <- pbeta(0.5, a1, b1) #Calculate probability of futility based on posterior
    
    threshold1_SampleSize <- 1 - lambda * (n1 / n2)^gamma #Calculate threshold, based on sample size so far
    
    if (probfut1_SampleSize > threshold1_SampleSize) {
      N[i] <- n1
    } else {
      N[i] <- n2
    } #Assign sample size for that trial to the each row of empty matrix. N takes 2 values: n1 if trial is stopped after interim analysis, and n2 otherwise.
  }
  
  return(mean(N))
  }

ChosenLambdaGamma 

#Apply ExpectedSampleSize function using these combinations of lambda and gamma

ExpectedSampleSize(lambda = 0.3, gamma = 0.1, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.1, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.15, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.2, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.25, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.3, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.35, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.4, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.45, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.5, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.55, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.6, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.65, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.7, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.75, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.8, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.85, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.9, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 0.95, n1 = 25, n2 =50)

ExpectedSampleSize(lambda = 0.35, gamma = 1, n1 = 25, n2 =50)

 

#Adjust the code for type I and type II errors to allow for 2 interim analyses. 

#Instead of n1 and n2, we now have n1, n2 and n3. where n2 = n3-n1 and n1 = n3-n2. n3 is the total number of participants. 

ProbofRejectedNull2Interims <- function(x, n1, n2, n3, theta){ #The trial now has 3 sets of participants therefore 3 arguments for sample size in the function.
  
  lambda <- x[1]
  gamma <- x[2]
  
  M_Stage1 <- 10^5 #We want 10^5 simulated trials, in order to have a large sample for Monte Carlo.
  
  y1 <- rbinom(M_Stage1, n1, theta)  #Generate a 'number of responses from first stage of trial' variable for each M in M_Stage1
  
  a1 <- 0.5 + y1 #Work out parameters for posterior, 0.5 has been used for a0
  b1 <- 0.5 + n1 - y1 #Work out parameters for posterior, 0.5 has been used for b0
  
  probfut1_2Interims <- pbeta(0.5, a1, b1) #Find probability of futility using the posterior. 
  
  threshold1_2Interims <- 1 - lambda*(n1/n3)^gamma
  
  #The trials that pass this stage will then go to the second stage 
  
  M_Stage2 <- sum(probfut1_2Interims <= threshold1_2Interims)
  
  #We now have M_Stage2 trials that are being assessed again. 
  
  y2 <- rbinom(M_Stage2, n2-n1, theta) #Generate a 'number of responses from second stage of trial' variable for each M in M_Stage2
  
  a2 <- 0.5 + y2 #Updating the posterior 
  b2 <- 0.5 + (n2 - (y1 + y2)) #Updating the posterior  
  
  probfut2_2Interims <- pbeta(0.5, a2, b2) #New probability of futility based on new posterior. 
  
  threshold2_2Interims <- 1 - lambda*(n2/n3)^gamma #New threshold to compare probability of futility for second stage of trial to.
  
  M_Stage3 <- sum(probfut2_2Interims <= threshold2_2Interims) #M_Stage3 trials pass the stage and recruit more patients.
  
  y3 <- rbinom(M_Stage3, n3-n2, theta) #Generate a 'number of responses from second stage of trial' variable for each M in M_Stage3
  
  a3 <- 0.5 + y3 #Updating the Posterior
  b3 <- 0.5 + (n3 - (y1 + y2 + y3)) #Updating the Posterior
  
  probfut3_2Interims <- pbeta(0.5, a3, b3) #New probability of futility based on new posterior. 
  
  threshold3_2Interims <- 1 - lambda*(n3/n3)^gamma #New threshold to compare probability of futility for second stage of trial to
  
  rejectednulls2Interims <- sum(probfut1_2Interims <= threshold1_2Interims & probfut2_2Interims <= threshold2_2Interims & probfut3_2Interims <= threshold3_2Interims) #The Null is rejected if all thresholds are not crossed at each step.
  
  return(rejectednulls2Interims/M_Stage1) #Finding proportion of trials that rejected the null
}

#Calculate the type I and type II errors using the trials we have simulated. 

TypeI_2Interims <- apply(LambGamCombo, 1, ProbofRejectedNull2Interims, n1 = 20, n2 = 35, n3 = 50, theta =0.5) #type I error is proportion of trials where null is rejected when theta = theta0.

RejectedUnderTheta1_2Interims <- apply(LambGamCombo, 1, ProbofRejectedNull2Interims, n1 =20, n2 =35, n3 = 50, theta =0.7)
TypeII_2Interims <- 1 - RejectedUnderTheta1 #Type II error is proportion of trials where the null was rejected, when theta = theta1, subtracted from 1.

TypeIAndTypeII_2Interims <- cbind(TypeI_2Interims, TypeII_2Interims)

ChosenCombos2Interims <- which(TypeIAndTypeII_2Interims[ ,1] <= 0.05 & TypeIAndTypeII_2Interims[ ,2] <= 0.2)

ChosenLambdaGamma2Interims <- LambGamCombo[ChosenCombos2Interims, ] #Extract values of lambda and gamma that give us desired error rates

nrow(ChosenLambdaGamma2Interims) #Calculate how many combinations of lambda and gamma will give desired error rates.

#Adapt function for expected sample size for 2 interim analyses. 

  ExpectedSampleSize_2Interims <- function(lambda, gamma, n1, n2, n3) {
  
  M <- 10^5 #M simulations
  
  N <- rep(NA, M) #Creating a vector in which we can store the expected sample sizes.
  
  for (i in 1:M) {
    
    theta <- rbeta(1, 0.5, 0.5) #Generate a value of theta from its prior, defined as a0 = 0.5 and b0 = 0.5
    
    y1 <- rbinom(1, n1, theta) #Generate a number of responses from first stage of trial
    
    a1 <- 0.5 + y1 #Calculate posterior parameters
    b1 <- 0.5 + n1 - y1 #Calculate posterior paraeters
    
    probfut1_2InterimsSS <- pbeta(0.5, a1, b1) #Calculate probability of futility using posterior distribution of theta
    
    threshold1_2InterimsSS <- 1 - lambda * (n1 / n3)^gamma #Calculate threshold to compare probfut1_2InterimsSS to
    
    y2 <- rbinom(1, n2-n1, theta) #Generate number of responses from second stage of trial
    
    a2 <- 0.5 + y2 #Calculate posterior parameters
    b2 <- 0.5 + (n2 - (y1 + y2)) #Calculate posterior paraeters
    
    probfut2_2InterimsSS <- pbeta(0.5, a2, b2) #Calculate probability of futility using updated posterior distribution of theta
    
    threshold2_2InterimsSS <- 1 - lambda * (n2 / n3)^gamma #Calculate threshold to compare probfut2_2InterimsSS to
    
    if (probfut1_2InterimsSS > threshold1_2InterimsSS) { #Stop trial if probability exceeds threshold, only n1 patients given treatment
      N[i] <- n1
    } 
    
    else{
    
    if (probfut2_2InterimsSS > threshold2_2InterimsSS) { #Stop trial if probability exceeds second threshold, n2 patients given treatment
      N[i] <- n2
    }
    else{
    N[i] <- n3 #First two thresholds are passed, so full sample size, n3 patients are given treatment. 
    }
    } 
  
  }
  
  
  return(mean(N)) # Return the estimated expected sample size, using monte Carlo Methods. 
}

ChosenLambdaGamma2Interims

nrow(ChosenLambdaGamma2Interims)


#We test out 3 values to see what the expected sample size seems to be generally 

ExpectedSampleSize_2Interims(lambda = 0.25 , gamma = 0.15, n1 = 20, n2 = 35, n3 = 50)

ExpectedSampleSize_2Interims(lambda = 0.3 , gamma = 1, n1 = 20, n2 = 35, n3 = 50)

ExpectedSampleSize_2Interims(lambda = 0.1 , gamma = 0.8, n1 = 20, n2 = 35, n3 = 50)











