# MedicalStatsChloeG

#We first want to define a function in R that will give us the minimum sample size n that will give us the desired type I and type II rates. This will give us an idea of what values to choose for the number of patients we recruit.  

#We adapt the type II error formula (for two sample, continuous variable trial) to apply to our binary, one sample case, and rearrange for n to find a formula that gives us a minimum sample size n. 

#For alpha = type 1 error, beta = type 2 error, theta0 = theta under null hypothesis, theta1 = theta under alternate hypothesis, we create a function:




