# -------------------------------------------------------------------------
#  1- Run Bayesian Approach with Weibull distribution.
# -------------------------------------------------------------------------  



## Scenario: "t ~ Weibull(2,4), c ~ Exp(0.02), p=0.20, n=200"

### 1-1 R codes:

{
  rm(list = ls())
  # Install packages:survival & R2openBUGS.
  library(survival)
  library(R2OpenBUGS)
  # Set working directory and modelfile.
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Generate Data
  set.seed(12345)
  n = 200  # n=100; 200; 300
  x = rep(0:1, c(0.50*n, 0.50*n))  # Weibull scale parameter related to x.
  table(x)
  shape = 2  # Shape: 0.5; 1; 2
  b = c(-3.2, 1) # set b1 and b2 with table 2 in paper.
  lambda = exp(b[1] + b[2]*x)   # Link the parameter to covariate x .
  summary(lambda)
  scale = lambda^(-1/shape)   # Since weibull formulla in winbugs is different to R, we need to convert
  # formula to get similar results.
  summary(scale)   # Mean scale parameter is near to 4.
  #Generate Observed time
  y = rweibull(n,shape, scale )    
  summary(y)
  range(y)
  # Generate censored time
  delta1 = rep(1,n)   # to make censored data
  cen = rexp(n,0.06)             # Censored time
  delta = as.numeric(y < cen)
  cenper = 1 - mean(delta); cenper   # Get percent of censoring
  # Merge observed and censored time.
  z = pmin(y,cen)  # to select observed time or censored time. Every one that is lesser than other.
  # make variable "t" as observed time and variable "c" as censored time to use in BUGS.
  t <- ifelse(delta == 1, z, NA)
  c <- ifelse(delta == 1, 0, z)
  # Run model in BUGS.
  modeltext = "model {
      for(i in 1:n){
        t[i] ~ dweib(shape,lambda[i])C(c[i], )
      	log(lambda[i]) <- b[1]+b[2]*x[i]
        cim[i] <- step(c[i]-1.0E-5)*pow(log(2)/lambda[i]+pow(c[i],shape), 1/shape)
        }
      	# priors
      	shape ~ dgamma(0.01,0.01)  # Non-informative prior
      	for(j in 1:2) {b[j]~dnorm(0,0.01)}		
      }
      "
  # write BUGS output into file.
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list(n = n, x = x, t = t, c = c)
  modelinit = list(list(b = rep(0,length(b)), shape = shape))
  param = c("shape","b","cim")
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 10
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  # output ----------------------------------------
  bugsOut$DIC
  # Which records is censored:
  ic = which(delta==0); ic; length(ic)
  # Dimension of output:
  dim(bugsOut$sims.array)
  # Describe censored simulations.
  bugsOut$summary[c(1:3,3+ic),c(1,2)] 
  # Describe parameter simulations:
  parsim1 = bugsOut$sims.array[,1,1:3]   #parameter simulation
  parsim1[1:5,]  # Only five rows of 10.000 simulation for parameters.
  # print median of simulations for every censor that replaced.
  bugsOut$median$cim[ic]  
  impsim = bugsOut$sims.array[,1,3+ic]  # 10000 imputations for every censor.
  impsim[1:3,]
  impsim[ic] 
  timp=t
  #
}


#-------------------------------------------------------------------------------
###1-2     Figure 3-21(H) in the thesis.

{
  # Kaplan-Meier Curve:
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("H                  t~Weibull(2,4), c~Exp(0.06), p=0.20, n=200") ) 
  
  
  # Curve with Median of Simulated Times
  # output ----------------------------------------
  # imputation      h=hat
  bh = bugsOut$mean$b; bh
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = exp(bh[1] + bh[2]*x); lambdah #every person has specific lambda because it has specific X.
  scaleh = lambdah^(-1/shapeh); scaleh
  # Compute median of Simulations.
  zmed = qweibull(.5*pweibull(cen,shapeh,scaleh, lower.tail = FALSE),shapeh, scaleh, lower.tail = FALSE)
  
  zimp = rep(NA,n)
  zimp[ic] = zmed[ic]
  zimp[-ic] = z[-ic]  # zimp = failure times+imputed censored times
  
  curve2 = survfit(Surv(zimp,delta1) ~ x); curve2     # Bayesian Imputation
  lines(curve2, mark.time = TRUE, col = "Blue", lty = 1)
  
  #Curve without Censored Times 
  tOC = z[delta==1]   #time omitting censored
  deltaOC = rep(1, length(tOC))
  curve3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); curve3       # Omitting_Censored
  lines(curve3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
}


# -------------------------------------------------------------------------


###1-3 Figure 3-23(B) in the thesis.


{  
  km1 = survfit(Surv(z,delta) ~ x); km1
  plot(km1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
       main = paste("B                  t~Weibull(2,4), c~Exp(0.06), p=0.20, n=200"))  #KM_Estimation
  
  # simulation 
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    #Sys.sleep(.5)
  }
  timp[ic] <- colMeans(impsim)
  kmmean = survfit(Surv(timp,delta1) ~ x)
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  lines(km1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
        main = paste("B                  t~Weibull(2,4), c~Exp(0.06), p=0.20, n=200"))  #KM_Estimation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
}  




# -------------------------------------------------------------------------
#     2- Run Bayesian Approach with Birnbaum-Saunders distribution.
# -------------------------------------------------------------------------  

##  Scenario: "t ~ BS(2,4), c ~ Exp(0.02), p=0.20, n=200"

### 2-1 R codes:


{
  rm(list = ls())
  library(survival)
  library(R2OpenBUGS)
  # Set working directory and modelfile.
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  
  # generate Data
  set.seed(12345)
  n = 200 # n=100; 200; 300
  x = rep(0:1, c(0.50*n, 0.50*n))  # BS scale parameter related to x.
  table(x)
  shape = 2    # Shape: 0.5; 1; 2                   
  b = c(1, 0.65)  # set b1 and b2 with table 4 in paper.
  lambda = exp(b[1] + b[2]*x)
  summary(lambda)
  scale = lambda             
  # Define rbn to generate numbers from BS distribution.
  rbn <- function(n, shape, scale){    # shape = a, scale = b
    x <- rnorm(n, 0, shape/2)
    t <- scale * (1 + 2 * x^2 + 2 * x * sqrt(1 + x^2))
    return(t)
  }
  #Generate Observed time
  y <- rbn(n, shape, lambda)
  # Generate censored time
  delta1 = rep(1,n)
  cen = rexp(n,0.021)   # repeat to get 20% censoring.
  delta = as.numeric(y < cen)
  cenper = 1 - mean(delta); cenper   # % censoring
  # Merge observed and censored time.
  
  z = pmin(y,cen)
  # make variable "t" as observed time and variable "c" as censored time to use in BUGS.
  t <- ifelse(delta == 1, z, NA)
  c <- ifelse(delta == 1, 0, z)
  # Run model in BUGS.
  # model---------------------------------------- this part is write as BUGs code
  modeltext = "model {
  for(i in 1:n){
	  t[i] ~ dbs(shape, lambda[i])C(c[i], )
	  log(lambda[i]) <- b[1]+b[2]*x[i]
  
  
  cim[i] <- step(c[i]-1.0E-5) * lambda[i] * pow( (shape/2) * 0.5875441*logit(0.5 + 0.5 * phi((1/shape) * 
            (sqrt(c[i]/lambda[i]) - sqrt(lambda[i]/c[i])))) + sqrt( 1 + pow( (shape/2) * 
            0.5875441 * logit(0.5 + 0.5 * phi((1/shape) * (sqrt(c[i]/lambda[i]) -
            sqrt(lambda[i]/c[i])))),2)) ,2)
  
  }
	# priors
	shape ~ dgamma(0.01,0.01)
	for(j in 1:2) {b[j]~dnorm(0,0.01)}		
  }
  "
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list(n = n, x = x, t = t, c = c)
  modelinit = list(list(b = rep(0,length(b)), shape = shape))
  
  param = c("shape","b","cim")
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 1
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  # output ----------------------------------------
  bugsOut$DIC
  ic = which(delta==0); ic
  dim(bugsOut$sims.array)
  bugsOut$summary[c(1:3,3+ic),c(1,2)]  #3+ic means that after third row show rows given ic(censored data).
  parsim = bugsOut$sims.array[,1,1:3]   #parameter simulation
  impsim = bugsOut$sims.array[,1,3+ic]  # imputation simulation: cim=exp(b1+b2*x), x=0,1
  mean(impsim[3,])
  
  small = as.numeric(impsim[ic]<z[ic])
  test = data.frame("True" = z[ic], "simulated" = impsim[ic], "small" = small, "ic" = ic)
  #View(test)
  smallimpper = mean(small); smallimpper   # Simulated<True
  timp = t
}  




# -------------------------------------------------------------------------


### 2-2 figure 3-25(H) in the thesis. 

{
  # Kaplan-Meier Curve:
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("H                            t~BS(2,4), c~Exp(0.02), p=0.20, n=200") )  #KM_Estimation
  # Curve with Median of Simulated Times
  # output ----------------------------------------
  # imputation      h=hat
  bh = bugsOut$mean$b; bh
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = exp(bh[1] + bh[2]*x); lambdah #every person has specific lambda because it has specific X.
  scaleh = lambdah; scaleh
  
  #install.packages("extraDistr")
  library(extraDistr)
  # Compute median of Simulations.
  zmed = qfatigue(.5*pfatigue(cen,shapeh,scaleh, mu = 0, lower.tail = FALSE),shapeh, scaleh,mu = 0, lower.tail = FALSE)
  zmed
  # Make a variable include median of simulations.
  zimp = rep(NA,n)
  zimp[ic] = zmed[ic]
  zimp[-ic] = z[-ic]  # zimp = failure times+imputed censored times
  
  
  #
  curve2 = survfit(Surv(zimp,delta1) ~ x); curve2
  lines(curve2, mark.time = TRUE, col = "Blue", lty = 1)
  
  #Curve without Censored Times 
  tOC = z[delta==1]
  deltaOC = rep(1, length(tOC))
  curve3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); curve3
  lines(curve3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"),
         lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
  
}



# -------------------------------------------------------------------------

### 2-3  figure 3-27(B) in the thesis.

{
  # Kaplan-Meier Curve:  
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1, lwd = 2, col = "black", 
       main = paste("B                            t~BS(2,4), c~Exp(0.02), p=0.20, n=200"))  #KM_Estimation
  
  # Curve for 10,000 Times Imputation
  timp=t
  impsim = bugsOut$sims.array[,1,3+ic]  
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    
  }
  # Curve for Imputations Mean
  timp[ic] <- colMeans(impsim)
  kmmean = survfit(Surv(timp,delta1) ~ x)
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  lines(curve1, mark.time = TRUE,lty = 1, lwd = 2, col = "black", 
        main = paste("t~BS(2,4), c~Exp(0.01), p=0.10, n=200"))  #KM_Estimation
  
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
  
  
}



# -------------------------------------------------------------------------

#      3-  Breast Cancer Dataset- time as Weibull
# -------------------------------------------------------------------------

## 3-1 R codes

{
  ## Import and Make dummy Variable
  breast <- read.table("Data_Paper1.txt", header = TRUE)
  #View(breast)
  breast$Stage2 <- as.numeric(breast$Stage == 2)
  breast$Stage3 <- as.numeric(breast$Stage == 3)
  breast$Stage4 <- as.numeric(breast$Stage == 4)
  breast$Grade2 <- as.numeric(breast$Grade == 2)
  breast$Grade3 <- as.numeric(breast$Grade == 3)
  breast$Education2 <- as.numeric(breast$Education == 2)
  breast$Education3 <- as.numeric(breast$Education == 3)
  breast$Education4 <- as.numeric(breast$Education == 4)
  breast$Education5 <- as.numeric(breast$Education == 5)
  #View(breast)
  
  write.table(breast, file = "breast-Full.txt", quote = FALSE, row.names = FALSE)
  
  # Description times
  summary(breast$time[breast$Status=="1"])
  summary(breast$time[breast$Status=="0"])
  mean(breast$time[breast$Status=="1"]); sd(breast$time[breast$Status=="1"])
  mean(breast$time[breast$Status=="0"]); sd(breast$time[breast$Status=="0"])
  table(breast$Status)
  prop.table(table(breast$Status))
  table(breast$AgeC)
  
  
  #rm(list = ls())
  library(survival)
  library(R2OpenBUGS)
  # path ----------------------------------------
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Import Data
  dim(breast)
  AgeC <- breast$AgeC   #Code 1, Age<40
  prop.table(table(x))
  t <- breast$t  #time based on month
  summary(t)
  length(t)
  c <- breast$c
  length(c)
  length(t[t == "NA"])/length(t)   # Percent of Censoring, 88 Censor, 40% 
  length(c[c == "0"])/length(c)   # Percent of Observed
  n = length(breast$ID); n
  z = breast$time   # Composed from Observed and Censored data
  delta = breast$delta   # delta=0 means Censoring
  ic = which(delta == "0")   # indicator censor
  length(ic)
  
  # model
  modeltext = "model{
    for (i in 1:n) {
      t[i] ~ dweib(shape,lambda[i])C(c[i],)
    	log(lambda[i]) <- b[1] + b[2] * Stage2[i] + b[3] * Stage3[i] + b[4] * Stage4[i] +
    	             b[5] * Grade2[i] + b[6] * Grade3[i] + b[7] * Education2[i] + b[8] * Education3[i] +
    	            b[9] * Education4[i] + b[10] * Education5[i] + b[11] * ER[i] + b[12] * PR[i] +
    	            b[13] * HER2[i] + b[14] * Age[i]
    	            
    	   cim[i] <- step(c[i]-1.0E-5)*pow(log(2)/lambda[i]+pow(c[i],shape), 1/shape)      #tmed(8)
      }
    	
    	# priors
    	shape ~ dgamma(0.01,0.01)
    	for (j in 1:14) {b[j]~dnorm(0,0.01)}		
    }
    "
  
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list( n = n, t = breast$t, c = breast$c, Stage2 = breast$Stage2, Stage3 = breast$Stage3, Stage4 = breast$Stage4,
                    Grade2 = breast$Grade2, Grade3 = breast$Grade3, Education2 = breast$Education2, 
                    Education3 = breast$Education3, Education4 = breast$Education4, Education5 = breast$Education5,
                    ER = breast$ER, PR = breast$PR, HER2 = breast$HER2, Age = breast$Age )
  
  modelinit = list(list(b = rep(1,14) , shape = 1))
  param = c("shape","b","cim")
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 1
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  
  # output ----------------------------------------
  bugsOut$DIC
  
  dim(bugsOut$sims.array)   #composed: alpha, lambda, 88 simulation,deviance = 91 columns.
  
  # Head censor simulations
  bugsOut$sims.array[1:5,1,16:103]  # Head: report 1 till 5 from 100 times censored times simulations.
  # Tail censor simulations
  bugsOut$sims.array[9996:10000,1, 16:103]  # Tail:  report 1 till 5 from 100 times censored times simulations.
  bugsOut$summary[1:15, c(1:2)]    # mean & sd parameters: alpha & lambda
  
  
  parsim = bugsOut$sims.array[,1,1:15]   #parameter simulation
  head(parsim)
  
  impsim = bugsOut$sims.array[,1,16:103]  # imputation simulation
  head(impsim)
  
  dim(parsim); dim(impsim)  # 88 censored data and 15 parameter(shape+14 b)
  
  test = data.frame("True" = z[ic], "simulated" = impsim[ic])
  #View(test)
  #
  timp = t   # 88 NA and 132 failure time.
}


#------------------------------------------------------

# 3-2 figure 3-29 in the thesis.

{
  km1 = survfit(Surv(z,delta) ~ breast$AgeC); km1
  plot(km1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("Posterior Estimate: Shape=1.24,Scale=0.001,DIC=1698"))  #KM_Estimation
  
  # output ----------------------------------------
  # imputation      h=hat
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = bugsOut$mean$lambda    # 220 lambda
  mean(lambdah)
  
  
  cen=c
  #zmed = qweibull(.5*pweibull(cen,shapeh,lambdah, lower.tail = FALSE),shapeh, lambdah, lower.tail = FALSE)
  #zmed
  #length(zmed[zmed=="Inf"])
  
  library(miscTools)
  zmed = colMedians(impsim)
  zmed
  #
  ic = which(delta==0); ic      #index censor to count number of censored case.
  length(ic)
  zimp = rep(NA,n)
  zimp[ic] = zmed[ic]
  zimp[-ic] = z[-ic]  # zimp = imputed censored times(1:88) + failure times(89:220)
  delta1 = rep(1,n)  # after impute, all of times are observed then we made delta1.
  #
  km2 = survfit(Surv(zimp,delta1) ~ breast$AgeC); km2     # Bayesian Imputation
  lines(km2, mark.time = TRUE, col = "Blue", lty = 1)
  
  
  tOC = z[delta==1]  # number of observed times
  deltaOC = rep(1, length(tOC))
  length(deltaOC)
  km3 = survfit(Surv(tOC, deltaOC) ~ breast$AgeC[delta==1]); km3       # Omitting_Censored
  lines(km3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
}


#-------------------------------------------------------------------------------


#   3-3 figure 3-28 in the thesis.



{
  km1 = survfit(Surv(z,delta) ~ breast$AgeC); km1
  plot(km1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
       main = paste("t~Weibull, p=0.40, n=220"))  #KM_Estimation
  
  # simulation 
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ breast$AgeC)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    #Sys.sleep(.5)
  }
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  lines(km1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
        main = paste("t~Weibull, p=0.40, n=220"))  #KM_Estimation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
  
}  






# -------------------------------------------------------------------------
#     4- Breast Cancer Data distributed as Birnbaum-Saunders
# -------------------------------------------------------------------------  

## 4-1 R codes

{
  
  ## Import and Make dummy Variable
  breast <- read.table("Data_Paper1.txt", header = TRUE)
  #View(breast)
  breast$Stage2 <- as.numeric(breast$Stage == 2)
  breast$Stage3 <- as.numeric(breast$Stage == 3)
  breast$Stage4 <- as.numeric(breast$Stage == 4)
  breast$Grade2 <- as.numeric(breast$Grade == 2)
  breast$Grade3 <- as.numeric(breast$Grade == 3)
  breast$Education2 <- as.numeric(breast$Education == 2)
  breast$Education3 <- as.numeric(breast$Education == 3)
  breast$Education4 <- as.numeric(breast$Education == 4)
  breast$Education5 <- as.numeric(breast$Education == 5)
  #View(breast)
  
  #rm(list = ls())
  library(survival)
  library(R2OpenBUGS)
  # path ----------------------------------------
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Import Data
  dim(breast)
  x <- breast$AgeC   #Code 1, Age<40
  prop.table(table(x))
  t <- breast$t  #time based on month
  summary(t)
  length(t)
  c <- breast$c
  length(c)
  length(t[t == "NA"])/length(t)   # Percent of Censoring, 88 Censor, 40% 
  length(c[c == "0"])/length(c)   # Percent of Observed
  n = length(x); n
  z = breast$z   # Composed from Observed and Censored data
  delta = breast$delta   # delta=0 means Censoring
  ic = which(delta == 0)   # indicator censor
  length(ic)
  age <- breast$AgeC
  # model
  modeltext = "model {
          for(i in 1:n){
      	  t[i] ~ dbs(shape, lambda[i])C(c[i], )
      		log(lambda[i]) <- b[1] + b[2] * Stage2[i] + b[3] * Stage3[i] + b[4] * Stage4[i] +
        	                  b[5] * Grade2[i] + b[6] * Grade3[i] + b[7] * Education2[i] + b[8] * Education3[i] +
        	                  b[9] * Education4[i] + b[10] * Education5[i] + b[11] * ER[i] + b[12] * PR[i] +
        	                  b[13] * HER2[i] + b[14] * Age[i]
        # conditional tmed
        cim[i] <- step(c[i]-1.0E-5) * lambda[i] * pow( (shape/2) * 0.5875441*logit(0.5 + 0.5 * phi((1/shape) * 
                  (sqrt(c[i]/lambda[i]) - sqrt(lambda[i]/c[i])))) + sqrt( 1 + pow( (shape/2) * 
                  0.5875441 * logit(0.5 + 0.5 * phi((1/shape) * (sqrt(c[i]/lambda[i]) -
                  sqrt(lambda[i]/c[i])))),2)) ,2)
        
  }
	# priors
	shape ~ dgamma(0.01,0.01)
	for(j in 1:14) {b[j]~dnorm(0,0.01)}		
 }
 "

  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list( n = n, Stage2 = breast$Stage2, Stage3 = breast$Stage3, Stage4 = breast$Stage4,
                    Grade2 = breast$Grade2, Grade3 = breast$Grade3, Education2 = breast$Education2, 
                    Education3 = breast$Education3, Education4 = breast$Education4, Education5 = breast$Education5,
                    ER = breast$ER, PR = breast$PR, HER2 = breast$HER2, Age = breast$Age , t = breast$t, c = breast$c)
  modelinit = list(list(b = rep(1,14) , shape = 1))
  param = c("shape","b","cim", "lambda")
  
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 1
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  
  # output ----------------------------------------
  bugsOut$DIC
  dim(bugsOut$sims.array)   #composed: alpha, lambda, 88 simulation,deviance = 91 columns.
  bugsOut$sims.array[1:5,1,16:103]  # report 1 till 5 from 88 times censored times simulations.
  bugsOut$summary[1:15, c(1:2)]    # mean & sd parameters: alpha & lambda
  
  
  parsim = bugsOut$sims.array[,1,1:15]   #parameter simulation: 100*2
  impsim = bugsOut$sims.array[,1,16:103]  # imputation simulation: 100*88
  
  test = data.frame("True" = z[ic], "simulated" = impsim[ic])
  #View(test)
  timp = t
}


#------------------------------------------------------


# 4-2 figure 3-31 in the thesis.

{
  km1 = survfit(Surv(z,delta) ~ breast$AgeC); km1
  plot(km1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("Posterior Estimate: Shape=1.22, Scale=145.21, DIC=1510"))  #KM_Estimation
  
  # output ----------------------------------------
  # imputation      h=hat
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = bugsOut$mean$lambda; lambdah
  
  cen=c
  # Corredted -->
  zmed = qweibull(.5*pweibull(cen,shapeh,lambdah, lower.tail = FALSE),shapeh, lambdah, lower.tail = FALSE)
  zmed
  
  
  #
  ic = which(delta==0); ic      #index censor to count number of censored case.
  length(ic)
  zimp <- rep(NA, n)
  zimp[ic] <- zmed[ic]
  zimp[-ic] <- z[-ic]  # zimp = failure times+imputed censored times
  zimp
  delta1 = rep(1,n)  # after impute, all of times are observed then we made delta1.
  #
  km2 = survfit(Surv(zimp,delta1) ~ x); km2     # Bayesian Imputation
  lines(km2, mark.time = TRUE, col = "Blue", lty = 1)
  
  
  tOC = z[delta==1]  # number of observed times
  deltaOC = rep(1, length(tOC))
  length(deltaOC)
  km3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); km3       # Omitting_Censored
  lines(km3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
}



#-------------------------------------------------------------------------------

## 4-3 figure 3-30 in the thesis.


{
  km1 = survfit(Surv(z,delta) ~ breast$AgeC); km1
  plot(km1, mark.time = TRUE,lty = 1, lwd=2, col = "black",
       main = paste("t~Birnbaum-Saunders, p=0.40, n=220"))  #KM_Estimation
  
  # simulation 
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    #Sys.sleep(.5)
  }
  timp[ic] <- colMeans(impsim)
  kmmean = survfit(Surv(timp,delta1) ~ x)
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  lines(km1, mark.time = TRUE,lty = 1, lwd=2, col = "black",
        main = paste("t~Birnbaum-Saunders, p=0.40, n=220"))  #KM_Estimation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
}  



# -------------------------------------------------------------------------
#  5- Mix Bayesian Network, package deal.
# -------------------------------------------------------------------------  


# 5-1 Simulate censored times in the breast cancer data to impute them.

{
  
  ## Import and Make dummy Variable
  breast <- read.table("Data_Paper1.txt", header = TRUE)
  #View(breast)
  breast$Stage2 <- as.numeric(breast$Stage == 2)
  breast$Stage3 <- as.numeric(breast$Stage == 3)
  breast$Stage4 <- as.numeric(breast$Stage == 4)
  breast$Grade2 <- as.numeric(breast$Grade == 2)
  breast$Grade3 <- as.numeric(breast$Grade == 3)
  breast$Education2 <- as.numeric(breast$Education == 2)
  breast$Education3 <- as.numeric(breast$Education == 3)
  breast$Education4 <- as.numeric(breast$Education == 4)
  breast$Education5 <- as.numeric(breast$Education == 5)
  #View(breast)
  
  #rm(list = ls())
  library(survival)
  library(R2OpenBUGS)
  # path ----------------------------------------
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Import Data
  dim(breast)
  x <- breast$AgeC   #Code 1, Age<40
  prop.table(table(x))
  t <- breast$t  #time based on month
  summary(t)
  length(t)
  c <- breast$c
  length(c)
  length(t[t == "NA"])/length(t)   # Percent of Censoring, 88 Censor, 40% 
  length(c[c == "0"])/length(c)   # Percent of Observed
  n = length(x); n
  z = breast$z   # Composed from Observed and Censored data
  delta = breast$delta   # delta=0 means Censoring
  ic = which(delta == 0)   # indicator censor
  length(ic)
  age <- breast$AgeC
  # model
  modeltext = "model {
          for(i in 1:n){
      	  t[i] ~ dbs(shape, lambda[i])C(c[i], )
      		log(lambda[i]) <- b[1] + b[2] * Stage2[i] + b[3] * Stage3[i] + b[4] * Stage4[i] +
        	                  b[5] * Grade2[i] + b[6] * Grade3[i] + b[7] * Education2[i] + b[8] * Education3[i] +
        	                  b[9] * Education4[i] + b[10] * Education5[i] + b[11] * ER[i] + b[12] * PR[i] +
        	                  b[13] * HER2[i] + b[14] * Age[i]
        
        cim[i] <- step(c[i]-1.0E-5) * lambda[i] * pow( (shape/2) * 0.5875441*logit(0.5 + 0.5 * phi((1/shape) * 
                  (sqrt(c[i]/lambda[i]) - sqrt(lambda[i]/c[i])))) + sqrt( 1 + pow( (shape/2) * 
                  0.5875441 * logit(0.5 + 0.5 * phi((1/shape) * (sqrt(c[i]/lambda[i]) -
                  sqrt(lambda[i]/c[i])))),2)) ,2)
        
  }
	# priors
	shape ~ dgamma(0.01,0.01)
	for(j in 1:14) {b[j]~dnorm(0,0.01)}		
  }
  "
  
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list( n = n, Stage2 = breast$Stage2, Stage3 = breast$Stage3, Stage4 = breast$Stage4,
                    Grade2 = breast$Grade2, Grade3 = breast$Grade3, Education2 = breast$Education2, 
                    Education3 = breast$Education3, Education4 = breast$Education4, Education5 = breast$Education5,
                    ER = breast$ER, PR = breast$PR, HER2 = breast$HER2, Age = breast$Age , t = breast$t, c = breast$c)
  modelinit = list(list(b = rep(1,14) , shape = 1))
  param = c("shape","b","cim", "lambda")
  
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 1
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  
  # output ----------------------------------------
  bugsOut$DIC
  dim(bugsOut$sims.array)   #composed: alpha, lambda, 88 simulation,deviance = 91 columns.
  bugsOut$sims.array[1:5,1,16:103]  # report 1 till 5 from 88 times censored times simulations.
  bugsOut$summary[1:15, c(1:2)]    # mean & sd parameters: alpha & lambda
  
  
  parsim = bugsOut$sims.array[,1,1:15]   #parameter simulation: 100*2
  impsim = bugsOut$sims.array[,1,16:103]  # imputation simulation: 100*88
  
  test = data.frame("True" = z[ic], "simulated" = impsim[ic])
  #View(test)
  timp = t
  
  #imputation      h=hat
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = bugsOut$mean$lambda; lambdah
  lambdamed = bugsOut$median$lambda; lambdamed
  
  cen=c
  zmed = qweibull(.5*pweibull(cen,shapeh,lambdah, lower.tail = FALSE),shapeh, lambdah, lower.tail = FALSE)
  zmed
  #View(breast)
  #Incorporate failure time with imputed censor time
  ic = which(delta==0); ic      #index censor to count number of censored case.
  length(ic)
  zimp = rep(NA,n)
  zimp[ic] = zmed[ic]
  zimp[-ic] = z[-ic]  
  zimp   # zimp = failure times+imputed censored times#View(breast)
}


# 5-2 Run Bayesian Network and DAG.


{
  # Structural Learning, Time -----> u = 0.5 (sqrt(t/beta) - sqrt(beta/t))
  
  # Structural Learning
  #Add covariates from Breast Cancer data.

  BStrans <- 0.5 * (sqrt(zimp/lambdah) - sqrt(lambdah/zimp))
  
  
  datBS <- cbind(Status = breast$Status, Time = BStrans, Age = breast$Age, Education = breast$Education, 
                 ER = breast$ER, PR = breast$PR, Her2 = breast$HER2 , Stage = breast$Stage, Grade = breast$Grade
  )
  datBS <- as.data.frame(datBS)
  #View(datBS)
  write.csv(datBS, file="datBS.csv") 
  
  # Multi-Collinearity
  model1 <- lm(zimp~ Age + Education + Stage + ER + PR + Her2 + Grade, data = datBS)
  #install.packages("olsrr")
  library(olsrr)  
  
  
  ols_vif_tol(model1)
  
  datBS$Status <- as.factor(datBS$Status)
  datBS$Stage <- as.factor(datBS$Stage)
  datBS$Grade <- as.factor(datBS$Grade)
  datBS$ER <- as.factor(datBS$ER)
  datBS$PR <- as.factor(datBS$PR)
  datBS$Her2 <- as.factor(datBS$Her2)
  datBS$Education <- as.factor(datBS$Education)
  datBS$Age <-as.numeric(datBS$Age)
  datBS$Time<- as.numeric(datBS$Time)
  #View(datBS)
  summary(datBS$Time); sd(datBS$Time)
  cor.test(datBS$Age, datBS$Time)
  reg1 <- lm(Time ~ Age, data = datBS)
  summary(reg1)
  plot(datBS$Age, datBS$Time, ylab = "Time to Death", xlab = "Age at Diagnosis cancer")
  text(65,800,"Correlation = -0.13, Sig. = 0.056")
  abline(lm(datBS$Time~datBS$Age))
  
  
  #                         Comment Deal: 
  #   a data frame, where the columns define the variables. 
  #   A continuous variable should have type numeric and 
  #   discrete variables should have type factor.
  
  # Fit Bayesian Network
  #
  #setwd("C:/Users/novingostar/Documents/R-studio")
  #install.packages("deal")
  library(deal)
  ## specify prior network
  
  # We have no prior knowledge about specific dependency relations, 
  # so for simplicity we use the empty DAG as the prior DAG and 
  # let the probability distribution of the discrete variables be uniform.
  
  ## specify prior network
  br.nw <- network(datBS)
  plot(br.nw)    # an empty DAG
  # Number of networks can be draw for nd and nc. 
  numbermixed(7,2)    
  # Adjust prior probability distribution br.nw:
  br.nw$nodes$Education$prob <- c(0.22, 0.22, 0.20, 0.20, 0.16)
  br.nw$nodes$ER$prob <- c(0.55, 0.45)
  br.nw$nodes$PR$prob <- c(0.55, 0.45)
  br.nw$nodes$Her2$prob <- c(0.45, 0.55)
  br.nw$nodes$Stage$prob <- c(0.20, 0.23, 0.27 ,0.30 )
  br.nw$nodes$Grade$prob <- c(0.25, 0.30, 0.45)
  
  # Non-adjust probability distribution br.nw:
  # br.nw$nodes$Education$prob <- c(0.20, 0.20, 0.20, 0.20, 0.20)
  # br.nw$nodes$ER$prob <- c(0.50, 0.50)
  # br.nw$nodes$PR$prob <- c(0.50, 0.50)
  # br.nw$nodes$Her2$prob <- c(0.50, 0.50)
  # br.nw$nodes$Stage$prob <- c(0.25, 0.25, 0.25, 0.25 )
  # br.nw$nodes$Grade$prob <- c(0.33, 0.33, 0.33)
  
  
  ## make joint prior distribution
  br.prior <- jointprior(br.nw, phiprior = "bottcher")
  br.prior
  # We do not allow arrows into Age and ER and set ban list.
  ## ban arrows towards Age and ER.
  
  mybanlist <- matrix(c(1,1,1,1,1,1,1,1,  2,2,2,2,2,2,2,2,  3,3,3,3,  4,4,4,4,4,  5,5,5,5,  6,6,6,6,  7,7,7,7,  8,8,8,8,8, 9,9,9,9,9,  # *-->
                        2,3,4,5,6,7,8,9,  1,3,4,5,6,7,8,9,  5,6,7,9,  3,5,6,7,9,  3,4,6,7,  3,4,5,7,  3,4,5,6,  3,4,5,6,7, 3,4,5,6,7), ncol=2)
  
  
  banlist(br.nw) <- mybanlist
  #mybanlist                     # banned edges indicated with red line.
  
  ## learn the initial network
  
  br.nw <- learn(br.nw, datBS, br.prior)$nw
  plot(br.nw)   # show banlist   
  
  
  ## Do structural search
  br.search <- autosearch(initnw = br.nw, data = datBS, prior = br.prior,
                          maxiter = 1000, trace=TRUE, removecycles = TRUE)
  
  
  # DAG to Paper ----->>
  plot.network(thebest1 <- br.search$nw, showban = FALSE, cexscale =11 ,arrowlength =0.20 , 
               font = 1,unitscale = 30,
               yr=c(350,0), xr = c(0,330))
}


# -------------------------------------------------------------------------
#               The END
# -------------------------------------------------------------------------



