library(sas7bdat)
setwd("C:/Users/SoYoung Kim/Downloads/test1")
OriginalData <-read.sas7bdat("ml_data_chap3.sas7bdat",debug=FALSE)

####################
# Extract residual #
####################
library(lme4)
fm1 <- lmer (lninj ~ 1+ HIV_Region + ethnic + HIV_Region*ethnic 
             + gender + gdmc_age + Edu_gp1 + Edu_gp2 + (ethnic | site), 
             OriginalData, REML=TRUE) # Analyze with the original data
Level2.resd <- ranef(fm1)$site # Extract level-2 residual
names(Level2.resd)<- c("u0", "u1")
Level2.resd$site <- as.numeric(row.names(Level2.resd))
Comp.resd <- as.numeric(resid(fm1)) # Extract composite residual
library(dplyr)
tmp <- full_join(OriginalData, Level2.resd, by="site") 
tmp$u_j <- tmp$u0 + (tmp$u1*tmp$ethnic)
tmp$e_ij <- Comp.resd - tmp$u_j # Calculate level-1 residual 


# Transform Level-2 residual 

# matrix "S" is the empirical var/cov matrix 
center_U <- Level2.resd
center_U$c_u0 <- Level2.resd$u0 - mean(Level2.resd$u0)
center_U$c_u1 <- Level2.resd$u1 - mean(Level2.resd$u1)
C_U <- center_U %>% select(c_u0, c_u1)
CU <- as.matrix(as.data.frame(lapply(C_U, as.numeric)))
SbyJ <- t(CU) %*% CU
J <- nrow (C_U)
S2 <- SbyJ/J
LS2 <- t(chol(S2))

# "Sigma" is the model estimated var/cov matrix
fm1tmp <- summary(fm1)
Sigma2 <- matrix(as.numeric(fm1tmp$varcor$site),2,byrow=T)
LSigma2 <- t(chol(Sigma2))
library(matlib)
A2 <- t(LSigma2%*%inv(LS2))
New_U <- CU%*%A2
tmp_U <- as.numeric(New_U)
New_U <- matrix(tmp_U,nrow(Level2.resd))
colnames(New_U) <- c("N_u0","N_u1")


# Transform Level-1 residual 

# matrix "S" is the empirical var/cov matrix 
tmp$c_eij <- tmp$e_ij - mean(tmp$e_ij)
C_E <- tmp %>% select(c_eij)
CE <- as.matrix(as.data.frame(lapply(C_E, as.numeric)))
SbyI <- t(CE) %*% CE
I <- nrow (C_E)
S1 <-SbyI/I 
LS1 <- chol(S1)

# "Sigma" is the model estimated var/cov matrix
Sigma1 <- getME(fm1,"sigma")^2
LSigma1 <- chol(Sigma1)
A1 <- LSigma1/LS1
New_E <- CE%*%A1



# Set parameters

# I = 9824; Sample size 
# J = 20; Group size
G00_hat <- fixef(fm1)[1] # effect of Intercept  =  3.59552039
G01_hat <- fixef(fm1)[2] # effect of HIV_Region =  0.53096931
G10_hat <- fixef(fm1)[3] # effect of ethnic     = -0.54761327
G11_hat <- fixef(fm1)[8] # effect of HIV_Region:ethnic = 0.33108555
B2_hat <- fixef(fm1)[4] # effect of gender     = -0.04521547
B3_hat <- fixef(fm1)[5] # effect of gdmc_age   =  0.01395336
B4_hat <- fixef(fm1)[6] # effect of Edu_gp1    = -0.08517052
B5_hat <- fixef(fm1)[7] # effect of Edu_gp2    = -0.04973010


##################
# Bootstrap Loop #
##################
nboot <- 500
for(i in 1:nboot) {
  SRS.E <- sample(New_E, I, replace=TRUE)
  x <-sample(nrow(New_U),J,replace=T)
  SRS.U <- New_U[x,]
  fm1frame<-fm1@frame 
  level1Data<-fm1frame[,c('site', 'lninj', 'ethnic', 'gender', 'gdmc_age', 'Edu_gp1', 'Edu_gp2')]
  level2Data<-fm1frame[!duplicated(fm1frame$site),c('site','HIV_Region')]
  Newlevel1Data <- cbind(level1Data, SRS.E)
  Newlevel2Data <- cbind(level2Data, SRS.U)
  Newlevel2Data$New_alpha = G00_hat + G01_hat*Newlevel2Data$HIV_Region + Newlevel2Data$N_u0 
  Newlevel2Data$New_beta1 = G10_hat + G11_hat*Newlevel2Data$HIV_Region + Newlevel2Data$N_u1  
  BootstrapData <- merge(Newlevel1Data, Newlevel2Data, "site")
  head(BootstrapData)
  class(BootstrapData)
  BootstrapData$New_y <- BootstrapData$New_alpha +
    BootstrapData$New_beta1*BootstrapData$ethnic + 
    B2_hat*BootstrapData$gender +
    B3_hat*BootstrapData$gdmc_age +
    B4_hat*BootstrapData$Edu_gp1 +
    B5_hat*BootstrapData$Edu_gp2 +
    BootstrapData$SRS.E 
  
  doboot <- lmer (New_y ~ 1+ HIV_Region + ethnic + HIV_Region*ethnic +
                    gender + gdmc_age + Edu_gp1 + Edu_gp2 + (ethnic | site), 
                  BootstrapData, REML=TRUE)
  smr.doboot <- summary(doboot)
  Reg.coef <- fixef(doboot)
  u_11 <- smr.doboot$varcor$site[1]
  u_21 <- smr.doboot$varcor$site[3]
  u_22 <- smr.doboot$varcor$site[4]
  sigma.sqr <- smr.doboot$sigma^2
  Resd_Cov <- c(u_11,u_21,u_22,sigma.sqr)
  names(Resd_Cov) <- c("u_11","u_21","u_22","sigma.sqr")
  
  if (doboot@beta[1] == summary(doboot)$coefficient[1]) {
    Fixed_p <- as.data.frame( t(data.frame(Reg.coef)) )
    Random_p <- as.data.frame( t(data.frame(Resd_Cov)) )
    if (i==1) {
      fixed_effects <- Fixed_p
    }  else {
      fixed_effects <- bind_rows (fixed_effects, Fixed_p)}
    if (i==1) {
      random_effects<-Random_p
    }  else {
      random_effects <- bind_rows(random_effects, Random_p)}
  }
}

# Estimates Result

Est.fixed <- c(
  mean(fixed_effects$`(Intercept)`), 
  mean(fixed_effects$HIV_Region),
  mean(fixed_effects$ethnic),
  mean(fixed_effects$gender),
  mean(fixed_effects$gdmc_age),
  mean(fixed_effects$Edu_gp1),
  mean(fixed_effects$Edu_gp2),
  mean(fixed_effects$`HIV_Region:ethnic`)
)
names(Est.fixed) <- names(fixed_effects)
boot_fixed <- data.frame(Est.fixed)

Est.random <- c(
  mean(random_effects$u_11),
  mean(random_effects$u_21),
  mean(random_effects$u_22),
  mean(random_effects$sigma.sqr)
)
names(Est.random) <- names(random_effects)
boot_random <- data.frame(Est.random)

boot_fixed
boot_random


