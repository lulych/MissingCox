#Load the packages to use
paquetes <- c("readxl", "dplyr", "survival", "hydroGOF", "broom", "VIM", 
              "missForest", "mice", "coxed")
lapply(paquetes, require, character.only = TRUE)

#-------------------------------------------------------------------------------
#General configuration
set.seed(127)
#(100,0.10)=276, (100,0.20)=123, (100,0.30)=127, (100,0.40)=No converge
#(200,0.10)=276, (200,0.20)=123, (200,0.30)=127, (200,0.40)=No converge
#(500,0.10)=276, (500,0.20)=123, (500,0.30)=127, (500,0.40)=777
#(1000,0.10)=276, (1000,0.20)=123, (1000,0.30)=127, (1000,0.40)=775
n <- 500 #sample size
sim <- 10000 #Number of iterations
prop <- 0.3 #Missing data proportion

cuali <- c(2,5) #Vector of qualitative variable indices
cuanti <- c(1,3,4) #Vector of indices of quantitative variables

#-------------------------------------------------------------------------------
#Prepare arrays to store simulation results
estim <- array(0,c(sim,6,12))
InfCI <- array(0,c(sim,6,12))
SupCI <- array(0,c(sim,6,12))
nrmse <- array(0,c(sim,length(cuanti),12))
MICE_nrmse <- array(0,c(5,length(cuanti),1))
meannrmse <- array(0,c(12,length(cuanti)))
ecm <- array(0,c(12,6))
cobertura <- array(0,c(12,6))
MICE_FCE <- array(0,c(5,length(cuali),1))
FCE <- array(0,c(sim,length(cuali),12))
meanFCE <- array(0,c(12,length(cuali)))
times <- array(NA,c(sim,12))
meantimes <- array(NA,c(1,12))

par <- c("B1","B2","B3","B4","B5","B6")
metodo <- c("CCA","RVI","K3-","k3+","K5-","k5+","k9-","k9+",
            "MF-","MF+","MICE-","MICE+")

dimnames(estim) <- list(c(1:sim),par,metodo)
dimnames(InfCI) <- list(c(1:sim),par,metodo)
dimnames(SupCI) <- list(c(1:sim),par,metodo)
dimnames(nrmse) <- list(c(1:sim),cuanti,metodo)
dimnames(MICE_nrmse) <- list(c(1:5),cuanti)
dimnames(meannrmse) <- list(metodo,cuanti)
dimnames(ecm) <- list(metodo,par)
dimnames(cobertura) <- list(metodo,par)
dimnames(MICE_FCE) <- list(c(1:5),cuali)
dimnames(FCE) <- list(c(1:sim),cuali, metodo)
dimnames(meanFCE) <- list(metodo, cuali)

colnames(times) <- metodo
colnames(meantimes) <- metodo

#-------------------------------------------------------------------------------
#Simulate data for the covariates and the times and states for prefixed beta's 
for (s in 1:sim){
X1 <- rnorm(n = n, mean = 0, sd = 1)
print (c(s,"X1 ok"))
X2 <- rbinom(n=n, size=1, prob=0.3)
print (c(s,"X2 ok"))
X3 <-  runif(n,0,10)
print (c(s,"X3 ok"))
X4 <- rexp(n, rate = 1)
print (c(s,"X4 ok"))
y <- c(1,2,3)
X5 <- colSums(rmultinom(prob=c(0.3,0.2,0.5),size=1, n=n)*y)
print (c(s,"X5 ok"))
Xa <- ifelse(X5==2,1,0)
print (c(s,"Xa ok"))
Xb <- ifelse(X5==3,1,0)
print (c(s,"Xb ok"))
simdata <- cbind(X1,X2,X3,X4,Xa,Xb)

my.hazard <- function(t){ 
  dnorm((log(t) - log(50))/log(10)) /
    (log(10)*t*(1 - pnorm((log(t) - log(50))/log(10))))
}

simdata2 <- sim.survdata(T=1000, X=simdata, num.data.frames = 1, censor=0.2,
                         beta=c(-0.10,0.25,-0.13,0.50,0.35,-0.15), hazard.fun=my.hazard)
v <- simdata2$data
print (c(s,"data sim ok"))

coef <- c(-0.10,0.25,-0.13,0.50,0.35,-0.15)
coefcorr <- matrix(t(coef),nrow=sim,ncol=6,byrow=T)

for (j in 1:n){
  ifelse(simdata2$data$Xa[j]==1, simdata2$data$X5[j] <- 2, 
         ifelse (simdata2$data$Xb[j]==1, simdata2$data$X5[j] <- 3, simdata2$data$X5[j] <- 1))
}
print (c(s,"simdata2 ok"))
cox <- coxph(Surv(y, failed) ~ X1+X2+X3+X4+Xa+Xb, data = simdata2$data, method="breslow")

data <- as.data.frame(cbind(simdata2$data[,c(1:4,9,7:8)]))
data$X2 <- as.factor(data$X2)
data$X5 <- as.factor(data$X5)

data2 <- data #Database with missing
miss <- prop*n #No. of data to delete per variable
miss <- round(miss,0)

#-------------------------------------------------------------------------------
#Missing data generation (NA's)
for (k in 1:5){ #For each variable
  i <- 1
  while (i < miss+1) { #Repeat until you reach no. of missing in each variable
    j <- sample(1:n, 1) #Pick the IDs randomly
    if (is.na(data2[j,k]) == 'FALSE') {data2[j,k] <- NA #Check that it is not NA
    i <- i+1}  #Keeps track of NA generated
  }}
print (c(s,"missing data ok"))

#-------------------------------------------------------------------------------
#Estimate under CCA
cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = data2, method="breslow")

resul<- tidy(cox)
estim[s,,1] <- t(resul[,2])
InfCI[s,,1] <- t(resul[,6])
SupCI[s,,1] <- t(resul[,7])

print (c(s,"CCA ok"))

#-------------------------------------------------------------------------------
#Representative value imputation 
impu <- data2
t <- proc.time()

#Median imputation
impu$X1[(is.na(impu$X1))] <- median(impu$X1[!is.na(impu$X1)])
impu$X3[(is.na(impu$X3))]  <- median(impu$X1[!is.na(impu$X3)])
impu$X4[(is.na(impu$X4))]  <- median(impu$X1[!is.na(impu$X4)])

#Mode Imputation
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

v <- subset(impu$X2, impu$X2!='NA')
modo <- getmode(v)
impu$X2[(is.na(impu$X2))] <- modo

v <- subset(impu$X5, impu$X5!='NA')
modo <- getmode(v)
impu$X5[(is.na(impu$X5))] <- modo

tiempo <- proc.time() - t
times[s,2] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = impu, method="breslow")

resul <- tidy(cox)
estim[s,,2] <- t(resul[,2])
InfCI[s,,2] <- t(resul[,6])
SupCI[s,,2] <- t(resul[,7])

nrmse[s,,2] <- hydroGOF::nrmse(impu[,cuanti], data[,cuanti])

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,2] <- ifelse(impu[j,i]== data[j,i],FCE[s,r,2],FCE[s,r,2]+1)
  }
  r <- r+1}
print (c(s,"RVI ok"))

#-------------------------------------------------------------------------------
#K-nearest neighbor imputation

# Option (-): Only covariates are considered in the calculation of the distance
# Option (+): the time of censorship or event and the Status are also considered in the calculation of the distance 

#For k=3

#Option (-)
t <- proc.time()
k3a <- kNN(data2, dist_var=c("X1","X2","X3","X4","X5"), k=3)
tiempo <- proc.time() - t
times[s,3] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k3a, method="breslow")

resul <- tidy(cox)

estim[s,,3] <- t(resul[,2])
InfCI[s,,3] <- t(resul[,6])
SupCI[s,,3] <- t(resul[,7])

nrmse[s,,3] <- hydroGOF::nrmse(k3a[,cuanti], data[,cuanti])
k3a <- k3a[,c(1:7)]

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,3] <- ifelse(k3a[j,i]== data[j,i],FCE[s,r,3],FCE[s,r,3]+1)
  }
  r <- r+1}
print (c(s,"K3- ok"))

#Option (+)
t <- proc.time()
k3b <- kNN(data2, k=3)
tiempo <- proc.time() - t
times[s,4] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k3b, method="breslow")

resul <- tidy(cox)
estim[s,,4] <- t(resul[,2])
InfCI[s,,4] <- t(resul[,6])
SupCI[s,,4] <- t(resul[,7])

nrmse[s,,4] <- hydroGOF::nrmse(k3b[,cuanti], data[,cuanti])

k3b <- k3b[,c(1:7)]
r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,4] <- ifelse(k3b[j,i]== data[j,i],FCE[s,r,4],FCE[s,r,4]+1)
  }
  r <- r+1}
print (c(s,"k3+ ok"))


#For k=5

#Option (-)
t <- proc.time()
k5a <- kNN(data2, dist_var=c("X1","X2","X3","X4","X5"), k=5)
tiempo <- proc.time() - t
times[s,5] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k5a, method="breslow")

resul <- tidy(cox)

estim[s,,5] <- t(resul[,2])
InfCI[s,,5] <- t(resul[,6])
SupCI[s,,5] <- t(resul[,7])

nrmse[s,,5] <- hydroGOF::nrmse(k5a[,cuanti], data[,cuanti])
k5a <- k5a[,c(1:7)]

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,5] <- ifelse(k5a[j,i]== data[j,i],FCE[s,r,5],FCE[s,r,5]+1)
  }
  r <- r+1}
print (c(s,"k5- ok"))

#Option (+)
t <- proc.time()
k5b <- kNN(data2, k=5)
tiempo <- proc.time() - t
times[s,6] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k5b, method="breslow")

resul <- tidy(cox)
estim[s,,6] <- t(resul[,2])
InfCI[s,,6] <- t(resul[,6])
SupCI[s,,6] <- t(resul[,7])

nrmse[s,,6] <- hydroGOF::nrmse(k5b[,cuanti], data[,cuanti])

k5b <- k5b[,c(1:7)]
r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,6] <- ifelse(k5b[j,i]== data[j,i],FCE[s,r,6],FCE[s,r,6]+1)
  }
  r <- r+1}
print (c(s,"k5+ ok"))

#For k=9

#Option (-)
t <- proc.time()
k9a <- kNN(data2, dist_var=c("X1","X2","X3","X4","X5"), k=9)
tiempo <- proc.time() - t
times[s,7] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k9a, method="breslow")

resul <- tidy(cox)

estim[s,,7] <- t(resul[,2])
InfCI[s,,7] <- t(resul[,6])
SupCI[s,,7] <- t(resul[,7])

nrmse[s,,7] <- hydroGOF::nrmse(k3a[,cuanti], data[,cuanti])
k9a <- k9a[,c(1:7)]

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,7] <- ifelse(k9a[j,i]== data[j,i],FCE[s,r,7],FCE[s,r,7]+1)
  }
  r <- r+1}
print (c(s,"k9- ok"))

#Option (+)
t <- proc.time()
k9b <- kNN(data2, k=9)
tiempo <- proc.time() - t
times[s,8] <- tiempo[3]

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = k9b, method="breslow")

resul <- tidy(cox)
estim[s,,8] <- t(resul[,2])
InfCI[s,,8] <- t(resul[,6])
SupCI[s,,8] <- t(resul[,7])

nrmse[s,,8] <- hydroGOF::nrmse(k9b[,cuanti], data[,cuanti])

k9b <- k9b[,c(1:7)]
r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,8] <- ifelse(k9b[j,i]== data[j,i],FCE[s,r,8],FCE[s,r,8]+1)
  }
  r <- r+1}
print (c(s,"k9+ ok"))

#----------------------------------------------------------------------------
#MissForest imputation

#Option (-)
t <- proc.time()
mf <- missForest(data2[,c(1:5)]) #To impute without considering status and days
tiempo <- proc.time() - t
times[s,9] <- tiempo[3]

#Join the imputed data with status and days
mf <- cbind.data.frame(mf$ximp,data2[,c(6,7)]) 

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = mf, method="breslow")

resul <- tidy(cox)
estim[s,,9] <- t(resul[,2])
InfCI[s,,9] <- t(resul[,6])
SupCI[s,,9] <- t(resul[,7])

nrmse[s,,9] <- hydroGOF::nrmse(mf[,cuanti], data[,cuanti])

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,9] <- ifelse(mf[j,i]== data[j,i],FCE[s,r,9],FCE[s,r,9]+1)
  }
  r <- r+1}
print (c(s,"Mf- ok"))

#Option (+)
t <- proc.time()
mfb <- missForest(data2)
tiempo <- proc.time() - t
times[s,10] <- tiempo[3]

mfb <- mfb$ximp #Keeps the imputed data and discards other results

cox <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = mfb, method="breslow")

resul <- tidy(cox)
estim[s,,10] <- t(resul[,2])
InfCI[s,,10] <- t(resul[,6])
SupCI[s,,10] <- t(resul[,7])

nrmse[s,,10] <- hydroGOF::nrmse(mfb[,cuanti], data[,cuanti])

r <- 1
for (i in cuali){
  for (j in 1:n){
    FCE[s,r,10] <- ifelse(mfb[j,i]== data[j,i],FCE[s,r,10],FCE[s,r,10]+1)
  }
  r <- r+1}
print (c(s,"Mf+ ok"))

#-------------------------------------------------------------------------------
#MICE imputation

#Option (-)
t <- proc.time()
impua <- mice(data2[,c(1:5)]) #To impute without considering status and days
tiempo <- proc.time() - t
times[s,11] <- tiempo[3]

analyses <- as.list(1:impua$m)
for (nim in 1:impua$m){
  compl <- cbind.data.frame(complete(impua,nim),data2[,c(6,7)])
  
  #Cox regression analysis for each imputed set
  analyses[[nim]] <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = compl, method="breslow")
  r <- 1
  for (i in cuali){
    for (j in 1:n){
      MICE_FCE[nim,r,] <- ifelse(compl[j,i]== data[j,i],MICE_FCE[nim,r,],MICE_FCE[nim,r,]+1)
    }
    r <- r+1}
  MICE_nrmse[nim,,] <- hydroGOF::nrmse(compl[,cuanti], data[,cuanti])
}
FCE[s,,11] <- colMeans(MICE_FCE)
MICE_FCE <- array(0,c(5,length(cuali),1))
nrmse[s,,11] <- colMeans(MICE_nrmse)
MICE_nrmse <- array(0,c(5,length(cuanti),1))
object <- list(call=call, call1=impua$call, nmis=impua$nmis, analyses=analyses)
oldClass(object) <- c("mira", "coxph")
coxmicea <- summary(pool(object))

estim[s,,11] <- t(coxmicea[,1])
InfCI[s,,11] <- t(coxmicea[,6])
SupCI[s,,11] <- t(coxmicea[,7])
print (c(s,"MICE- ok"))

#Option (+)
t <- proc.time()
impub <- mice(data2) #To impute considering status and days
tiempo <- proc.time() - t
times[s,12] <- tiempo[3]

analyses <- as.list(1:impub$m)
for (nim in 1:impub$m){
  compl <- cbind.data.frame(complete(impub,nim))
  
  #Cox regression analysis for each imputed set
  analyses[[nim]] <- coxph(Surv(y, failed) ~ X1+as.factor(X2)+X3+X4+as.factor(X5), data = compl, method="breslow")
  r <- 1
  for (i in cuali){
    for (j in 1:n){
      MICE_FCE[nim,r,] <- ifelse(compl[j,i]== data[j,i],MICE_FCE[nim,r,],MICE_FCE[nim,r,]+1)
    }
    r <- r+1}
  MICE_nrmse[nim,,] <- hydroGOF::nrmse(compl[,cuanti], data[,cuanti])
}
FCE[s,,12] <- colMeans(MICE_FCE)
MICE_FCE <- array(0,c(5,length(cuali),1))
nrmse[s,,12] <- colMeans(MICE_nrmse)
MICE_nrmse <- array(0,c(5,length(cuanti),1))

object <- list(call=call, call1=impub$call, nmis=impub$nmis, analyses=analyses)
oldClass(object) <- c("mira", "coxph")
coxmiceb <- summary(pool(object))

estim[s,,12] <- t(coxmiceb[,1])
InfCI[s,,12] <- t(coxmiceb[,6])
SupCI[s,,12] <- t(coxmiceb[,7])
print (c(s,"MICE+ ok"))
}

#-------------------------------------------------------------------------------
#Results evaluation

#Mean square error matrix by estimated parameter and imputation method
for (j in 1:6){
  for (i in 1:12){
    ecm[i,j] <- hydroGOF::mse(estim[,j,i], coefcorr[,j])
  }}
ECM <- round(ecm,5)

#Matrix of means of the NRMSE by quantitative variable and imputation method

for (i in 1:12){
  meannrmse[i,] <- colMeans (nrmse[,,i])
}

#Matrix with CI percentage that cover the correct estimate, by parameter and method

for (s in 1:sim){
  for (j in 1:6){
    for (i in 1:12){
      cobertura[i,j] <- ifelse(coefcorr[1,j]>InfCI[s,j,i] & coefcorr[1,j]<SupCI[s,j,i],
                               cobertura[i,j]+1,cobertura[i,j])  
    }
  }
}

meanFCE <- colMeans(FCE)/(prop*n)*100
meantimes <- colMeans(times)

for (i in 1:12){
  a <- colMeans (estim[,,i])
  print(a)
}