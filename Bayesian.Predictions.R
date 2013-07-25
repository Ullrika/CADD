
##############################

Bayesian <- function(Xtrain,Ytrain,Xext,Yext,blasso.out,burnin,size=10000){

###### Fit a Bayesian model
#blasso.out <- blasso(Xtrain,Ytrain, thin = 10,T = 5000, RJ = TRUE, verb = 0)

## make predictions
N <- 1 #number of samples for each iteration
MCMC.sample.size <- min(2000,blasso.out$T-100)
#burnin <- blasso.out$T-MCMC.sample.size
pred <- array(0,c(MCMC.sample.size,N))

d2 <- dim(Xext)[1]
if(length(d2)==0){
d2 <- length(Xext)
## replace missing values in Xext with zero
Xext[!is.finite(Xext)] <- 0
MCMC.PRED <- array(0,c(d2,N*MCMC.sample.size))
for(k in 1:d2){
for(i in 1:MCMC.sample.size){

pred[i,] <- rnorm(N,blasso.out$mu[burnin+i]+Xext[k]%*%blasso.out$beta[burnin+i],sqrt(blasso.out$s2[burnin+i]))
}
MCMC.PRED[k,] <- array(pred)
}
}else{

MCMC.PRED <- array(0,c(d2,N*MCMC.sample.size))
## replace missing values in Xext with zero
for(k in 1:dim(Xext)[1]){
Xext[k,!is.finite(as.matrix(Xext[k,]))] <- 0
}
for(k in 1:d2){
for(i in 1:MCMC.sample.size){
pred[i,] <- rnorm(N,blasso.out$mu[burnin+i]+blasso.out$beta[burnin+i,]%*%Xext[k,],sqrt(blasso.out$s2[burnin+i]))
}
MCMC.PRED[k,] <- array(pred)
}
}



## calculate log likelihood score 
log.lik.score <- 0
for(k in 1:d2){
log.lik.score <- log.lik.score + log(dnorm(Y[k],mean(MCMC.PRED[k,]),sqrt(var(MCMC.PRED[k,]))))
}


## calculate empirical coverages
conf <- seq(0.99,0.01,by = -0.01)
hit <- conf
for(tt in 1:length(conf)){
ss <- 0
for(k in 1:d2){
INT <- quantile(MCMC.PRED[k,],probs = c((1-conf[tt])/2,1-(1-conf[tt])/2))
ss <- ss + (Yext[k] > INT[1])*(Yext[k] < INT[2])
}
hit[tt] <- ss/d2
}


#############
## do the same thing for the training data
d1 <- dim(Xtrain)[1]
if(length(d1)==0){
d1 <- length(Xtrain)
MCMC.PRED.TRAIN <- array(0,c(d1,N*MCMC.sample.size))
for(k in 1:d1){
for(i in 1:MCMC.sample.size){
pred[i,] <- rnorm(N,blasso.out$mu[burnin+i]+Xtrain[k]%*%blasso.out$beta[burnin+i],sqrt(blasso.out$s2[burnin+i]))
}
MCMC.PRED.TRAIN[k,] <- array(pred)
}
}else{

MCMC.PRED.TRAIN <- array(0,c(d1,N*MCMC.sample.size))
for(k in 1:d1){
for(i in 1:MCMC.sample.size){
pred[i,] <- rnorm(N,blasso.out$mu[burnin+i]+blasso.out$beta[burnin+i,]%*%as.matrix(Xtrain[k,]),sqrt(blasso.out$s2[burnin+i]))
}
MCMC.PRED.TRAIN[k,] <- array(pred)
}
}

## calculate log likelihood scores assuming Gaussian distribution (I have an alteratnive with empirical density function
log.lik.score.TRAIN <- 0
for(k in 1:d1){
log.lik.score.TRAIN <- log.lik.score.TRAIN + log(dnorm(Y[k],mean(MCMC.PRED.TRAIN[k,]),sqrt(var(MCMC.PRED.TRAIN[k,]))))
}


## calculate empirical coverages
conf <- seq(0.99,0.01,by = -0.01)
hit.train <- conf
for(tt in 1:length(conf)){
ss <- 0
for(k in 1:d1){
INT <- quantile(MCMC.PRED.TRAIN[k,],probs = c((1-conf[tt])/2,1-(1-conf[tt])/2))
ss <- ss + (Ytrain[k] > INT[1])*(Ytrain[k] < INT[2])
}
hit.train[tt] <- ss/d1
}


####
#y.PRED.ext <- apply(MCMC.PRED,1,'mean')
#plot(y.PRED.ext,Yext)
#plot(conf,hit.train)

return(list(MCMC.PRED=MCMC.PRED,MCMC.PRED.TRAIN = MCMC.PRED.TRAIN,
hit.rate=cbind(conf,hit.train,hit),LOGLIK=c(log.lik.score.TRAIN,log.lik.score),
y.PRED.ext=apply(MCMC.PRED,1,'mean'),y.PRED.train = apply(MCMC.PRED.TRAIN,1,'mean'))
)

}


########

clean.QSAR.data <- function(Xtrain,Ytrain,Xext,Yext){ 
ok <- complete.cases(Xtrain)
Xtrain <- Xtrain[ok,]
Ytrain <- Ytrain[ok]

## scale data before modelling
temp <- Xtrain
temp.ext <- Xext
S <- scale(temp)
center <- attr(S,'scaled:center')
scl <- array(attr(S,'scaled:scale'))
S.ext <- Xext
for(k in 1:dim(S)[2]){
S.ext[,k] <- (temp.ext[,k] - center[k])/scl[k]
}
ok <- complete.cases(t(S))
Xtrain <- S[,ok]
Xext <- S.ext[,ok]

center <- attr(S,'scaled:center')[ok]
scl <- attr(S,'scaled:scale')[ok]

SY <- scale(Ytrain)
Ytrain <- SY
Yext <- (Yext  - attr(SY,'scaled:center'))/attr(SY,'scaled:scale')

Ycenter <- attr(SY,'scaled:center')
Yscl <- attr(SY,'scaled:scale')

## replace missing values in Xext with mean values of the training data
for(k in 1:dim(Xext)[2]){
Xext[!is.finite(as.matrix(Xext[,k])),k] <- mean(Xtrain[,k])
}

return(list(Xtrain=Xtrain,Xext=Xext,Ytrain=Ytrain,Yext=Yext,
SCALE=cbind(center,scl),YSCALE=c(Ycenter,Yscl)))
}


