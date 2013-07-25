setwd('where ever')

#install.packages('QSARdata')
library('QSARdata')

#install.packages('monomvn')
library('monomvn')

#install.packages('bootstrap')
library('bootstrap')

#install.packages('KernSmooth')
library('KernSmooth') 
#install.packages('FNN')
library(FNN)


source('Bayesian.Predictions.R')
source('APU_CADD.R')

load('AquaticToxdata.Rdata')
str(Xtrain)  ## study the X-matrix

blasso.out <- blasso(Xtrain,Ytrain, thin = 10,T = 500, RJ = TRUE, verb = 0)
str(blasso.out) ## look at the output from the blasso
## study values on coefficients
head(blasso.out$beta) ## some are zero some dont
image(blasso.out$beta) ## too large to see anything
image(blasso.out$beta[,1:100]) ## some are always zero, some take values for a while
image(blasso.out$beta[,sample.int(dim(Xtrain)[2], size = 100)]) ## just another way of selection which descriptors to look at

## study the simulation
plot(blasso.out$llik,xlab = 'iteration',ylab = 'log likelihood') ##the likelihood for the parameters increase and then saturates
burnin <- 100  # these iterations are thrown away
## note 500 iterations is very small, increasing the number of simulations takes time. Here I do another thing to speed up the simulation

## find which descriptors that from the previous simulation have non-zero values at a certain degree
s <- summary(blasso.out,burnin = 100)
str(s)
plot(s$bn0,xlab = 'descriptor index',ylab = 'proportion of times being non-zero')
select <- s$bn0 > 0.5  
Xtrain.sel <- Xtrain[,select]
Xext.sel <- Xext[,select]

## running this takes time so load the result file instead

## run the Bayesian lasso again with 5000 iterations
#blasso.out <- blasso(Xtrain.sel,Ytrain, thin = 10,T = 5000, RJ = TRUE, verb = 0)

## derive the predictions for the training and external data set using a function in the file Bayesian.R
#qB <- Bayesian(Xtrain.sel,Ytrain,Xext.sel,Yext,blasso.out,burnin=1000,size = 2000)

#save(blasso.out,Xtrain.sel,Xext.sel,Ytrain,Yext,qB,file = 'blasso.AquaticTox.5000.Rdata')

load('blasso.AquaticTox.5000.Rdata')

plot(blasso.out$llik,xlab = 'iteration',ylab = 'log likelihood')
s <- summary(blasso.out,burnin = 1000)
plot(s$bn0,xlab = 'descriptor index',ylab = 'proportion of times being non-zero') ## this time descriptors have more or less influence in the model

## explore the predictive distribtions 
str(qB)
id <- 4 ## select a compound
plot(qB$MCMC.PRED[id,]) ## is seems to be stable over time
hist(qB$MCMC.PRED[id,],main= paste('predictive distribution for compound ',eval(id))) 
plot(range(qB$MCMC.PRED[id,]),c(0,1),col='white',xlab = 'predicted value',ylab = 'cdf',main= paste('predictive distribution for compound ',eval(id)))
Fn <- ecdf(qB$MCMC.PRED[id,]) ## this derives the empirical distributin function
lines(sort(qB$MCMC.PRED[id,]),Fn(sort(qB$MCMC.PRED[id,])))

## now we produce a prediction interval with chosen confience level
conf.lev = 0.95
conf.bounds <- quantile(qB$MCMC.PRED[id,],c((1-conf.lev)/2,1-(1-conf.lev)/2))

## add the prediction inteval in the plot
xval <- c(seq(min(qB$MCMC.PRED[id,]),conf.bounds[1],length.out =50))
polygon(c(min(xval),xval,max(xval)),c(0,Fn(xval),0),col='gray')
xval <- c(seq(max(qB$MCMC.PRED[id,]),conf.bounds[2],length.out =50))
polygon(c(max(xval),xval,min(xval)),c(0,Fn(xval),0),col='gray')
mtext(paste(eval(100*conf.lev),'% pred.int=','[',eval(round(conf.bounds,2)[1]),',',eval(round(conf.bounds,2)[2]),']'))

## add the observed experimental value in the plot
lines(rep(Yext[id],2),c(0,1),col='red')


## derive empirical coverage plots

## a function to caluculate hit rates for a sample with predictive distributions MCMC and observed values Y
emp.cov <- function(MCMC,Y){
d2 <- length(Y)
conf <- seq(0.99,0.01,by = -0.01)
hit <- conf
for(tt in 1:length(conf)){
ss <- 0
for(k in 1:d2){
INT <- quantile(MCMC[k,],probs = c((1-conf[tt])/2,1-(1-conf[tt])/2))
ss <- ss + (Y[k] > INT[1])*(Y[k] < INT[2])
}
hit[tt] <- ss/d2
}
as.data.frame(list(conf=conf,hit=hit))
}

EC.ext <- emp.cov(qB$MCMC.PRED,Yext)  ## the empirical coverage for the external data set
str(EC.ext)
plot(c(0,1),c(0,1),xlim=c(0,1),ylim=c(0,1),xlab='confidence',ylab='hit rate',
main='empirical coverage',col='white')
lines(EC.ext[,'conf'],EC.ext[,'hit'],col=2)
abline(0,1)


## compare with the training data set
EC.train <- emp.cov(qB$MCMC.PRED.TRAIN,Ytrain) 
lines(EC.train[,'conf'],EC.train[,'hit'],col=3,lty=2)  ##predictive distributions are too wide

## derive log likelihood scores
## let us assume the predictive distributions to be Gaussian - if not we can use empirical distributions for this
dnorm(Yext[id],mean(qB$MCMC.PRED[id,]),sqrt(var(qB$MCMC.PRED[id,])))

log.lik.score <- function(MCMC,Y){
d2 <- length(Y)
s <- 0
for(i in 1:d2){
s <- s + log(dnorm(Y[i],mean(MCMC[i,]),sqrt(var(MCMC[i,]))))
}
s
}


log.lik.score(qB$MCMC.PRED,Yext)
log.lik.score(qB$MCMC.PRED.TRAIN,Ytrain)


### let us take another method to compare with 
## bootstrap based on pls
## prepare data for pls
X <- rbind(Xtrain,Xext)
X <- as.matrix(unname(X))
y <- unname(c(Ytrain,Yext))
data <- as.data.frame(list(y=y,X=X))
split <- c(rep(1,length(Ytrain)),rep(2,length(Yext)))

## specify bootstrap function
theta <- function(x,data,split){
sel.data <- data[x,]
sel.split <- split[x]==1
y <- sel.data[,'y']
X <- as.matrix(unname(sel.data[,-1]))
pls.out <- plsr(y ~ X,data = sel.data,subset=sel.split,ncomp=5)
X <- as.matrix(unname(data[,-1]))
predict(pls.out,newdata=X,type='response')[,1,5]
}

## perform the bootsrap - to save time I use 200 samples
#results <- bootstrap(x=which(split==1),200,theta,data,split)
#save(results,data,split,file='bootstrap.AquaTox.Rdata')
load('bootstrap.AquaTox.Rdata')
str(results)
boot.train <- results$thetastar[split==1,] ## extract the predicted values for training data
boot.ext <- results$thetastar[split==2,]

## calcuate empirical coverages
EC.boot.train <- emp.cov(boot.train,Ytrain) 
EC.boot.ext <- emp.cov(boot.ext,Yext) 

## plot and compare with the Bayesian assessment
plot(c(0,1),c(0,1),xlim=c(0,1),ylim=c(0,1),xlab='confidence',ylab='hit rate',
main='empirical coverage',col='white')
lines(EC.ext[,'conf'],EC.ext[,'hit'],col=2)
abline(0,1)
lines(EC.train[,'conf'],EC.train[,'hit'],col=2,lty=2)  ##predictive distributions are too wide
lines(EC.boot.ext[,'conf'],EC.boot.ext[,'hit'],col=3)
lines(EC.boot.train[,'conf'],EC.boot.train[,'hit'],col=3,lty=2)

nam <- c('Blasso.Ext','Blasso.Train','Bootstrap.Ext','Bootstrap.Train')
legend('topleft',c('1:1',nam),lty = c(1,1,2,1,2),col=c(1,2,2,3,3),bty='n')

## plot the log likelihood scores for comparison
LL <- c(log.lik.score(qB$MCMC.PRED,Yext),log.lik.score(boot.ext,Yext))
nam <- c('Blasso.Ext','Bootstrap.Ext')
bp <- barplot(LL,xlab = 'log likelihood score',horiz = TRUE,names.arg = FALSE,
ylab = 'Assessment of predictive error')
text(rep(0,length(bp)),bp,labels = nam,pos = 2)

## are the pls and blasso good models?

pls.out <- plsr(y ~ X,data=data,subset=(split==1),ncomp=10)
pls.out.org <- pls.out
#plot(array(pls.out$validation$PRESS))
X <- rbind(Xtrain,Xext,Xext)
X <- as.matrix(unname(X))
pred.pls <- predict(pls.out,newdata=X,type='response')[,1,4]
scores.pls <- predict(pls.out,newdata=X,type='scores')
split2 <- c(rep(1,length(Ytrain)),rep(2,length(Yext)),rep(3,length(Yext))) 

par(mfrow = c(1,2))
plot(Ytrain,pred.pls[split2==1],ylab = 'predicted',xlab = 'observed',main = 'training data')
points(Ytrain,qB$y.PRED.train,col = 2)
abline(0,1)
mtext(paste('R2_pls =', eval(round(1 - var(Ytrain-pred.pls[split2==1])/var(Ytrain),2)),
'R2_Blasso =', eval(round(1 - var(Ytrain-qB$y.PRED.train)/var(Ytrain),2))))

plot(Yext,pred.pls[split2==2],ylab = 'predicted',xlab = 'observed',main = 'test data')
points(Yext,qB$y.PRED.ext,col=2)
abline(0,1)
mtext(paste('R2_pls =', eval(round(1 - var(Yext-pred.pls[split2==2])/var(Yext),2)),
 'R2_Blasso =',eval(round(1 - var(Yext-qB$y.PRED.ext)/var(Yext),2))))

####
## A third way - assess predictive distriubtion from judgment of predictive reliability

## I cheat and simplify descriptor space by using the components from the pls
D.space <- rbind(scores.pls)
## prepare for input to the coming function
pred <- pred.pls
actual <- c(Ytrain,Yext,Yext)
split2 <- c(rep(1,length(Ytrain)),rep(2,length(Yext)),rep(3,length(Yext)))  
STD.train <- sqrt(apply(boot.train,1,'var'))
STD.ext <- sqrt(apply(boot.ext,1,'var'))
STD <- c(STD.train,STD.ext,STD.ext)

## run Local Assessment of Predictive errors assuming the predictive distribution to be Gaussian
q <- APU_CADD(actual,pred,split2,D.space,STD=STD,ASS = 'EXT',
CP = c('euclidean','leverage','ADdens','std'),DIST = 'Gaussian',plotta = TRUE,
pred.new = pred[1:2],D.space.new = D.space[1:2,],STD.new = STD[1:2])


## run Local Assessment of Predictive errors assuming the predictive distribution to not change in shape of the AD, but no specific shape
q.NP <- APU_CADD(actual,pred,split2,D.space,STD=STD,ASS = 'EXT',
CP = c('euclidean','leverage','ADdens','std'),DIST = 'nonparam',plotta = TRUE,
pred.new = pred[1:2],D.space.new = D.space[1:2,],STD.new = STD[1:2])

## identify the best assessment 
ind <- which(q$LOGLIK[1,]==max(q$LOGLIK[1,]))
names(ind)  ## the name of the best approach
LL <- c(log.lik.score(qB$MCMC.PRED,Yext),q$LOGLIK[1,ind],q$LOGLIK[1,1])

## plot comparision between Blasso and null model having equal errors from PRESS
nam <- c('Blasso',colnames(q$LOGLIK)[ind],'equal')
par(mfrow = c(1,1))
bp <- barplot(LL,xlab = 'log likelihood score',horiz = TRUE,names.arg = FALSE,
ylab = 'Assessment of predictive error',main='Relative comparision')
text(rep(0,length(bp)),bp,labels = nam,pos = 2)




