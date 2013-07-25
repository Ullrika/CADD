data(AquaticTox)
### split into training and external data set
sel.ind <- sample.int(n=length(AquaticTox_Outcome$Activity))
size <- round(length(AquaticTox_Outcome$Activity)*0.60)
sel <- rep(1,length(sel.ind))
sel[sel.ind > size] <- 2
Ytrain <- AquaticTox_Outcome$Activity[sel == 1]
Xtrain <- AquaticTox_Dragon[sel == 1,-1]
Yext <- AquaticTox_Outcome$Activity[sel == 2]
Xext <- AquaticTox_Dragon[sel == 2,-1]
dd <- clean.QSAR.data(Xtrain,Ytrain,Xext,Yext)  ### this function removes NaNs, Infs and missing data
Xtrain <- dd$Xtrain
Xext <- as.matrix(dd$Xext)
Ytrain <- dd$Ytrain
Yext <- dd$Yext
SCALE <- dd$SCALE
YSCALE <- dd$YSCALE
save(Xtrain,Xext,Ytrain,Xext,SCALE,YSCALE,file = 'AquaticToxdata.Rdata')
