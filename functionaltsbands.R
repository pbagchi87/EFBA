library(fields)
library(viridis)
library(MASS)
library(cmvnorm)
library(lattice)

#function to simulate functional white noise data
fws.sim <- function(nb=15,gsz=1000,Ts){
  #The functional white noise data data X_t(u) = c_1f_1(u) + ... + c_df_d(u)
  #f_j's are j-th fourier basis
  #c_j ~ N(0,e^((j-1)/20)), independent, 
  #the variance should increase with increasing j, 
  #to ensure more randomness of the process
  
  # create fourier basis
  #nb <- 15; #number of basis functions, needs to be odd
  #gsz <- 1000; #grid size
  x <- seq(0,1,length.out=gsz)
  fb <- matrix(NA,nrow=gsz,ncol=nb)
  fb[,1] <- 1;
  fb[,seq(2,nb-1,by=2)] <- sqrt(2)*sapply(1:floor(nb/2),FUN=function(n) sin(2*pi*n*x));
  fb[,seq(3,nb,by=2)] <- sqrt(2)*sapply(1:floor(nb/2),FUN=function(n) cos(2*pi*n*x));
  
  # draw coefficients for fourier basis
  covmat <- diag(exp(((1:nb)-1)/20)); # covariance matrix for coefficients
  fcf <- sqrt(covmat)%*%matrix(rnorm(nb*Ts),ncol=Ts);
  
  # draw one time point realization for process
  fwn <- t(fcf)%*%t(fb);
  #plot(x,fwn,type='l')
  return(t(fwn));
}

#function to compute multitaper power spectrum estimator for stationary series
fmts <- function(X,K){
  
  #initialize parameters
  Ts <- nrow(X); #total length of time series
  R <- ncol(X); #number of elements in functional time series
  freq <- seq(from=0,by=1/Ts,length.out=floor(Ts/2)+1); #frequencies
  Fs <- floor(Ts/2)+1; #number of frequencies
  tapers <- outer(1:Ts,1:K,FUN=function(n,k) sqrt(2/(Ts+1))*sin(pi*k*n/(Ts+1))); #sine tapers
  
  #take FFT of tapered series
  ft.tmp <- array(apply(X=X,MARGIN=2,FUN=function(x) mvfft(tapers*x)),dim=c(Ts,K,R))[1:Fs,,];
  
  #compute single-taper auto and cross spectra estimates
  ft.tmp2 <- array(NA,dim=c(Fs,K,R^2))
  for (i in 1:K){
    ft.tmp2[,i,] <- apply(X=ft.tmp[,i,],MARGIN=2,FUN=function(x) sqrt(Conj(x)*ft.tmp[,i,]));
  }
  
  #take the mean across tapers to get multitaper estimates and reformat array
  ft.tmp3 <- array(data=t(apply(X=ft.tmp2,MARGIN=c(1,3),FUN=mean)),dim=c(R,R,Fs));
  
  return(list(freq=freq,mtspec=ft.tmp3));
}

#function to partition series into segments of size N
partN <- function(X,N){
  #initialize parameters
  Ts <- nrow(X); #total length of time series
  R <- ncol(X); #number of elements in functional time series
  B <- floor(Ts/N); #number of partition segments
  if(Ts%%N!=0){
    tmp <- floor(Ts%%N/2+0.5);
    if(Ts%%N%%2!=0){			
      Xarr <- array(data=X[tmp:(Ts-tmp),],dim=c(N,B,R));
      warning(paste("Length of X is not a multiple of N.  First ",tmp-1," and last ",tmp," observations have been discarded.",sep=""));
    }else if(Ts%%N%%2==0){
      Xarr <- array(data=X[(tmp+1):(Ts-tmp),],dim=c(N,B,R));
      warning(paste("Length of X is not a multiple of N.  First ",tmp," and last ",tmp," observations have been discarded.",sep=""));
    }
    Ts <- N*B; #update T
  } else{
    Xarr <- array(data=X,dim=c(N,B,R));
  }
  return(Xarr);
}

#function to linearly detrend a data series
detrend <- function(vec,std){
  #remove linear trend for each interval
  xmat <- cbind(matrix(1,length(vec),1),seq(1,length(vec),length=length(vec)));
  linfit <- solve(t(xmat)%*%xmat)%*%t(xmat)%*%vec;	
  vec <- vec-xmat%*%linfit;
  
  #standardized variance
  if (std){
    vec <- vec/sd(vec);
  }
  
  #return
  return(vec);
}

#function to compute local multitaper estimates
lfmts <- function(Xarr.dt,K){
  #initialize parameters
  N <- dim(Xarr.dt)[1];
  B <- dim(Xarr.dt)[2];
  R <- dim(Xarr.dt)[3];
  Fs <- floor(N/2)+1; #number of frequencies
  freq <- seq(from=0,by=1/N,length.out=Fs); #vector frequencies
  
  #compute local multitaper estimates and organize into array BxFsxRxR
  Xarr.fmts <- aperm(array(data=apply(X=Xarr.dt,MARGIN=2,FUN=function(x) fmts(x,K)$mtspec),
                           dim=c(R,R,Fs,B)),
                     c(4,3,1,2));
  return(list(freq=freq,lfmts=Xarr.fmts));
}

#function to find covariance matrix of demeaned spectral estimates
Q.test <- function(Qts,Qint,X,B,R,K,omega0,delta,ndraw=50000){	
  
  #covariance of G under null
  covg <- array(NA,dim(X)[c(1,1,3,4,3,4)]);
  
  for (b in 1:B){
    for (c in 1:B){
      for (i in 1:R){
        for (j in 1:R){
          for (k in 1:R){
            for (l in 1:R){
              if (b==c){
                covg[b,b,i,j,k,l] <- (1-(2/B))*X[b,omega0+delta,i,k]*X[b,omega0+delta,j,l]+
                  (1/B^2)*sum(X[,omega0+delta,i,k]*X[,omega0+delta,j,l])+
                  (1/(delta^2))*(1-(2/B))*sum(X[b,omega0:(omega0+delta-1),i,k]*X[b,omega0:(omega0+delta-1),j,l])+
                  (1/(B^2*delta^2))*sum(X[,omega0:(omega0+delta-1),i,k]*X[,omega0:(omega0+delta-1),j,l]);
              } else{
                covg[b,c,i,j,k,l] <- (-1/B)*X[b,omega0+delta,i,k]*X[b,omega0+delta,j,l]+
                  (-1/B)*X[c,omega0+delta,i,k]*X[c,omega0+delta,j,l]+
                  (1/B^2)*sum(X[,omega0+delta,i,k]*X[,omega0+delta,j,l])+
                  (-1/(delta^2*B))*sum(X[b,omega0:(omega0+delta-1),i,k]*X[b,omega0:(omega0+delta-1),j,l])+
                  (-1/(delta^2*B))*sum(X[c,omega0:(omega0+delta-1),i,k]*X[c,omega0:(omega0+delta-1),j,l])+
                  (1/(B^2*delta^2))*sum(X[,omega0:(omega0+delta-1),i,k]*X[,omega0:(omega0+delta-1),j,l]);
              }
            }
          }
        }
      }
    }
  }
  
  #rearrange into matrix
  covg.mat <- matrix(NA,nrow=B*R^2,ncol=B*R^2);
  for (b in 1:B){
    for (c in 1:B){
      covg.mat[(((b-1)*R^2)+1):(b*R^2),(((c-1)*R^2)+1):(c*R^2)] <- matrix(covg[b,c,,,,],nrow=R^2,byrow=FALSE)
    }
  }
  
  #draw from GP
  gpdraws <- rcmvnorm(ndraw,rep(0,B*R^2),covg.mat)
  
  #compute Q(tau,sigma) and Qint for each draw
  Qts.nulldist <- array(apply(gpdraws,1,function(x) (1/K)*rowSums(matrix(Mod(x)^2,ncol=B,byrow=FALSE))),dim=c(R,R,ndraw))
  Qint.nulldist <- apply(gpdraws,1,function(x) sum(Mod(x)^2*(1/R^2)*(1/K)))
  Qts.pval <- matrix(rowSums(apply(Qts.nulldist,3,function(x) x>Qts))/ndraw,nrow=R); #DOUBLE CHECK ARRAYS CONFORMED APPROPRIATELY
  Qint.pval <- sum(Qint.nulldist>Qint)/ndraw;
  
  #plot density under null
  x11();plot(density(Qint.nulldist));abline(v=Qint,col='red');print(Qint)
  #x11();plot(density(Qts.nulldist[1,]))
  
  return(list(Qts.pval=Qts.pval,Qint.pval=Qint.pval))
}

#simulate functional white noise
nb <- 7;
gsz <- 20;
Ts <- 3000;
test <- fws.sim(nb=nb,gsz=gsz,T=Ts);

x11();plot(seq(0,1,length.out=gsz),test[,2],type='l');

x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=t(test), 
                 axes = TRUE, col = inferno(256), 
                 main = 'Functional White Noise',xlab='t',ylab='f_t(x)',xaxs="i"); 

#compute multitaper power spectrum estimator
X <- t(test); 
# K <- 7;
# mtspec <- fmts(X,K=K)
# 
# #auto spectra plot
# plot(mtspec$freq,mtspec$mtspec[2,2,],type='l')
# 
# #coherence plot
# plot(mtspec$freq,Mod(mtspec$mtspec[1,4,])^2/(mtspec$mtspec[1,1,]*mtspec$mtspec[4,4,]),type='l',ylim=c(0,1))


#partition series into segments of size N
B=30;
N=Ts/B;
K=max(1,floor(.05*(N+1)-1));

Xarr <- partN(X=X,N=N); #NxBxR

#linearly detrend data in each partition segment
Xarr.dt <- apply(X=Xarr,MARGIN=c(2,3),FUN=function(x) detrend(x,std=FALSE));

#local multitaper estimate
Xarr.fmts <-lfmts(Xarr.dt,K=K);

x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[1,,2,2],type='l')
x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[4,,2,2],type='l')
x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[10,,2,2],type='l')

x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=X,
                 axes = TRUE, col = inferno(256),
                 main = 'Functional White Noise',xlab='t',ylab='f_t(x)',xaxs="i");

Fs <- dim(Xarr.fmts$lfmts)[2];
x11();image.plot(x=(1:B)*(Ts/B),y=Xarr.fmts$freq,z=Re(Xarr.fmts$lfmts[,,1,1]), 
                 axes = TRUE, col = inferno(256), 
                 main = 'Multitaper Autospectrum',xlab='Time',ylab='Hz',xaxs="i"); 

#compute demeaned time-varying spectrum
mu.lfmts <- apply(X=Xarr.fmts$lfmts,MARGIN=c(2,3,4),FUN=mean);
Xarr.g <- sweep(Xarr.fmts$lfmts,c(2,3,4),mu.lfmts);
x11();image.plot(x=(1:B)*(Ts/B),y=Xarr.fmts$freq,z=Re(Xarr.g[,,1,1]), 
                 axes = TRUE, col = inferno(256), 
                 main = 'Demeaned Multitaper Autospectrum',xlab='Time',ylab='Hz',xaxs="i"); 

#compute Qts and Qint
startf <- 2;
endf <- length(Xarr.fmts$freq)-1;
N.srch <- endf-startf; #number of frequencies searching over to find b	
R <- dim(Xarr.g)[3];
B <- dim(Xarr.g)[1];
Qb <- array(data=NA,dim=c(B,N.srch,R,R));
for (l in 1:B){
  for (i in 1:R){
    for (j in 1:R){
      for (k in (startf+1):endf){
        Qb[l,k-startf,i,j] <- Xarr.g[l,k,i,j] - rowMeans(Xarr.g[l,startf:(k-1),i,j,drop=FALSE]) 
      }
    }
  }
}

Qts <- array(apply(apply(X=Qb,MARGIN=c(1,2),FUN=function(x) Mod(x)^2),c(1,3),sum),dim=dim(Qb)[c(4,3,2)])
Qint <- apply(X=Qb,MARGIN=2,FUN=function(x) sum(Mod(x)^2*(1/R^2)))

x11();plot(x=Xarr.fmts$freq[(startf+1):endf],y=Qint,type='l')



ndraw <- 5000
pval.list <- list();
for (i in 1:(endf-startf)){

  i <- which(Xarr.fmts$freq==0.25)-startf
  print(i);
  pval.list[[i]] <- Q.test(Qts=Qts[,,i],Qint=Qint[i],
                           X=Xarr.fmts$lfmts,
                           B=B,R=R,K=K,
                           omega0=startf,delta=i,
                           ndraw=ndraw);
}

plot(Xarr.fmts$freq[(startf+1):endf],
     unlist(lapply(pval.list, `[[`, 2)),type='l')

Qtsarr <- array(unlist(lapply(pval.list, `[[`, 1)),dim=c(R,R,endf-startf))






































#code to simulated banded structure nonstationary functional TS
nb <- 5;
gsz <- 10;
Ts <- 5000;

#low frequencies
test <- fws.sim(nb=nb,gsz=gsz,T=Ts);
X <- t(test); #100 elements in functional time series observed over 5000 time points (500 x 100)
Ts <- nrow(X);
Fs <- floor(Ts/2)+1;
f <- seq(from=0,by=1/Ts,length.out=floor(Ts/2)+1);
f <- c(f,rev(f[c(-1,-which(f==0.5))])); 
dft <- mvfft(X)/Ts;
dft[which(f>0.15),] <- 0; 
fwn1.lf <- Re(mvfft(dft,inverse=TRUE));
fwn1.lf <- fwn1.lf/sd(fwn1.lf);

# mtspec <- fmts(fwn1.lf,K=10)
# x11();plot(seq(0,1,length.out=gsz),fwn1.lf[2,],type='l');
# x11();plot(1:Ts,fwn1.lf[,2],type='l');
# x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=fwn1.lf,
#                  axes = TRUE, col = inferno(256),
#                  main = 'Functional White Noise',xlab='t',ylab='f_t(x)',xaxs="i");
# x11();plot(mtspec$freq,mtspec$mtspec[2,2,],type='l')

#middle frequencies
test <- fws.sim(nb=nb,gsz=gsz,T=Ts);
X <- t(test); #100 elements in functional time series observed over 5000 time points (500 x 100)
Ts <- nrow(X);
Fs <- floor(Ts/2)+1;
f <- seq(from=0,by=1/Ts,length.out=floor(Ts/2)+1);
f <- c(f,rev(f[c(-1,-which(f==0.5))])); 
dft <- mvfft(X)/Ts;
dft[which(f<=0.15 | f>0.35),] <- 0;
fwn1.mf <- Re(mvfft(dft,inverse=TRUE));
fwn1.mf <- fwn1.mf/sd(fwn1.mf);

# mtspec <- fmts(fwn1.mf,K=10)
# x11();plot(seq(0,1,length.out=gsz),fwn1.mf[2,],type='l');
# x11();plot(1:Ts,fwn1.mf[,2],type='l');
# x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=fwn1.mf, 
#                  axes = TRUE, col = inferno(256), 
#                  main = 'Functional White Noise',xlab='t',ylab='f_t(x)',xaxs="i"); 
# x11();plot(mtspec$freq,mtspec$mtspec[2,2,],type='l')

#high frequencies
test <- fws.sim(nb=nb,gsz=gsz,T=Ts);
X <- t(test); #100 elements in functional time series observed over 5000 time points (500 x 100)
Ts <- nrow(X);
Fs <- floor(Ts/2)+1;
f <- seq(from=0,by=1/Ts,length.out=floor(Ts/2)+1);
f <- c(f,rev(f[c(-1,-which(f==0.5))])); 
dft <- mvfft(X)/Ts;
dft[which(f<=0.35),] <- 0; 
fwn1.hf <- Re(mvfft(dft,inverse=TRUE));
fwn1.hf <- fwn1.hf/sd(fwn1.hf);

# mtspec <- fmts(fwn1.hf,K=10)
# x11();plot(seq(0,1,length.out=gsz),fwn1.hf[2,],type='l');
# x11();plot(1:Ts,fwn1.hf[,2],type='l');
# x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=fwn1.hf, 
#                  axes = TRUE, col = inferno(256), 
#                  main = 'Functional White Noise',xlab='t',ylab='f_t(x)',xaxs="i"); 
# x11();plot(mtspec$freq,mtspec$mtspec[2,2,],type='l')


#nonstationary 3 segments linear
#combine
coef1 <- seq(from=sqrt(100),to=1,length.out=Ts);
coef2 <- seq(from=1,to=1,length.out=Ts);
coef3 <- seq(from=1,to=sqrt(100),length.out=Ts);
X.3bL <- coef1*fwn1.lf*sqrt(.3) + coef2*fwn1.mf*sqrt(.4) + coef3*fwn1.hf*sqrt(.3);

x11();plot(seq(0,1,length.out=gsz),X.3bL[1,],type='l');
x11();plot(1:Ts,X.3bL[,2],type='l');
x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=X.3bL,
                 axes = TRUE, col = inferno(256),
                 main = 'Nonstationary Functional Banded',xlab='t',ylab='f_t(x)',xaxs="i");

#partition series into segments of size N
N=500
Xarr <- partN(X=X.3bL,N=N); #NxBxR

#linearly detrend data in each partition segment
Xarr.dt <- apply(X=Xarr,MARGIN=c(2,3),FUN=function(x) detrend(x,std=FALSE));

#local multitaper estimate
K=max(1,floor(.05*(N+1)-1));
Xarr.fmts <-lfmts(Xarr.dt,K=K);
  
x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[1,,2,2],type='l')
x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[4,,2,2],type='l')
x11();plot(Xarr.fmts$freq,Xarr.fmts$lfmts[10,,2,2],type='l')

x11();image.plot(x=1:Ts,y=seq(0,1,length.out=gsz),z=X.3bL,
                 axes = TRUE, col = inferno(256),
                 main = 'Nonstationary Functional Banded',xlab='t',ylab='f_t(x)',xaxs="i");

B <- dim(Xarr.fmts$lfmts)[1];
Fs <- dim(Xarr.fmts$lfmts)[2];
x11();image.plot(x=(1:B)*(Ts/B),y=Xarr.fmts$freq,z=Re(Xarr.fmts$lfmts[,,1,1]), 
                 axes = TRUE, col = inferno(256), 
                 main = 'Multitaper Autospectrum',xlab='Time',ylab='Hz',xaxs="i"); 
abline(h=c(0.15,0.35),col='green');

#compute demeaned time-varying spectrum
mu.lfmts <- apply(X=Xarr.fmts$lfmts,MARGIN=c(2,3,4),FUN=mean);
Xarr.g <- sweep(Xarr.fmts$lfmts,c(2,3,4),mu.lfmts);
x11();image.plot(x=(1:B)*(Ts/B),y=Xarr.fmts$freq,z=Re(Xarr.g[,,1,1]), 
                 axes = TRUE, col = inferno(256), 
                 main = 'Demeaned Multitaper Autospectrum',xlab='Time',ylab='Hz',xaxs="i"); 
abline(h=c(0.15,0.35),col='green');

##############################################
## PICK UP FROM HERE
###########################################

#compute Qts and Qint
startf <- 2;
endf <- length(Xarr.fmts$freq)-1;
N.srch <- endf-startf; #number of frequencies searching over to find b	
R <- dim(Xarr.g)[3];
B <- dim(Xarr.g)[1];
Qb <- array(data=NA,dim=c(B,N.srch,R,R));
for (l in 1:B){
  for (i in 1:R){
    for (j in 1:R){
      for (k in (startf+1):endf){
        Qb[l,k-startf,i,j] <- Xarr.g[l,k,i,j] - rowMeans(Xarr.g[l,startf:(k-1),i,j,drop=FALSE]) 
      }
    }
  }
}

Qts <- array(apply(apply(X=Qb,MARGIN=c(1,2),FUN=function(x) Mod(x)^2),c(1,3),sum),dim=dim(Qb)[c(4,3,2)])
Qint <- apply(X=Qb,MARGIN=2,FUN=function(x) sum(Mod(x)^2*(1/R^2)))

x11();plot(x=Xarr.fmts$freq[(startf+1):endf],y=Qint,type='l')

ndraw <- 5000
pval.list <- list();
for (i in 1:(endf-startf)){
  print(i);
  pval.list[[i]] <- Q.test(Qts=Qts[,,i],Qint=Qint[i],
                           X=Xarr.fmts$lfmts,
                           B=B,R=R,K=K,
                           omega0=startf,delta=i,
                           ndraw=ndraw);
}

plot(Xarr.fmts$freq[(startf+1):endf],
        unlist(lapply(pval.list, `[[`, 2)),type='l')

Qtsarr <- array(unlist(lapply(pval.list, `[[`, 1)),dim=c(R,R,endf-startf))




#add Q function code
###function to find frequency partition point using Hochberg step up rule
eba.b <- function(X.dm,f,startf,endf,covg,alpha) {
  
  #initialize data container for results
  out <- vector(mode='numeric',length=6);
  names(out) <- c('bhat.idx','bhat','bhat.Q','bhat.pval','bhat.pval.th','bhat.sig');
  
  #Obtain vector of Q values for different endpoints
  N.srch.b <- endf-startf; #number of frequencies searching over to find b	
  X.dm.cumrowsum <- apply(X=X.dm[startf:(endf-1),,drop=FALSE],MARGIN=2,FUN=cumsum);
  Qb <- apply(X=(X.dm[(startf+1):endf,,drop=FALSE]-(1:N.srch.b)^-1*X.dm.cumrowsum)^2,MARGIN=1,FUN=sum); 	
  #x11();plot(f[(startf+1):endf],Qb,type="l",xlab="Hz",ylab="Q",mgp=c(3,1,0));
  
  #Obtain corresponding vector of p-values for different endpoints
  bidx <- (1:length(Qb))+1;
  #clusterExport(c1, 'covg', envir=environment())
  sigma2 <- pmax(sapply(X=bidx,FUN=function(idx) covg[idx,idx,] 
                        + (-2/(idx-1))*colSums(matrix(covg[idx,1:(idx-1),],nrow=idx-1))
                        + (1/(idx-1)^2)*colSums(matrix(covg[1:(idx-1),1:(idx-1),],nrow=(idx-1)^2))),1e-16);
  #clusterExport(c1, c('sw','sigma2','Qb'), envir=environment())
  pval.b <- sapply(X=1:length(Qb),FUN=function(x) 1-sw(sigma2[,x],Qb[x]));
  #x11();plot(f[(startf+1):endf],pval.b,type="l",xlab="Hz",ylab="p-value",mgp=c(3,1,0));
  
  #Hochberg step up procedure for determining significance
  pval.ord <- sort(pval.b);
  pval.th <- alpha/((endf-startf)-(1:(endf-startf))+1);
  pval.rej <- as.numeric(pval.ord<pval.th);
  
  
  stp <- 0;
  i<-length(pval.rej);
  while (stp == 0 & i>=1){
    if (pval.rej[i]==1){
      pval.rej[1:i] <- 1;
      stp <- 1;
    } else{
      i <- i-1;
    }
  }
  
  pval.idx <- which(pval.b %in% pval.ord[as.logical(pval.rej)])[1]; #smallest significant frequency 
  
  if(is.na(pval.idx)){#if no rejections, set pval index to frequency with smallest p-value
    pval.idx <- which(pval.b==min(pval.b));
    b.sig.ind <- 0;
  } else {
    b.sig.ind <- 1;
  }
  
  bhat <- startf+pval.idx;
  
  #output results for bhat
  out[c('bhat.idx','bhat','bhat.Q','bhat.pval','bhat.pval.th','bhat.sig')] <- c(bhat,f[bhat],Qb[pval.idx],
                                                                                pval.b[pval.idx],pval.th[which(pval.ord==pval.b[pval.idx])],b.sig.ind);
  return(out);
  
}




###############################################################################
###############################################################################
###############################################################################
###############################################################################
# ADAPT LATER
###############################################################################
###############################################################################
###############################################################################
###############################################################################


#nonstationary 3 segments sinusoidal
#low frequency series
wn1 <- rnorm(T,0,1);
dft <- fft(wn1)/T; #discrete Fourier transform (DFT)
dft[which(f>0.15)] <- 0; #0 out frequencies above f=0.15
wn1.new <- Re(fft(dft,inverse=TRUE));
wn1.new <- wn1.new/sd(wn1.new);

#middle band frequency series
wn2 <- rnorm(T,0,1);
dft <- fft(wn2)/T; #discrete Fourier transform (DFT)
dft[which(f<=0.15 | f>0.35)] <- 0; #0 out frequencies not between 0.15 < f <= 0.35
wn2.new <- Re(fft(dft,inverse=TRUE));
wn2.new <- wn2.new/sd(wn2.new);

#high frequency series
wn3 <- rnorm(T,0,1);
dft <- fft(wn3)/T; #discrete Fourier transform (DFT)
dft[which(f<=0.35)] <- 0; #0 out frequencies <= 0.35
wn3.new <- Re(fft(dft,inverse=TRUE));
wn3.new <- wn3.new/sd(wn3.new);  

#combine
coef1 <- sqrt(9)*sin(4*pi*seq(0,1,length=T));
coef2 <- sqrt(9)*cos(4*pi*seq(0,1,length=T));
coef3 <- sqrt(9)*cos(8*pi*seq(0,1,length=T)); 
X.3bS <- coef1*wn1.new*sqrt(.3) + coef2*wn2.new*sqrt(.4) + coef3*wn3.new*sqrt(.3)+rnorm(T,0,1);
