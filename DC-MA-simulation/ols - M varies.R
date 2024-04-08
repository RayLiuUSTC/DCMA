rm( list = ls ( all = TRUE))#remove all variables
#options(warn = 0)
setwd("~/pic M varies new")

library(snowfall)

md=3+2+3*2 #five methods:PD,DC_{ew},Full,DC_{opt} and CSL
MM<-c(seq(from=10,to=100,length.out=10)[-10],
      seq(from=100,to=1000,length.out=10)) #the setting for number of sites
vaR<-c(0.5,1,123) #the setting for variance across different sites
ns<-c(50,100,1000) #the setting for local sample size

mseM <-array(0,c(length(MM),md,length(ns),length(vaR))) #MSE(mean):site, method, local sample size and variance
mseS <-array(0,c(length(MM),md,length(ns),length(vaR))) #MSE(sd):site, method, local sample size and variance
estM <-array(0,c(length(MM),md,length(ns),length(vaR))) #est(mean):site, method, local sample size and variance
estS <-array(0,c(length(MM),md,length(ns),length(vaR))) #est(sd):site, method, local sample size and variance
msename <- c('partial-data','equal-weight','full-data.AIC','optimal-weight.solnp','optimal-weight.osqp',paste0("CSL-ave-a3",1:3),paste0("CSL-cen-a3",1:3))
timecostM <-array(0,c(length(MM),3,length(ns),length(vaR)))
wM <- list()
length(wM) <- length(MM)*length(ns)*length(vaR)
dim(wM) <- c(length(MM),length(ns),length(vaR)) #weights(solnp, more sparse): each element of the array store a list of the weights ditributed to M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m
w1M <- list()
length(w1M) <- length(MM)*length(ns)*length(vaR)
dim(w1M) <- c(length(MM),length(ns),length(vaR)) #weights(solnp, more sparse): each element of the array store a list of the weights ditributed to M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m

by=0 #difference parameter, used for 'seq' function, seq(from,by,length.out), to set local sample size and allow it different across different site
p0<-50 #the number of potential covariates
pM<-c(15,30) #the setting for the number of covariates on hand
alphaM<-c(0.5,1,1.5) #the setting for the decay rate of coefficient
G<-100  #the setting for replicates
numcpu <- 25 #the setting for the amount of CPUs requested for the cluster

extractAIC.lm <- function (fit, scale = 0, k = 2, ...) 
{
  n <- length(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- deviance.lm(fit)
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  sigma <- RSS/fit$df.residual
  c(edf, dev + k * edf, sigma)
}

deviance.lm <- function (object, ...) {
  return(sum(weighted.residuals(object)^2, na.rm = TRUE))
}

simfunc <- function(g){
  set.seed(1234+g)
  
  XX<-list()
  YY<-list()
  MU<-list()
  
  
  for (j in 1:M){#DGP for M sites
    nm<-nn[j]#the sample size of jth site
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0)#generate the covariates at the jth site
    Xm<-X0m[,1:p]#take the first p covariates as the ones on hand
    vm<-rnorm(nm,0,sig[j])#generate the errors of the jth site
    mum<-X0m%*%beta0#the true mean vector of the jth site
    Ym<-mum+vm#the response Y_m
    rm(X0m)
    
    for (i in 1:p){#for the jth site, estimate total p nested candidate models
      ols <- RcppEigen::fastLmPure(X=as.matrix(Xm[,1:i]), y=Ym)#without intercept already
      betahat[,i,j][1:i] <- ols$coefficients
      temp.info <- extractAIC.lm(ols)
      GCV[j,i] <- temp.info[2]
      sigmap[j,i] <- temp.info[3]
    }
    bestGCV <- which.min(GCV[j,])
    df[j] <- bestGCV
    sigma[j] <- sigmap[j,bestGCV]
    BetaHat[,j]<-betahat[,bestGCV,j]#at jth site, store the coefficients estimator of the candidate model with minimum AIC into matrix BetaHat
    
    XX[[j]]<-   Xm#store the design matrix on hand at jth site into XX list
    YY[[j]]<-   Ym#store the response vector at jth site into YY list
    MU[[j]]<-   mum#store the mean vector at jth site into MU list
  }
  
  ee <- matrix(data = FALSE,nrow = M,ncol = M)
  
  for (j in 1:M){#validate the prediction model picked at the jth site by the data of other M-1 sites
    Xm<-XX[[j]]
    Ym<-YY[[j]]
    em<-matrix(0,nn[j],M)##em is a matrix to record the residuals of M prediction models at the data of jth site
    em<-matrix(Ym,ncol = M,nrow = nn[j])-Xm%*%BetaHat
    ee <- ee + t(em)%*%em
  }
  #ee<-Reduce("+",emm)#don't use the command reduce here, because once M>700, even when nm=50, the program will run out of memory
  #emm size: 700 lists, with dims 700*700, uses 2.7GB memories
  
  a1 <- ee#step4:\sum_{m=1}^M Q_m
  a2<-sigma*df #compute \psai in step4
  w0 <- matrix(1,nrow=M,ncol=1)/M #equal weight vector
  
  ###########################################################################
  eq=function(w){sum(w)-1} #The equality constraint function
  obj = function(w){ crossprod(w,a1)%*%w+2*crossprod(w,a2)}#The optimaized objective:C(w)
  
  ######################################################################
  result1 <- solnp(matrix(1,nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result2 <- solnp(matrix(seq(from=0.01,by=0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result3 <- solnp(w0,obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result4 <- solnp(matrix(seq(from=0.01*M,by=-0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  bestobj <- which.min(c(tail(result1$values,n=1),tail(result2$values,n=1),tail(result3$values,n=1),tail(result4$values,n=1)))
  
  bestresult <- get(paste0("result",bestobj))
  w = bestresult$pars#note that in result, the argument 'UB' is not considered, probably because it is enough for the restraint that has been already considered
  betahatOMA_solnp=BetaHat%*%w;#the proposed estimator of beta: OPTW
  
  #######osqp solver
  
  A=rbind(rep(1,M),diag(1,M))
  l=c(1,rep(0,M))
  u=c(1,rep(1,M))
  res <- solve_osqp(P=a1*2,q=a2*2,A=A,l=l,u=u,pars = osqpSettings(verbose = FALSE))
  w1 <- res$x
  w1[w1<0] <- 0
  w1 <- w1/sum(w1)
  betahatOMA_osqp=BetaHat%*%w1;#the proposed estimator of beta: OPTW
  
  
  betahatEMA=BetaHat%*%w0;#the equal-weight estimator of beta: EW
  
  
  X<-c() #matrix(0,sum(nn),p)
  Y<-c() #matrix(0,sum(nn),1)
  if (vaR[vari]==123){#for the heteroscadatic situation, vari==123,sig<-sqrt(0.5+(1:M)/M)
    for (j in 1:M){
      X<-rbind(X,sigma[j]^(-0.5)*XX[[j]])
      Y<-rbind(Y,sigma[j]^(-0.5)*YY[[j]])
    }
    for (i in 1:p){#estimate the coefficients of nested candidate models after pooling the data from all sites
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      temp.infof <- extractAIC.lm(olsf)
      GCVf[1,i] <- temp.infof[2]#store the AIC value
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestAICf <- which.min(GCVf[1,])
    betahatfAIC<-betahatfms[,bestAICf]
  }  else  {
    for (j in 1:M){
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }
    for (i in 1:p){#estimate the coefficients of nested candidate models after pooling the data from all sites
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      temp.infof <- extractAIC.lm(olsf)
      GCVf[1,i] <- temp.infof[2]#store the AIC value
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestAICf <- which.min(GCVf[1,])
    betahatfAIC<-betahatfms[,bestAICf]
  }
  
  
  
  #########
  gradient <- function(beta,Y,X){
    (-2*crossprod(X,Y)+2*crossprod(X)%*%beta)/length(Y)
  }
  
  Ltilde = function(beta,gradc1=gradc,gradg1=gradg){crossprod(YY[[Centersite.index]]-XX[[Centersite.index]]%*%beta)/length(YY[[Centersite.index]]) - crossprod(gradc1-gradg1,beta)}
  #########
  
  #time.ols_large_per_site <- vector(length = M)
  #for (j in 1:M){
  #  ptm <- proc.time()
  #  ols <- RcppEigen::fastLmPure(X=as.matrix(XX[[j]]), y=YY[[j]])#without intercept already
  #  betahat[,p,j][1:p] <- ols$coefficients
  #  time.ols_large_per_site[j] <- (proc.time()-ptm)[3]
  #}
  #time.ols_large_per_site.max <- max(time.ols_large_per_site)
  
  Centersite.index <- floor(M/2)#centeral site
  betainitial2 <- betahat[,p,Centersite.index]#if you take betainitial2 into func 'gradient', you will find the column corresponding to centersite.index of the matrix gradM below will be 0!
  
  betainitial <- apply(X=betahat[,p,],MARGIN = 1,FUN = mean)
  
  Iterative <- function(n,betainit,method=rep(1,n),Center=Centersite.index,X=XX,Y=YY){
    betatildeN <- matrix(nrow = length(betainit),ncol = n)
    for (i in 1:n) {
      if(i>1){betainit=betatildeN[,(i-1)]}
      gradM <- mapply(FUN=gradient,Y=Y,X=X,MoreArgs = list(beta=betainit))
      gradg <- apply(X=gradM,MARGIN = 1,FUN=mean)
      if(method[i]==1){
        hessiancinv <- solve(2*crossprod(X[[Center]])/nrow(X[[Center]]))
        betatildeH <- betainit-hessiancinv%*%gradg
        betatildeN[,i] <- betatildeH
      } else {#exact est
        gradc <- gradient(beta = betainit,Y=Y[[Center]],X=X[[Center]])
        betatildeE = solnp(matrix(rep(1/p,p),p,1),Ltilde,control=list(trace=0),gradc=gradc,gradg=gradg)$pars
        betatildeN[,i] <- betatildeE
      }
    }
    return(betatildeN)
  }
  
  betatildeN <- Iterative(n=3,betainit = betainitial)#the simulation like Jordan's
  betatildeN2 <- Iterative(n=3,betainit = betainitial2)#centersite est of beta as the initial input
  
  #betatildeNE <- Iterative(n=3,betainit = betainitial,method = c(1,1,2))
  #betatildeNE2 <- Iterative(n=3,betainit = betainitial,method = c(2,1,1))
  
  
  
  Yst<-c();Xst<-c();must<-c()
  for (j in 1:M){
    #generate the predictive data with the same size as the training data. Note that here the same variable names are used as the training part, thus the training data will be covered.
    nm<-nn[j]
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0)
    Xm<-X0m[,1:p]  
    vm<-rnorm(nm,0,sig[j])
    mum<-X0m%*%beta0
    Xst<-rbind(Xst,Xm)#pool the design matrices from all sites altogether
    must<-rbind(must,mum)
    rm(X0m)
  }
  
  storeloss[1] <-mean((must-Xst%*%BetaHat[,M])^2)#prediction:PD
  storeloss[2] <-mean((must-Xst%*%betahatEMA)^2)#prediction:DC_{ew}
  storeloss[3] <-mean((must-Xst%*%betahatfAIC)^2)#prediction:Full
  storeloss[4] <-mean((must-Xst%*%betahatOMA_solnp)^2)#prediction:DC_{opt}
  storeloss[5] <-mean((must-Xst%*%betahatOMA_osqp)^2)#prediction:DC_{opt}
  for (i in 1:3) {
    storeloss[5+i] <- mean((must-Xst%*%betatildeN[,i])^2)#prediction:CSL(ave init)
  }
  for (i in 1:3) {
    storeloss[8+i] <- mean((must-Xst%*%betatildeN2[,i])^2)#prediction:CSL(centersite init)
  }
  
  #time.pred_CSL3 <- vector(length = 3)
  #for (i in 1:3) {
  #  ptm16 <- proc.time()
  #  storeloss[12+i] <- mean((must-Xst%*%betatildeN3[,i])^2)#prediction:CSL(ave init,method=2,2,2)
  #  time.pred_CSL3[i] <- (proc.time()-ptm16)[3]
  #}
  #
  #time.pred_CSL4 <- vector(length = 3)
  #for (i in 1:3) {
  #  ptm17 <- proc.time()
  #  storeloss[15+i] <- mean((must-Xst%*%betatildeN4[,i])^2)#prediction:CSL(ave init,method=2,1,1)
  #  time.pred_CSL4[i] <- (proc.time()-ptm17)[3]
  #}
  
  esterror[1] <-sum((beta0[1:p]-BetaHat[,M])^2)#prediction:PD
  esterror[2] <-sum((beta0[1:p]-betahatEMA)^2)#prediction:DC_{ew}
  esterror[3] <-sum((beta0[1:p]-betahatfAIC)^2)#prediction:Full
  esterror[4] <-sum((beta0[1:p]-betahatOMA_solnp)^2)#prediction:DC_{opt}
  esterror[5] <-sum((beta0[1:p]-betahatOMA_osqp)^2)#prediction:DC_{opt}
  for (i in 1:3) {
    esterror[5+i] <-sum((beta0[1:p]-betatildeN[,i])^2)#prediction:CSL(ave init)
  }
  for (i in 1:3) {
    esterror[8+i] <-sum((beta0[1:p]-betatildeN2[,i])^2)#prediction:CSL(centersite init)
  }
  
  return(list(storeloss=storeloss,esterror=esterror,w=w,w1=w1))
}

sfInit(parallel = T,cpus = numcpu)

sfLibrary(Rsolnp)
sfLibrary(osqp)
sfExport('md')
sfExport('MM')
sfExport('vaR')
sfExport('ns')
sfExport('mseM')
sfExport('by')
sfExport('p0')
sfExport('G')
sfExport('extractAIC.lm')
sfExport('deviance.lm')



for (alphaind in 1:length(alphaM)) {
  alpha <- alphaM[alphaind]
  for (pind in 1:length(pM)) {
    p <- pM[pind]
    #generate the coefficients
    beta0<-sqrt(2*alpha)*c((1:p0)^(-alpha-0.5))#(runif(p0))
    
    sfExport('p')
    sfExport('alpha')
    sfExport('beta0')
    
    for (vari in 1:length(vaR)) { #start with \sigma^2=0.5
      for (nsi in 1:length(ns)) {
        n1=ns[nsi] #start with n_m=50
        for (mi in 1:length(MM)) {
          M<-MM[mi] #start with M=10
          nn<-seq(n1, by=by,length.out=M) #distribute sample size for every site
          
          betahat<-array(0,c(p,p,M))#c(p,p,M): the coefficients length of betahat, the number of nested candidate models at each site, the number of sites
          BetaHat<-array(0,c(p,M))#c(p,M): the betahat of the nested model picked by AIC, the number of sites
          betahatfms<-matrix(0,nrow=p,ncol=p)#nrow=p,ncol=p: the coefficients length of betahat, the number of nested candidate models
          #em<-matrix(0,n,M)
          storeloss <-vector("numeric",length=md) #prediciton loss for five methods
          esterror <-vector("numeric",length=md) #error of beta for five methods
          #storelossb <-matrix(0,G,md)
          GCV<-matrix(0,M,p) #IC value for local prediction model: the number of sites, the number of covariates on hand
          sigmap<-matrix(0,M,p) #variance: the number of sites, the number of covariates on hand
          GCVf<-matrix(0,1,p)#IC value for full data model: the number of sites, the number of covariates on hand
          #note that 'sig' represents standard deviation which will be used in rnorm(mean=0,sd=sig)
          if (vaR[vari]==0.5){sig<-sqrt(0.5)*rep(1,M)} else #constant variance: \sigma^2=0.5
            if(vaR[vari]==1){sig<-1*rep(1,M)} else #constant variance: \sigma^2=1
            {sig<-sqrt(0.5+(1:M)/M)}#varied variance: \sigma^2=0.5+(1:M)/M
          #note that 'sigma' represents variance which will be used in criteria C(w)
          sigma<-matrix(0,M,1)#the matrix to store the variance ests of M local prediction models
          df<-matrix(0,M,1)#the matrix to store the df ests of M local prediction models
          
          sfExport('vari')
          sfExport('M')
          sfExport('nn')
          sfExport('betahat')
          sfExport('BetaHat')
          sfExport('betahatfms')
          sfExport('storeloss')
          sfExport('esterror')
          sfExport('GCV')
          sfExport('GCVf')
          sfExport('sig')
          sfExport('sigma')
          sfExport('sigmap')
          sfExport('df')
          
          
          ptm <- proc.time()
          Store <- sfSapply(1:100,simfunc)
          timecost <- proc.time()-ptm
          
          StoreLoss <- matrix(data = unlist(Store[1,]),ncol = md,byrow = T)
          EstError <- matrix(data = unlist(Store[2,]),ncol = md,byrow = T)
          w.record <- matrix(data = unlist(Store[3,]),ncol = MM[mi],byrow = T)
          w1.record <- matrix(data = unlist(Store[4,]),ncol = MM[mi],byrow = T)
          
          
          mseM[mi,,nsi,vari]<-apply(StoreLoss,2,mean)#MSE: M, methods, n_m, variance
          mseS[mi,,nsi,vari]<-apply(StoreLoss,2,sd)#MSE: M, methods, n_m, variance
          estM[mi,,nsi,vari]<-apply(EstError,2,mean)#est(beta): M, methods, n_m, variance
          estS[mi,,nsi,vari]<-apply(EstError,2,sd)#est(beta): M, methods, n_m, variance
          wM[mi,nsi,vari]<-list(w.record)
          w1M[mi,nsi,vari]<-list(w1.record)
          timecostM[mi,,nsi,vari]<-timecost[1:3]
          
          
          print(paste0("sitesamplesize=",n1," sitenumbers=",M," variance=",vari))
          print(timecost)
          #the best way to avoid crush: save your data!
          save.image(file = paste("pic_M_varies_CSL123_","linear_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",numcpu,".Rdata",sep=""))    }
      }
    }
    
    
    
    mse<-c()#organize the mse table
    mse.sd <- c()
    est<-c()
    est.sd<-c()
    timecost <- c()
    
    for (vari in 1:length(vaR)) {
      for (nsi in 1:length(ns)) {
        mse<-rbind(mse,mseM[,,nsi,vari])
        mse.sd<-rbind(mse.sd,mseS[,,nsi,vari])
        est<-rbind(est,estM[,,nsi,vari])
        est.sd<-rbind(est.sd,estS[,,nsi,vari])
        timecost<-rbind(timecost,timecostM[,,nsi,vari])
      }
    }
    colnames(mse)<-msename
    colnames(mse.sd)<-msename
    colnames(est)<-msename
    colnames(est.sd)<-msename
    colnames(timecost)<-c('user.all','system.all','elapsed.all')
    
    cvaR<-rep(vaR,each=length(MM)*length(ns))
    cnm<-rep(rep(ns,each=length(MM)),length(vaR))
    cM<-rep(MM,length(vaR)*length(ns))
    
    mse<-cbind(cvaR,cnm,cM,mse,mse.sd,est,est.sd)
    tim<-cbind(cvaR,cnm,cM,timecost)
    
    
    write.csv(mse,paste("pic_M_varies_CSL123_","linear_mse_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",numcpu,".csv",sep=""),row.names=T)
    write.csv(tim,paste("pic_M_varies_CSL123_","linear_time_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",numcpu,".csv",sep=""),row.names=T)
    
    
  }
}
sfStop()

#geom_rug

#estVar Calculate the Y-covariance matrix of multiple response regression
#SSD Calculate the residual sum of squares and residual inner product of Y in multi response regression

