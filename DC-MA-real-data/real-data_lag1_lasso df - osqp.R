rm( list = ls ( all = TRUE))#remove all variables
#options(warn = 0)
setwd("~/DCMA real data new")

library(snowfall)
library(Rsolnp)
library(osqp)

###########################################################################
eq=function(w){sum(w)-1} #The equality constraint function
obj = function(w){ t(w)%*%a1%*%w+2*t(w)%*%a2}#The optimaized objective:C(w)

datacleanfit <- function(j){
  dataj<-data0[[j]]
  X0m<-dataj[,-c(1,2,3,4,5,6,15,18)]#No.1~5 represents serial number, year, day, and hour; No.6 is the response variable PM2.5; No.15 is the rainfall situation; No.18 is the site
  #Originally, wind direction is included as a dummy variable
  #Xm<-scale(as.matrix(cbind(X0m[,-9],model.matrix(~wd,data=(X0m))[,-1],X0m[,-9]^2,X0m[,-9]^3)))
  #Including 27 explanatory variables of polynomial order 1, 2, and 3, # No.9 is attribute variable about wind direction information
  Xm<-as.matrix(cbind(X0m[,-9],X0m[,-9]^2,X0m[,-9]^3))#extract the covariates
  Ym<-log(as.matrix(dataj[,c(6)]))#extract the response and perform log transformation
  datajlag1 <- data.frame(Ym[-1],Xm[-nrow(Xm),]) 
  datajlag1.c <- na.omit(datajlag1)#Exclude observations with missing values. '.c' represents 'complete'.
  
  Xm<-scale(datajlag1.c[,-1])#standardize the covariates
  Ym<-scale(datajlag1.c[,1])#standardize the log response
  XX<-   Xm#store the design matrix at the jth site
  YY<-   Ym#store the response at the jth site
  Xm<-Xm[ntr,]#find the training samples
  Ym<-Ym[ntr,]#find the training samples
  
  #Should intercept be added to fit OLS model when next the estimates will be the input of adalasso ?
  #lasso.mod3 = glmnet(Xm, Ym, alpha=1,intercept=F,penalty.factor = penalty.f2)
  
  penalty.f=abs(lm(Ym ~ Xm)$coefficients[-1])^(-1)
  lasso.mod = glmnet(Xm, Ym, alpha=1,intercept=F,penalty.factor = penalty.f)
  lasso.mod2 <- glmnet(Xm, Ym, alpha=1,intercept=F)
  cv.lasso = cv.glmnet(Xm, Ym, alpha=1,intercept=F,penalty.factor = penalty.f)
  cv.lasso2 = cv.glmnet(Xm, Ym, alpha=1,intercept=F)
  betahat<-(coef(lasso.mod, s=cv.lasso$lambda.min,intercept=F)[-1,])
  betahat2<-(coef(lasso.mod2, s=cv.lasso2$lambda.min,intercept=F)[-1,])
  ############################################
  #enet.mod = glmnet(Xm, Ym,intercept=F) 
  #cv.enet = cv.glmnet(Xm, Ym,intercept=F)
  #betahat[,2,j]<-coef(enet.mod, s=cv.enet$lambda.min,intercept=F)[-1,]
  
  
  nn<-length(Ym)#Record the total number of samples at the jth site
  return(list(XX=XX,YY=YY,betahat=betahat,nn=nn,betahat2=betahat2))
}

Qldfsigma <- function(j){
  Xm<-XX[[j]][ntr,]
  Ym<-YY[[j]][ntr,]
  em<-matrix(0,nn[j],M)#an error matrix to record each site's model estimating error of the data from the jth site
  for (i in 1:M){
    em[,i]<-Ym-Xm[,]%*%betahat[[i]]
  }      
  #YYm=matrix(Ym,nn[j],nn[j])+delta*diag(nn[j])#df_m is calculated by bruce force:by adding 0.001Δ's variation on y
  #dft<-matrix(0,nn[j],1)
  #for (t in 1:nn[j]) {
  #  lasso.mod = glmnet(Xm, YYm[,t], alpha=1,intercept=F,penalty.factor = Penalty.f[,j]) 
  #  betahatt<-coef(lasso.mod, s=CV.lasso[j])[-1,1]#Remove the first intercept term and keep the first one in lambda. min (usually only one)
  #  dft[t]<-(Xm[t,]%*%(betahatt-betahat[[j]]))/delta#Originally, Xm% *% (...) was supposed to be taken, and then compute the trace of it. But after directly multiplying Xm [t,], there is no need to take trace
  #}
  emm<-   t(em)%*%em#step3: computing Q_l
  #df<-apply(dft,2,sum)
  df <- sum(betahat[[j]]!=0)#sum(betahat[[j]]!=0)#Approximating degrees of freedom with non zero coefficients (only for calculating variance estimation)
  sigma<-crossprod(Ym-Xm[,]%*%betahat[[j]])/(nn[j]-df)#Variance estimation: denominator: sample size minus degrees of freedom
  return(list(emm=emm,df=df,sigma=sigma))
}

CSLcv <- function(j){
  idxtrain <- idx!=j
  Ytrain <- YYc[idxtrain]
  Xtrain <- XXc[idxtrain,]
  idxtest <- idx==j
  Ytest <- YYc[idxtest]
  Xtest <- XXc[idxtest,]
  
  gradc <- gradient(beta = betainitial,Y=Ytrain,X=Xtrain)
  gradg <- (apply(X=gradM[,-Centersite.index],MARGIN = 1,FUN=sum)*n+gradc*(n-n_K))/(N-n_K)
  
  betatilde <- Variable(rows = p,cols = 1)
  #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
  #thus we use the formula 't(gradc-gradg)%*%betatilde' instead. 'Variable' (S4 type variable, CVXR)
  Ltilde <- sum((Ytrain-Xtrain%*%betatilde)^2)/(n-n_K) - t(gradc-gradg)%*%betatilde
  lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
  
  beta_vals <- sapply(lambda_vals,
                      function (lambda) {
                        obj <- Ltilde + lassoPenalty(betatilde, lambda)
                        rm(Ltilde)
                        prob <- Problem(Minimize(obj))
                        rm(obj)
                        result <- solve(prob)
                        result$getValue(betatilde)
                      })
  
  #dd <- Xtest%*%beta_vals
  #View(sweep(-dd,1,Ytest,"+"))
  ee_pred_cv <- Ytest - Xtest%*%beta_vals
  return(ee_pred_cv)
}


gradient <- function(beta,Y,X){
  (-2*crossprod(X,Y)+2*crossprod(X)%*%beta)/length(Y)
}

#Ltilde = function(beta,gradc1=gradc,gradg1=gradg){crossprod(as.matrix(YY[[Centersite.index]])-as.matrix(XX[[Centersite.index]])%*%beta)/length(YY[[Centersite.index]]) - crossprod(gradc1-gradg1,beta)}


CSL_cv_lambda <- function (lambda) {
  obj <- Ltilde + lassoPenalty(betatilde, lambda)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  result$getValue(betatilde)
}

######################################################################

set.seed(1234)
#N<-2000
md=3+2+1 #five methods:PD,DC_{ew},Full,DC_{opt} and CSL
ns<-c(5000,10000,20000)#training sample size for each local site to ensure ns training samples per window period
wdn<-c(0,500,1000,1500,2000,2500,3000,3500,4000)#start point of sliding window
G<-1#replication of experiment
M=12#Number of stations, air quality monitoring stations
p=27#the number of covariates
#delta=0.01#perturbation setting
mseM <-array(0,c(length(ns),md,length(wdn)))#MSE:training sample size at each site, methods, start point of sliding window
betaM <-array(0,c(length(ns),md,p,length(wdn)))#beta:training sample size at each site, methods, parameter dimension, start point of sliding window
lambdaM <-array(0,c(length(ns),length(wdn)))#lambda:training sample size at each site, start point of sliding window
idx.lambdaM <-array(0,c(length(ns),length(wdn)))#lambda location:training sample size at each site, start point of sliding window
time.sf.saM <- array(0,c(length(ns),length(wdn)))#time.sf.sa:training sample size at each site, start point of sliding window

data0<-list()
data0[[1]]<-read.csv("PRSA_Data_Aotizhongxin_20130301-20170228.csv")
data0[[2]]<-read.csv("PRSA_Data_Changping_20130301-20170228.csv")
data0[[3]]<-read.csv("PRSA_Data_Dingling_20130301-20170228.csv")
data0[[4]]<-read.csv("PRSA_Data_Dongsi_20130301-20170228.csv")
data0[[5]]<-read.csv("PRSA_Data_Guanyuan_20130301-20170228.csv")
data0[[6]]<-read.csv("PRSA_Data_Gucheng_20130301-20170228.csv")
data0[[7]]<-read.csv("PRSA_Data_Huairou_20130301-20170228.csv")
data0[[8]]<-read.csv("PRSA_Data_Nongzhanguan_20130301-20170228.csv")
data0[[9]]<-read.csv("PRSA_Data_Shunyi_20130301-20170228.csv")
data0[[10]]<-read.csv("PRSA_Data_Tiantan_20130301-20170228.csv")
data0[[11]]<-read.csv("PRSA_Data_Wanliu_20130301-20170228.csv")
data0[[12]]<-read.csv("PRSA_Data_Wanshouxigong_20130301-20170228.csv")





#betahat<-list()
#betahatd<-list()
#BetaHat<-array(0,c(p,M))
storeloss <-matrix(0,G,md)#predictive loss: one experiment, methods
#sigma<-matrix(0,M,1)#standard deviation of M sites
#df<-matrix(0,M,1)#df of M sites
#nn<-matrix(0,M,1)#sample size of M sites




for (wdi in 1:length(wdn)){
  wd=wdn[wdi]# start from the first start point 0 of sliding window 
  for(nsi in 1:length(ns)){
    ni=ns[nsi]#start from the n_m=5000 
    
    
    ntr<-(wd+1):(wd+ni)#the sub-index of training samples: sample(1:30000)[1:nt]
    #XX<-list()##list XX: store design marices of each site
    #YY<-list()
    #emm<-list()
    
    if(nsi==3){
      sfInit(parallel = T,cpus = 3)
    }
    if(nsi==2){
      sfInit(parallel = T,cpus = 6)
    }
    if(nsi==1){
      sfInit(parallel = T,cpus = 12)
    }
    sfLibrary(glmnet)
    sfExport('data0')
    sfExport('ntr')
    myresult <- sfLapply(1:M,datacleanfit)
    
    XX <- lapply(myresult,function(x){x$XX})
    YY <- lapply(myresult,function(x){x$YY})
    #Penalty.f <- sapply(myresult,function(x){x$penalty.f})
    #CV.lasso <- sapply(myresult,function(x){x$cv.lasso})
    betahat <- lapply(myresult,function(x){x$betahat})
    nn <- sapply(myresult,function(x){x$nn})
    betahat2 <- lapply(myresult,function(x){x$betahat2})
    
    rm(myresult)
    
    sfExport('XX','YY','nn','M','betahat')
    
    myresult2 <- sfLapply(1:M,Qldfsigma)
    
    sfRemove('ntr','data0','XX','YY','nn','M','betahat')
    
    sfStop()
    
    emm <- lapply(myresult2,function(x){x$emm})
    df <- sapply(myresult2,function(x){x$df})
    #df.forsigma <- sapply(myresult2,function(x){x$df.forsigma})
    sigma <- sapply(myresult2,function(x){x$sigma})
    
    rm(myresult2)
    
    
    if(nsi==3){
      sfInit(parallel = T,cpus = 3)
    }
    if(nsi==2){
      sfInit(parallel = T,cpus = 6)
    }
    if(nsi==1){
      sfInit(parallel = T,cpus = 12)
    }
    #CSL method
    Centersite.index <- floor(M/2)
    #l0betahat <- apply(betahat[,2,],MARGIN = 2,FUN = function(x){sum(x!=0)})
    #supp.ave <- mean(l0betahat)
    #betainitial <- betahat2[[Centersite.index]]#if you take betainitial2 into func 'gradient', you will find the column corresponding to centersite.index of the matrix gradM below will be 0!
    betainitial <- apply(matrix(unlist(betahat2),nrow = p),1,mean)
    
    C <- seq(1e-04,10,length.out=100)
    #C <- seq(1e-04,1,length.out=100)
    #C <- c(seq(1e-04,1,length.out=50),seq(1.1,10,length.out=50))
    s <- p
    pairCs <- expand.grid(C=C,s=s)
    N <- sum(nn)
    n <- ni
    lambdaproduce <-function(C,s){2*C*sqrt(log(p)/N)+2*C^2*s*(log(p)/sqrt(n*N)+log(p)/n)} 
    lambda_vals <- do.call(what = lambdaproduce,args = as.list(pairCs))
    #View(matrix(lambda_vals))
    #vals_range <- range(lambda_vals)
    #lambda_vals2 <- c(seq(vals_range[1],vals_range[2],length.out=20),vals_range[1]+seq(0.01,0.1,0.01))
    YYg <- lapply(YY,function(x){x[ntr]})
    XXg <- lapply(XX,function(x){x[ntr,]})
    gradM <- mapply(FUN=gradient,Y=YYg,X=XXg,MoreArgs = list(beta=betainitial))
    
    beta_vals <- matrix(nrow = p,ncol = nrow(pairCs))
    
    Kfold <- 10
    n_K <- floor(n/Kfold)
    idx <- kronecker(1:Kfold,rep(1,n_K))
    #beta_vals_cv <- array(dim = c(p,Kfold,nrow(pairCs)))
    #beta_vals_cv <- matrix(nrow = p,ncol = Kfold)
    ee_pred_cv <- matrix(nrow = n,ncol = nrow(pairCs))##n=n_K*Kfold
    XXc <- XXg[[Centersite.index]]
    YYc <- YYg[[Centersite.index]]
    
    #it's time-consuming to  use 'sapply' and 'sfLapply' simultaneously
    sfExport('gradient','Centersite.index','betainitial','idx','XXc','YYc','gradM','n','n_K','N','lambda_vals','p','ntr')
    sfLibrary(CVXR)
    
    ptm <- proc.time()
    myresult3 <- sfLapply(1:Kfold,CSLcv)#18G,93%cpu in average
    time.sf.sa <- proc.time()-ptm#110.23s
    sfStop()
    ee_pred_cv <- Reduce(rbind,myresult3)
    
    # ptm <- proc.time()
    # sfLibrary(CVXR)
    # for (ii in 1:Kfold) {
    #   #there will be a problem if without [ntr,](the samples need to be copied to each slave of 'snowfall')
    #   for example, in Problem(Minimize(obj)), produce an object of 6.5GB(the matrix size is 7*31708^2/1024^3)
    #   idxtrain <- idx!=ii
    #   Ytrain <- YY[[Centersite.index]][ntr,][idxtrain]
    #   Xtrain <- XX[[Centersite.index]][ntr,][idxtrain,]
    #   idxtest <- idx==ii
    #   Ytest <- YY[[Centersite.index]][ntr,][idxtest]
    #   Xtest <- XX[[Centersite.index]][ntr,][idxtest,]
    #   
    #   gradc <- gradient(beta = betainitial,Y=Ytrain,X=Xtrain)
    #   gradg <- (apply(X=gradM[,-Centersite.index],MARGIN = 1,FUN=sum)*n+gradc*(n-n_K))/(N-n_K)
    #   
    #   betatilde <- Variable(rows = p,cols = 1)
    #   #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
    #   #Ltilde1 <- function(betatilde){sum((Ytrain-Xtrain%*%betatilde)^2)/(n-n_K) - t(gradc-gradg)%*%betatilde}
    #   Ltilde <- sum((Ytrain-Xtrain%*%betatilde)^2)/(n-n_K) - t(gradc-gradg)%*%betatilde
    #   lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
    #   
    #   #sfExport('gradient','Centersite.index','betainitial','idx','XX','YY','gradM','n','n_K','N','lambda_vals','p')
    #   sfExport('Ltilde')
    #   sfExport('lassoPenalty')
    #   sfExport('betatilde')
    #   
    #   beta_vals <- sfSapply(lambda_vals,CSL_cv_lambda)
    #   
    #   #slightly slow
    #   #beta_vals <- sapply(lambda_vals,
    #   #                    function (lambda) {
    #   #                      obj <- Ltilde + lassoPenalty(betatilde, lambda)
    #   #                      prob <- Problem(Minimize(obj))
    #   #                      result <- solve(prob)
    #   #                      result$getValue(betatilde)
    #   #                    })
    #   
    #   #dd <- Xtest%*%beta_vals
    #   #View(sweep(-dd,1,Ytest,"+"))
    #   ee_pred_cv[idxtest,] <- Ytest - Xtest%*%beta_vals
    # }#18-22G,37%-100%cpu
    # time.for.sf <- proc.time()-ptm#167.72s
    
    idx_lambda <- which.min(apply(ee_pred_cv,2,function(x){mean(x^2)}))
    if(idx_lambda==1|idx_lambda==100){warning("the range of lambda is small!")}
    cv.CSL <- lambda_vals[idx_lambda]
    
    gradg <- apply(X=gradM,MARGIN = 1,FUN=mean)
    gradc <- gradM[,Centersite.index]
    
    betatilde <- Variable(rows = p,cols = 1)
    #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
    Ltilde <- sum((YYc-XXc%*%betatilde)^2)/n - t(gradc-gradg)%*%betatilde
    lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
    beta_vals <- sapply(cv.CSL,
                        function (lambda) {
                          obj <- Ltilde + lassoPenalty(betatilde, lambda)
                          prob <- Problem(Minimize(obj))
                          rm(obj)
                          result <- solve(prob)
                          result$getValue(betatilde)
                        })
    betatilde <- beta_vals
    betatildeN <- round(beta_vals,digits = 6)
    
    sfStop()
    
    #DCMA method
    
    
    ee <- Reduce("+",emm)
    
    a1 <- ee#step4:\sum_{m=1}^M Q_m
    a2<-sigma*df #step4: computing \psai
    w0 <- matrix(1,nrow=M,ncol=1)/M #%
    
    result1 <- solnp(matrix(1,nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
    result2 <- solnp(matrix(seq(from=0.01,by=0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
    result3 <- solnp(w0,obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
    result4 <- solnp(matrix(seq(from=0.01*M,by=-0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
    bestobj <- which.min(c(tail(result1$values,n=1),tail(result2$values,n=1),tail(result3$values,n=1),tail(result4$values,n=1)))
    
    bestresult <- get(paste0("result",bestobj))
    w = bestresult$pars
    
    
    A=rbind(rep(1,M),diag(1,M))
    l=c(1,rep(0,M))
    u=c(1,rep(1,M))
    res <- solve_osqp(P=a1*2,q=a2*2,A=A,l=l,u=u,pars = osqpSettings(verbose = FALSE))
    w1 <- res$x
    w1[w1<0] <- 0
    w1 <- w1/sum(w1)
    
    #pooling betahat altogehter
    BetaHat <- Reduce(cbind,betahat)
    
    betahatOMA_solnp=BetaHat%*%w;#the estimator of beta:DC_{opt}
    betahatOMA_osqp=BetaHat%*%w1;#the estimator of beta:DC_{opt}
    betahatEMA=BetaHat%*%w0;#the estimator of beta:DC_{ew}
    #betahatMA1=BetaHat%*%w1;
    
    
    X <- Reduce(rbind,lapply(XX,function(x){x[ntr,]}))
    Y <- Reduce(rbind,lapply(YY,function(x){as.matrix(x[ntr,])}))
    
    penalty.f=abs(lm(Y ~ X)$coefficients[-1])^(-1)
    lasso.mod = glmnet(X, Y, alpha=1,intercept=F,penalty.factor = penalty.f) 
    cv.lasso = cv.glmnet(X, Y, alpha=1,intercept=F,penalty.factor = penalty.f)
    betahatf<-(coef(lasso.mod, s=cv.lasso$lambda.min,intercept=F)[-1,])
    #betahatf<-solve(SXX)%*%SXY
    
    Xst <- Reduce(rbind,lapply(XX,function(x){x[-c(1:max(ntr)),]}))
    Yst <- Reduce(rbind,lapply(YY,function(x){as.matrix(x[-c(1:max(ntr)),])}))
    
    storeloss[1] <-mean((Yst-Xst%*%BetaHat[,M])^2)
    storeloss[2] <-mean((Yst-Xst%*%betahatEMA)^2)
    storeloss[3] <-mean((Yst-Xst%*%betahatf)^2)
    storeloss[4] <-mean((Yst-Xst%*%betahatOMA_solnp)^2)
    storeloss[5] <-mean((Yst-Xst%*%betahatOMA_osqp)^2)
    storeloss[6] <-mean((Yst-Xst%*%betatilde)^2)
    
    #Xst1 <- XX[[Centersite.index]][-c(1:max(ntr)),]
    #Yst1 <- YY[[Centersite.index]][-c(1:max(ntr)),]
    #storeloss52 <- mean((Yst1-Xst1%*%betatilde)^2)
    
    mseM[nsi,,wdi]<-storeloss#because G=1, it can be left out
    idx.lambdaM[nsi,wdi] <- idx_lambda
    lambdaM[nsi,wdi] <- cv.CSL
    time.sf.saM[nsi,wdi] <- time.sf.sa[3]
    betaM[nsi,1,,wdi] <- BetaHat[,M]
    betaM[nsi,2,,wdi] <- betahatEMA
    betaM[nsi,3,,wdi] <- betahatf
    betaM[nsi,4,,wdi] <- betahatOMA_solnp
    betaM[nsi,5,,wdi] <- betahatOMA_osqp
    betaM[nsi,6,,wdi] <- betatilde
    
    save.image(file = paste("adalasso_real_data",".Rdata",sep="")) 
    print(paste0("wdi=",wdi,",nsi=",nsi))
  }
}
sfStop()

mse<-c()
for (wdi in 1:length(wdn)) {
  mse<-rbind(mse,mseM[,,wdi])#Reorganize array to a matrix: start point of sliding window, 3D to 2D
}
colnames(mse)<-c("PD","DC-ew","Full","DC-opt_solnp","DC-opt_osqp","CSL")
#colnames(mse)<-c("PD","DC-ew","Full","DC-opt_osqp","CSL")
rownames(mse)<-rep(ns,length(wdn))
write.csv(mse,paste("real_data_mse",".csv",sep=""),row.names=T)
#write.csv(msebM,paste("msebM_G_",G,"_p0_",p0,"_p_",p,"_N_",N,".csv",sep=""),row.names=T)

mse <- mse[,-4]#delete the "DC-opt_solnp" column
#for (ni in 1:nn) {
lthick1=matrix(1,wdn,1)###line width
col1=c("orange","blue","red","black","purple")#c(2,4,6,1)##line color
ltype1=c(6,4,2,1,3)##line type
ylim1=c(0.1,2)
#layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(4,4,0.5))

figname <- paste("real_data_mse",".eps",sep="")
#savePlot(file=figname,type="eps",dev.cur())
postscript(file = figname)#"save the picture below to location ...", use with 'dev.off'
par(mfcol=c(2,3))#drawn in an 2-by-3 array on the device by columns (mfcol)
#pdf(file=figname,width = 6, height = 4 )

for (i in 1:length(ns)) {#ns:5000,10000,20000 Three training sample sizes for plotting
  mse1=mse[(rownames(mse)==ns[i]),]#e.g.pick the mse that satisfies ns=5000 out, render a cumulative predictive mse plot about wdn
  main1 = bquote(paste('n','=',.(ns[i])))###
  par(mar=c(3.8, 3.8, 1.9, 0.5))###
  matplot(wdn,mse1,type="l",main=main1,xlab=expression(T[0]),ylab="MSPE",lty=ltype1,lwd=lthick1,ylim=,xaxt='n',col=col1,cex.lab=1.3,mgp=c(2.2,1,0))
  axis(1,wdn)
  #if(ni==1) 
  legend('topleft',c("PD",expression(DC[ew]),"Full",expression(DC[opt]),"CSL"), lty=ltype1,col=col1,bg='white',cex=1.2,ncol=2)##cex controls the size of legend
  matplot(wdn,cumsum(as.data.frame(mse1)),type="l",main=main1,xlab=expression(T[0]),ylab="cMSPE",lty=ltype1,lwd=lthick1,ylim=,xaxt='n',col=col1,cex.lab=1.3,mgp=c(2.2,1,0))
  axis(1,wdn)
  legend('topleft',c("PD",expression(DC[ew]),"Full",expression(DC[opt]),"CSL"), lty=ltype1,col=col1,bg='white',cex=1.2,ncol=2)##cex controls the size of legend
}



dev.off()

Nlag1_no.omit <- Reduce(sum,lapply(YY,nrow))
N <- Reduce(sum,lapply(data0,nrow))
#save.image(file = paste("real_data_PM25_adplasso_mse2",".Rdata",sep=""))
