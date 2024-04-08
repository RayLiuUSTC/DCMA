rm( list = ls ( all = TRUE))#remove all variables
#options(warn = 0)
setwd("~/lasso new center opt")

library(snowfall)

md=6+100 #five methods:PD,DC_{ew},Full,DC_{opt},CSL_{cv} and CSL_{oracle}
MM<-c(10,100) #the setting for number of sites
vaR<-c(0.5,1,123) #the setting for variance across different sites
ns<-c(50,100,1000) #the setting for local sample size
mseM <-array(0,c(length(MM),md,length(ns),length(vaR))) #MSE:site, method, local sample size and variance
timename <- c("p.max","EMA.max","f","OMA_solnp.max","OMA_osqp.max",
              "CSLcv.max","CSLall.max","CSLall_perlambda.max",
              'lasso_est_per_site.max', 'em_df_sigma_per_site.max','center_agg',
              'center_opt_solnp','center_opt_osqp','pred_OMA_solnp','pred_OMA_osqp',
              'center_opt_solnp.round','center_opt_osqp.round')
timeM <- array(0,c(length(MM),length(timename),length(ns),length(vaR)),dimnames = list(NULL,timename,NULL,NULL)) #Time Record:site, method, local sample size and variance
timecostM <-array(0,c(length(MM),3,length(ns),length(vaR)))
wM <- list()
length(wM) <- length(MM)*length(ns)*length(vaR)
dim(wM) <- c(length(MM),length(ns),length(vaR)) #weights(solnp, more sparse): each element of the array store a list of the weights ditributed to M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m
w1M <- list()
length(w1M) <- length(MM)*length(ns)*length(vaR)
dim(w1M) <- c(length(MM),length(ns),length(vaR)) #weights(osqp, more dense, coincide with the variance setting): each element of the array store a list of the weights ditributed to M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m
dfM <- list()
length(dfM) <- length(MM)*length(ns)*length(vaR)
dim(dfM) <- c(length(MM),length(ns),length(vaR)) #dfs: each element of the array store a list of the df estimator of the prediction model of M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m#MSE:站点数,存一张列表(模拟G次的MM个站点的自由度分布),站点样本量,方差
sigmaM <- list()
length(sigmaM) <- length(MM)*length(ns)*length(vaR)
dim(sigmaM) <- c(length(MM),length(ns),length(vaR)) #sigmas: each element of the array store a list of the estimates of variance of M sites(G=100 replicates) with loacl sample size n_m and error variance \sigma^2_m#MSE:站点数,存一张列表(模拟G次的MM个站点的方差估计分布),站点样本量,方差

by=0 #difference parameter, used for 'seq' function, seq(from,by,length.out), to set local sample size and allow it different across different site
p0<-50 #the number of potential covariates
pM<-c(15,30) #the setting for the number of covariates on hand
alphaM<-c(0.5,1,1.5) #the setting for the decay rate of coefficient
G<-100 #the setting for replicates
#delta <- 0.01#the setting for perturbation in derivative computing (computing df according to the definition)求导数近似间距
cpunum <- 25 #the setting for the amount of CPUs requested for the cluster


simfunc <- function(g){
  set.seed(1234+g)
  
  XX<-list() #save the on hand design matrix at each site
  YY<-list()
  MU<-list()
  emm<-list()
  
  time.lasso_est_per_site <- vector(length=M)
  
  for (j in 1:M){#DGP for M sites
    nm<-nn[j]#the sample size of jth site
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0) #generate the covariates at the jth site
    Xm<-X0m[,1:p] #take the first p covariates as the ones on hand
    vm<-rnorm(nm,0,sig[j]) #generate the errors of the jth site
    mum<-X0m%*%beta0 #the true mean vector of the jth site
    Ym<-mum+vm #the response Y_m
    
    ptm <- proc.time()
    lasso.mod2 = glmnet(Xm, Ym, alpha=1,intercept=F) #here compute lasso estimates for the initial value of method CSL
    #we apply 'standardize=TRUE' here.
    #10 folds by default. note that 'glmnet' partition the samples randomly rather than divide from the start with equal length.
    cv.lasso2 = cv.glmnet(Xm, Ym, alpha=1,intercept=F)#here compute hyperparameter in lasso for constructing surrogate likelihood of method CSL
    betahat[,2,j]<-(coef(lasso.mod2, s=cv.lasso2$lambda.min,intercept=F)[-1,]) #store the lasso estimates of beta at jth site
    time.lasso_est_per_site[j] <- (proc.time()-ptm)[3]
    
    XX[[j]]<-   Xm#store the design matrix on hand at jth site into XX list
    YY[[j]]<-   Ym#store the response vector at jth site into YY list
    MU[[j]]<-   mum#store the mean vector at jth site into MU list
  }
  
  time.lasso_est_per_site.max <- max(time.lasso_est_per_site)
  
  time.em_df_sigma_per_site <- vector(length=M)
  for (j in 1:M){#validate the prediction model picked at the jth site by the data of other M-1 sites
    Xm<-XX[[j]]
    Ym<-YY[[j]]
    
    ptm <- proc.time()
    em<-matrix(0,nn[j],M)#em is a matrix to record the residuals of M prediction models at the data of jth site
    for (i in 1:M){
      em[,i]<-Ym-Xm[,]%*%betahat[,2,i]
    }      
    emm[[j]]<-   t(em)%*%em#compute Q_l in step3
    df[j] <- sum(betahat[,2,j]!=0)#approximate the df with the non-zero elements of lasso estimates
    sigma[j]<-t(Ym-Xm[,]%*%betahat[,2,j])%*%(Ym-Xm[,]%*%betahat[,2,j])/(nn[j]-df[j])#variance estimate: the denominator is sample size minus df
    time.em_df_sigma_per_site[j] <- (proc.time()-ptm)[3]
  }
  time.em_df_sigma_per_site.max <- max(time.em_df_sigma_per_site)
  
  
  time.center_agg <- c()
  ptm <- proc.time()
  
  ee<-matrix(0,M,M)
  for (j in 1:M){
    ee<- ee + emm[[j]]#compute \sum_{m=1}^M Q_m in step4
  }
  a1 <- ee#step4:\sum_{m=1}^M Q_m
  a2<-sigma*df #compute \psai in step4
  w0 <- matrix(1,nrow=M,ncol=1)/M #the initial value of weights
  time.center_agg <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  ###########################################################################
  eq=function(w){sum(w)-1} #The equality constraint function
  obj = function(w){ t(w)%*%a1%*%w+2*t(w)%*%a2} #The optimaized objective:C(w)
  ######################################################################
  result1 <- solnp(matrix(1,nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result2 <- solnp(matrix(seq(from=0.01,by=0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result3 <- solnp(w0,obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result4 <- solnp(matrix(seq(from=0.01*M,by=-0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  bestobj <- which.min(c(tail(result1$values,n=1),tail(result2$values,n=1),tail(result3$values,n=1),tail(result4$values,n=1)))
  
  bestresult <- get(paste0("result",bestobj))
  w = bestresult$pars#note that in result, the argument 'UB' is not considered, probably because it is enough for the restraint that has been already considered
  #note that in result, the argument 'UB' is not considered, probably because it is enough for the restraint that has been already considered
  #trace=0: do not present the objective value and parameters of each iteration
  BetaHat=betahat[,2,]#extract lasso estimates of M local models from betahat to construct DC_{opt}
  betahatOMA_solnp=BetaHat%*%w;#DC_{opt}
  time.center_opt_solnp <- (proc.time()-ptm)[3]
  
  center_opt_solnp.iter.round <- result1$outer.iter+result2$outer.iter+result3$outer.iter+result4$outer.iter
  
  ptm <- proc.time()
  
  A=rbind(rep(1,M),diag(1,M))
  l=c(1,rep(0,M))
  u=c(1,rep(1,M))
  res <- solve_osqp(P=a1*2,q=a2*2,A=A,l=l,u=u,pars = osqpSettings(verbose = FALSE))
  w1 <- res$x
  w1[w1<0] <- 0
  w1 <- w1/sum(w1)
  BetaHat=betahat[,2,]#extract lasso estimates of M local models from betahat to construct DC_{opt}
  betahatOMA_osqp=BetaHat%*%w1;#the proposed estimator of beta: OPTW
  time.center_opt_osqp <- (proc.time()-ptm)[3]
  
  center_opt_osqp.iter.round <- res$info$iter
  
  
  
  time.center_ema <- c()
  ptm <- proc.time()
  betahatEMA=BetaHat%*%w0; #the proposed estimator of beta: OPTW
  time.center_ema <- (proc.time()-ptm)[3]
  
  time.full_est <- c()
  ptm <- proc.time()
  X<-c() 
  Y<-c() 
  if (vaR[vari]==123){#for the heteroscadatic situation, vari==123,sig<-sqrt(0.5+(1:M)/M)
    for (j in 1:M){
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }     
    #estimate the coefficients using lasso after pooling the data from all sites
    lasso.mod = glmnet(X, Y, alpha=1,intercept=F) 
    cv.lasso = cv.glmnet(X, Y, alpha=1,intercept=F)
    betahatf<-(coef(lasso.mod, s=cv.lasso$lambda.min,intercept=F)[-1,])
  } else { 
    for (j in 1:M){
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }
    #have the potential to improve the mse if considering the variance imformation: heteroscedastic adalasso
    #but we do not implement this temporarily
    lasso.mod = glmnet(X, Y, alpha=1,intercept=F) 
    cv.lasso = cv.glmnet(X, Y, alpha=1,intercept=F)
    betahatf<-(coef(lasso.mod, s=cv.lasso$lambda.min,intercept=F)[-1,])
  }
  time.full_est <- (proc.time()-ptm)[3]
  
  
  
  ######CSL method
  
  ptm <- proc.time()
  #######################
  gradient <- function(beta,Y,X){
    (-2*crossprod(X,Y)+2*crossprod(X)%*%beta)/length(Y)
  }
  #######################
  
  Centersite.index <- floor(M/2)
  betainitial <- betahat[,2,Centersite.index]#if you take betainitial2 into func 'gradient', you will find the column corresponding to centersite.index of the matrix gradM below will be 0!
  #betainitial2 <- apply(X=betahat[,2,],MARGIN = 1,FUN = mean)
  time.CSL_pre <- (proc.time()-ptm)[3]
  
  
  time.CSL_gradM <- vector(length = M)
  gradM <- matrix(nrow=p,ncol = M)
  for (j in 1:M) {
    ptm1 <- proc.time()
    gradM[,j] <- gradient(beta=betainitial,Y=YY[[j]],X=XX[[j]])
    time.CSL_gradM[j] <- (proc.time()-ptm1)[3]
  }
  time.CSL_gradM.max <- max(time.CSL_gradM)
  
  time.CSL_cv <- c()
  ptm2 <- proc.time()
  C <- seq(1e-04,1,length.out=100)
  s <- p
  pairCs <- expand.grid(C=C,s=s)
  N <- sum(nn)
  n <- n1
  lambdaproduce <-function(C,s){2*C*sqrt(log(p)/N)+2*C^2*s*(log(p)/sqrt(n*N)+log(p)/n)} 
  lambda_vals <- do.call(what = lambdaproduce,args = as.list(pairCs))
  ##########################
  CSLcv <- function(j){
    idxtrain <- idx!=j
    Ytrain <- YYc[idxtrain]
    Xtrain <- XXc[idxtrain,]
    idxtest <- idx==j
    Ytest <- YYc[idxtest]
    Xtest <- XXc[idxtest,]
    
    gradc <- gradient(beta = betainitial,Y=Ytrain,X=Xtrain)
    #only use the samples from central site to select hyperparameter, and do not implement cv on other site
    #thus when computing 'gradg' below, the sample size will differ across sites
    gradg <- (apply(X=gradM[,-Centersite.index],MARGIN = 1,FUN=sum)*n+gradc*(n-n_K))/(N-n_K)
    
    betatilde <- Variable(rows = p,cols = 1)
    #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
    #thus we use the formula 't(gradc-gradg)%*%betatilde' instead. 'Variable' (S4 type variable, CVXR)
    Ltilde <- sum((Ytrain-Xtrain%*%betatilde)^2)/(n-n_K) - t(gradc-gradg)%*%betatilde
    lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
    
    beta_vals <- sapply(lambda_vals,
                        function (lambda) {
                          obj <- Ltilde + lassoPenalty(betatilde, lambda)
                          #rm(Ltilde)
                          prob <- Problem(Minimize(obj))
                          rm(obj)
                          result <- solve(prob)
                          result$getValue(betatilde)
                        })
    ee_pred_cv <- Ytest - Xtest%*%beta_vals
    return(ee_pred_cv)
  }
  ########################
  
  ##CSLcv
  Kfold <- 10
  n_K <- floor(n/Kfold)
  idx <- kronecker(1:Kfold,rep(1,n_K))
  ee_pred_cv <- matrix(nrow = n,ncol = nrow(pairCs))#n=n_K*Kfold
  XXc <- XX[[Centersite.index]]
  YYc <- YY[[Centersite.index]]
  
  for (ii in 1:Kfold) {
    ee_pred_cv[idx==ii,] <- CSLcv(ii)
  }
  
  idx_lambda <- which.min(apply(ee_pred_cv,2,function(x){mean(x^2)}))
  if(idx_lambda==1|idx_lambda==100){warning("the range of lambda is small!")}
  cv.CSL <- lambda_vals[idx_lambda]
  time.CSL_cv <- (proc.time()-ptm2)[3]
  
  
  time.CSL_center_opt_est <- c()
  ptm <- proc.time()
  gradg <- apply(X=gradM,MARGIN = 1,FUN=mean)
  gradc <- gradM[,Centersite.index]
  
  betatilde <- Variable(rows = p,cols = 1)
  #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
  Ltilde <- sum((YYc-XXc%*%betatilde)^2)/n - t(gradc-gradg)%*%betatilde
  lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
  beta_val_cv <- sapply(cv.CSL,
                        function (lambda) {
                          obj <- Ltilde + lassoPenalty(betatilde, lambda)
                          prob <- Problem(Minimize(obj))
                          rm(obj)
                          result <- solve(prob)
                          result$getValue(betatilde)
                        })
  betatilde_cv <- beta_val_cv
  betatildeN_cv <- round(beta_val_cv,digits = 6)
  
  time.CSL_center_opt_est <- (proc.time()-ptm)[3]
  
  
  ############ all lambda
  
  time.CSL_all_lambda_beta_est <- c()
  ptm <- proc.time()
  betatilde <- Variable(rows = p,cols = 1)
  #note that when define 'Ltilde' below, function 'crossprod' can not be used because of the 'Variable' type of beta that we have defined
  Ltilde <- sum((YY[[Centersite.index]]-XX[[Centersite.index]]%*%betatilde)^2)/n - t(gradc-gradg)%*%betatilde
  lassoPenalty <- function(betatilde,lambda)lambda*p_norm(betatilde,1)
  beta_vals <- sapply(lambda_vals,
                      function (lambda) {
                        obj <- Ltilde + lassoPenalty(betatilde, lambda)
                        prob <- Problem(Minimize(obj))
                        result <- solve(prob)
                        result$getValue(betatilde)
                      })
  betatildeN <- round(beta_vals,digits = 6)
  time.CSL_all_lambda_beta_est <- (proc.time()-ptm)[3]
  time.CSL_all_lambda_beta_est.ave <- time.CSL_all_lambda_beta_est/length(lambda_vals)
  
  ############
  
  
  Yst<-c();Xst<-c();must<-c()
  for (j in 1:M){
    #generate the predictive data with the same size as the training data. Note that here the same variable names are used as the training part, thus the training data will be covered.
    nm<-nn[j]
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0)
    Xm<-X0m[,1:p]  
    vm<-rnorm(nm,0,sig[j])
    mum<-X0m%*%beta0
    Ym<-mum+vm
    Xst<-rbind(Xst,Xm)#pool the design matrices from all sites altogether
    must<-rbind(must,mum)
  }
  
  ptm <- proc.time()
  storeloss[1] <-mean((must-Xst%*%BetaHat[,M])^2)#prediction:PD
  time.pred_p <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[2] <-mean((must-Xst%*%betahatEMA)^2)#prediction:DC_{ew}
  time.pred_EMA <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[3] <-mean((must-Xst%*%betahatf)^2)#prediction:Full
  time.pred_f <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[4] <-mean((must-Xst%*%betahatOMA_solnp)^2)#prediction:DC_{opt}
  time.pred_OMA_solnp <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[5] <-mean((must-Xst%*%betahatOMA_osqp)^2)#prediction:DC_{opt}
  time.pred_OMA_osqp <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[6] <-mean((must-Xst%*%betatildeN_cv)^2)#prediction:CSL_{cv}
  time.pred_CSLcv <- (proc.time()-ptm)[3]
  
  ptm <- proc.time()
  storeloss[7:106] <-apply((matrix(must,nrow = length(must),ncol = 100)-Xst%*%betatildeN)^2,MARGIN = 2,FUN = mean)#prediction:CSL_{all}
  time.pred_CSLall <- (proc.time()-ptm)[3]
  time.pred_CSLall.ave <- time.pred_CSLall/length(lambda_vals)
  
  #summary of the computational time of all five methods
  
  time.p.max <- time.lasso_est_per_site[M] + time.pred_p
  time.EMA.max <- time.lasso_est_per_site.max + time.center_ema + time.pred_EMA
  time.f <- time.full_est + time.pred_f
  time.OMA.max <- time.lasso_est_per_site.max + time.em_df_sigma_per_site.max + time.center_agg + time.center_opt_solnp + time.pred_OMA_solnp
  time.OMA2.max <- time.lasso_est_per_site.max + time.em_df_sigma_per_site.max + time.center_agg + time.center_opt_osqp + time.pred_OMA_osqp
  time.CSL.max <- time.lasso_est_per_site[Centersite.index] + time.CSL_pre + time.CSL_gradM.max+time.CSL_cv + time.CSL_center_opt_est + time.pred_CSLcv
  time.CSLall.max <- time.lasso_est_per_site[Centersite.index] + time.CSL_pre + time.CSL_gradM.max + time.CSL_all_lambda_beta_est + time.pred_CSLall
  time.CSLallperlambda.max <- time.lasso_est_per_site[Centersite.index] + time.CSL_pre + time.CSL_gradM.max + time.CSL_all_lambda_beta_est.ave + time.pred_CSLall.ave
  
  #storelossb[g,1] <-mean((beta0-betahat[,1])^2)
  #storelossb[g,2] <-mean((beta0-betahatEMA)^2)
  #storelossb[g,3] <-mean((beta0-betahatf)^2)      
  #storelossb[g,4] <-mean((beta0-betahatMA)^2)
  
  timeTable <- c(time.p.max,time.EMA.max,time.f,time.OMA.max,time.OMA2.max,
                 time.CSL.max,time.CSLall.max,time.CSLallperlambda.max,
                 time.lasso_est_per_site.max, time.em_df_sigma_per_site.max, time.center_agg,
                 time.center_opt_solnp, time.center_opt_osqp, time.pred_OMA_solnp,time.pred_OMA_osqp,
                 center_opt_solnp.iter.round,center_opt_osqp.iter.round)
  
  return(list(storeloss=storeloss,timeTable=timeTable,
              w=w,w1=w1,df=df,sigma=sigma))
  #betahatP=BetaHat[,M],betahatEMA=betahatEMA,betahatf=betahatf,betahatOMA=betahatOMA,betatildeN_cv=betatildeN_cv,betatildeN=betatildeN))
}


sfInit(parallel = T,cpus = cpunum)

sfLibrary(glmnet)
sfLibrary(Rsolnp)
sfLibrary(CVXR)
sfLibrary(osqp)
sfExport('md')
sfExport('MM')
sfExport('vaR')
sfExport('ns')
sfExport('mseM')
sfExport('by')
sfExport('p0')
sfExport('G')

for (alphaind in 1:length(alphaM)) {
  alpha <- alphaM[alphaind]
  for (pind in 1:length(pM)) {
    p <- pM[pind]
    #generate the coefficients
    beta0<-sqrt(2*alpha)*c((1:p0)^(-alpha-0.5))#(runif(p0))
    
    sfExport('p')
    sfExport('alpha')
    sfExport('beta0')
    
    for (vari in 1:length(vaR)) {#start with \sigma^2=0.5
      for (nsi in 1:length(ns)) {
        n1=ns[nsi]#start with n_m=50
        for (mi in 1:length(MM)) {
          M<-MM[mi] #start with M=10
          nn<-seq(n1, by=by,length.out=M)#distribute sample size for every site
          
          betahat<-array(0,c(p,2,M))#c(p,2,M): the coefficients length of betahat; two kinds of ests - adalasso and lasso; the number of sites
          BetaHat<-array(0,c(p,M))
          #em<-matrix(0,n,M)
          storeloss <-vector("numeric",length=md)#prediciton loss for five methods
          #storelossb <-matrix(0,G,md)
          #note that 'sig' represents standard deviation which will be used in rnorm(mean=0,sd=sig)
          if (vaR[vari]==0.5){sig<-sqrt(0.5)*rep(1,M)} else#constant variance: \sigma^2=0.5
            if (vaR[vari]==1){sig<-1*rep(1,M)} else#constant variance: \sigma^2=1
            {sig<-sqrt(0.5+(1:M)/M)}#varied variance: \sigma^2=0.5+(1:M)/M
          #note that 'sigma' represents variance which will be used in criteria C(w)
          sigma<-matrix(0,M,1)#the matrix to store the variance ests of M local prediction models
          df<-matrix(0,M,1)#the matrix to store the df ests of M local prediction models
          
          
          sfExport('vari')
          sfExport('M')
          sfExport('nn')
          sfExport('n1')
          sfExport('betahat')
          sfExport('BetaHat')
          sfExport('storeloss')
          sfExport('sig')
          sfExport('sigma')
          sfExport('df')
          
          Store1 <- sfSapply(1:cpunum,simfunc)#warm start for 'proc.time'
          ptm <- proc.time()
          Store <- sfSapply(1:G,simfunc)
          timecost.whole <- proc.time()-ptm
          timecost.whole <- round(timecost.whole,3)
          
          
          StoreLoss <- matrix(data = unlist(Store[1,]),ncol = md,byrow = T)
          time.record <- matrix(data = unlist(Store[2,]),ncol = length(timename),byrow = T)
          w.record <- matrix(data = unlist(Store[3,]),ncol = MM[mi],byrow = T)
          w1.record <- matrix(data = unlist(Store[4,]),ncol = MM[mi],byrow = T)
          df.record <- matrix(data = unlist(Store[5,]),ncol = MM[mi],byrow = T)
          sigma.record <- matrix(data = unlist(Store[6,]),ncol = MM[mi],byrow = T)
          
          
          mseM[mi,,nsi,vari]<-apply(StoreLoss,2,mean)
          wM[mi,nsi,vari]<-list(w.record)
          w1M[mi,nsi,vari]<-list(w1.record)
          dfM[mi,nsi,vari]<-list(df.record)
          sigmaM[mi,nsi,vari]<-list(sigma.record)
          timeM[mi,,nsi,vari]<-apply(time.record,2,mean)
          timecostM[mi,,nsi,vari]<-timecost.whole[1:3]
          
          #apply(storeloss,2,mean)
          #msebM[mi,]<-apply(storelossb,2,mean)
          
          print(paste0("sitesamplesize=",n1," sitenumbers=",M," variance=",vari))
          print(paste0(names(timecost.whole[1]),"=",timecost.whole[1],",",names(timecost.whole[2]),"=",timecost.whole[2],",",names(timecost.whole[3]),"=",timecost.whole[3]))
          #the best way to avoid crush: save your data!
          save.image(file = paste("lasso_mse","_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",cpunum,".Rdata",sep=""))
        }
      }
    }
    #rownames(mseM)<-MM
    
    mse<-c()#organize the mse table
    tim<-c()#organize the time table
    timecost <- c()
    for (vari in 1:length(vaR)) {
      for (nsi in 1:length(ns)) {
        mse<-rbind(mse,mseM[,,nsi,vari])
        tim<-rbind(tim,timeM[,,nsi,vari])
        timecost<-rbind(timecost,timecostM[,,nsi,vari])
      }
    }
    colnames(mse)<-c('partial-data','equal-weight','full-data','optimal-weight.solnp','optimal-weight.osqp','CSL-cv',paste0("CSL-cen-lambda",1:100))
    colnames(tim)<-timename
    colnames(timecost)<-c('user.all','system.all','elapsed.all')
    
    cvaR<-rep(vaR,each=length(MM)*length(ns))
    cnm<-rep(rep(ns,each=length(MM)),length(vaR))
    cM<-rep(MM,length(vaR)*length(ns))
    
    mse<-cbind(cvaR,cnm,cM,mse)
    tim<-cbind(cvaR,cnm,cM,tim,timecost)
    
    
    write.csv(mse,paste("lasso_mse","_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",cpunum,".csv",sep=""),row.names=T)
    write.csv(tim,paste("lasso_time","_G_",G,"_p0_",p0,"_p_",p,"_alpha_",alpha,"_cpu_",cpunum,".csv",sep=""),row.names=T)
    #write.csv(msebM,paste("msebM_G_",G,"_p0_",p0,"_p_",p,"_N_",N,".csv",sep=""),row.names=T)
    
  }
}
sfStop()

