rm( list = ls ( all = TRUE))#清除所有变量
#options(warn = 0)
setwd("D:\\simulation ols")
#setwd("D:/Rcode/DMA/2022-01-07")

library(snowfall)

#N<-2000
md=3+3+3*4 #四种方法
MM<-c(10,100,150) #站点数量
vaR<-c(0.5,1,123) #方差取值
ns<-c(50,100,1000,10000) #局部站点样本数量
mseM <-array(0,c(length(MM),md,length(ns),length(vaR))) #MSE:站点数,方法,站点样本量,方差
timename <- c("p.ave","p.max","EMA.ave","EMA.max","f.large","f.aic","f.bic","OMA.ave","OMA.max",
              paste0("CSL.ave",1:3),paste0("CSL.max",1:3),paste0("CSL2.ave",1:3),paste0("CSL2.max",1:3),
              paste0("CSL3.ave",1:3),paste0("CSL3.max",1:3),paste0("CSL4.ave",1:3),paste0("CSL4.max",1:3))
timeM <- array(0,c(length(MM),2*3+3+3*4*2,length(ns),length(vaR)),dimnames = list(NULL,timename,NULL,NULL))
#msebM <-matrix(0,length(MM),md)
#n1=50
by=0 #间隔数,seq参数,seq(from,by,length.out)
p0<-50 #潜在协变量个数
pM<-c(15,20,30) #候选模型使用的变量个数(误设定)
alphaM<-c(0.5,1,1.5) #系数衰减速度
G<-100 #实验重复次数
numcpu <- 10

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

gauss.R <- function(R,tQY){
  Rncol <- ncol(R)
  beta <- vector(length = Rncol)
  for (i in 1:Rncol) {
    beta[Rncol-i+1] <- tQY[Rncol-i+1]/R[(Rncol-i+1),(Rncol-i+1)]
    tQY[1:(Rncol-i+1)] <- tQY[1:(Rncol-i+1)]-(R[,(Rncol-i+1)]*beta[Rncol-i+1])[1:(Rncol-i+1)]
  }
  return(beta)
}

simfunc <- function(g){
  set.seed(1234+g)
  #重复实验
  XX<-list()
  #XX0<-list()
  YY<-list()
  MU<-list()
  emm<-list()
  #X0[,]<-matrix(rnorm(N*(p0),0,1),N,p0)
  #X[,]<-X0[,1:p]  
  #e<-rnorm(N,0,1);
  #mu=X[,2]+(2*X[,3]-1)^2+sin(2*pi*X[,4])/(2-sin(2*pi*X[,4]))+(0.1*sin(2*pi*X[,5])+0.2*cos(2*pi*X[,5])+0.3*(sin(2*pi*X[,5]))^2+0.4*(cos(2*pi*X[,5]))^3+0.5*(sin(2*pi*X[,5]))^3)#+0.5*(sin(2*pi*U[,1]))^2#
  #mu<-X0%*%beta0;
  #Y<-as.matrix(mu)+e;
  #betahatf<-solve(t(X)%*%X)%*%t(X)%*%Y
  
  time.ols_sigma_df_est_per_site <- vector(length=M)
  for (j in 1:M){#分每个站点进行讨论(DGP)
    nm<-nn[j]#第M个站点的样本量
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0)#第M个站点协变量的生成,均为N(0,1)
    Xm<-X0m[,1:p]#只取前15个变量作为on hand变量
    vm<-rnorm(nm,0,sig[j])#第M个站点的误差生成,均为N(0,sig)
    mum<-X0m%*%beta0#第M个站点真实的均值向量
    Ym<-mum+vm#第M个站点的响应变量Y
    #betahat[,j]<-solve(t(Xm)%*%Xm)%*%t(Xm)%*%Ym
    ptm1 <- proc.time()
    for (i in 1:p){#对第m个站点的数据,估计相应的估计量(nested)
      ols <- RcppEigen::fastLmPure(X=as.matrix(Xm[,1:i]), y=Ym)#without intercept already
      #emi <- ols$residuals
      temp.info <- extractAIC.lm(ols)
      GCV[j,i] <- temp.info[2]
      betahat[,i,j][1:i] <- ols$coefficients
    }
    sigma[j] <- temp.info[3]
    bestGCV <- which.min(GCV[j,])
    df[j] <- bestGCV
    BetaHat[,j]<-betahat[,bestGCV,j]#将第j个站点中,最小GCV值对应的系数估计作为候选模型系数,存储到BetaHat矩阵
    time.ols_sigma_df_est_per_site[j] <- (proc.time()-ptm1)[3]
    
    XX[[j]]<-   Xm#将第j个站点的on hand设计矩阵储存到XX列表
    YY[[j]]<-   Ym#将第j个站点的生成响应变量储存到YY列表
    MU[[j]]<-   mum#将第j个站点的均值向量储存到MU列表
  }
  time.ols_sigma_df_est_per_site.max <- max(time.ols_sigma_df_est_per_site)
  time.ols_sigma_df_est_per_site.ave <- mean(time.ols_sigma_df_est_per_site)
  
  time.em_per_site <- vector(length=M)
  for (j in 1:M){#利用每个站点挑选出来的模型对其他站点的数据进行预测
    Xm<-XX[[j]]
    Ym<-YY[[j]]
    
    ptm2 <- proc.time()
    em<-matrix(0,nn[j],M)#创建一个误差矩阵,每个站点的模型对第m个站点的数据进行估计的误差
    
    em<-matrix(Ym,ncol = M,nrow = nn[j])-Xm%*%BetaHat#也可用BetaHat
    
    emm[[j]]<-   t(em)%*%em#计算step3:Q_l
    time.em_per_site[j] <- (proc.time()-ptm2)[3]
  }#为什么要绕一个大圈?帽子矩阵除以0.01又还原,去掉零算muhatii和补零算muhat两者有区别吗
  time.em_per_site.max <- max(time.em_per_site)
  time.em_per_site.ave <- mean(time.em_per_site)
  
  time.center_opt <- c()
  ptm3 <- proc.time()

  ee<-Reduce("+",emm)

  a1 <- ee#step4:\sum_{m=1}^M Q_m
  a2<-sigma*df #计算step4:\psai
  #a3 <- t(rbind(matrix(1,nrow=1,ncol=M),diag(M),-diag(M)))
  #a4 <- rbind(1,matrix(0,nrow=M,ncol=1),matrix(-1,nrow=M,ncol=1))#注意:这里有每个自变量的UB限制
  w0 <- matrix(1,nrow=M,ncol=1)/M #等权重向量
  #QP <- solve.QP(a1+10^(-10)*diag(dim(a1)[1]),a2,a3,a4,1)#防止a1不是正定的!
  # w <- QP$solution
  # w <- as.matrix(w)
  # w <- w*(w>0)
  # w <- w/sum(w)     
  
  ###########################################################################
  eq=function(w){sum(w)-1} ###约束条件
  obj = function(w){ crossprod(w,a1)%*%w+2*crossprod(w,a2)}#准则的优化目标
  
  ######################################################################
  
  result1 <- solnp(matrix(1,nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result2 <- solnp(matrix(seq(from=0.01,by=0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result3 <- solnp(w0,obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  result4 <- solnp(matrix(seq(from=1,by=-0.01,length.out=M),nrow=M,ncol=1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))
  bestobj <- which.min(c(tail(result1$values,n=1),tail(result2$values,n=1),tail(result3$values,n=1),tail(result4$values,n=1)))
  
  w = get(paste0("result",bestobj))$pars#未设置UB,应该是觉得限制已足够
  betahatOMA=BetaHat%*%w;#最终beta的估计:文章方法OPTW
  time.center_opt <- (proc.time()-ptm3)[3]
  
  #a2<-0*matrix(c(sigma*df),M,1) 
  
  #w1 = solnp(matrix(rep(1/M,M),M,1),obj,eqfun=eq,eqB=0,LB=matrix(rep(0,M),M,1),control=list(trace=0))$pars
  
  
  time.center_ema <- c()
  ptm4 <- proc.time()
  betahatEMA=BetaHat%*%w0;#beta的估计:等权重EW
  #betahatMA1=BetaHat%*%w1;
  time.center_ema <- (proc.time()-ptm4)[3]
  
  #full.largest
  X<-c() #matrix(0,sum(nn),p)
  Y<-c() #matrix(0,sum(nn),1)
  #mu<-c() #matrix(0,sum(nn),1)
  ptm5 <- proc.time()
  if (vaR[vari]==123){#M个站点变化方差,vari==123,sig<-sqrt(0.5+(1:M)/M)
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,sigma[j]^(-0.5)*XX[[j]])
      Y<-rbind(Y,sigma[j]^(-0.5)*YY[[j]])
    }
    olsf <- RcppEigen::fastLmPure(X=as.matrix(X), y=Y)#without intercept already
    betahatfmax<-olsf$coefficients
    }  else  { 
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }
    olsf <- RcppEigen::fastLmPure(X=as.matrix(X), y=Y)#without intercept already
    betahatfmax<-olsf$coefficients
    }
  time.full_est_largest <- (proc.time()-ptm5)[3]
  
  #SXX<-0
  #SXY<-0
  #if (vaR[vari]==123){#M个站点变化方差,vari==123,sig<-sqrt(0.5+(1:M)/M)
  #  for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
  #    #X<-rbind(X,XX[[j]])
  #    #Y<-rbind(Y,YY[[j]])
  #    #mu<-rbind(mu,MU[[j]])
  #    #利用最大模型对方差进行估计:sigma[j]<-t(Ym-Xm[,]%*%betahat[,p,j])%*%(Ym-Xm[,]%*%betahat[,p,j])/(nn[j]-p)
  #    #将第j个站点的on hand设计矩阵储存到XX列表:XX[[j]]<-Xm
  #    SXX<-SXX+sigma[j]^(-1)*crossprod(XX[[j]])#X'\Sigma^{-1} X
  #    SXY<-SXY+sigma[j]^(-1)*crossprod(XX[[j]],YY[[j]])#X'\Sigma^{-1} Y
  #  }
  #  SXXQR <- qr(SXX)
  #  SXXQR.Q <- qr.Q(SXXQR)
  #  SXXQR.R <- qr.R(SXXQR)
  #  betahatf1<-gauss.R(SXXQR.R,crossprod(SXXQR.Q,SXY))}  else { 
  #    SXX <- Reduce("+",lapply(XX,crossprod))
  #    SXY <- Reduce("+",Map(function(x,y){crossprod(x=x,y=y)},XX,YY))
  #    SXXQR <- qr(SXX)
  #    SXXQR.Q <- qr.Q(SXXQR)
  #    SXXQR.R <- qr.R(SXXQR)
  #    betahatf1<-gauss.R(SXXQR.R,crossprod(SXXQR.Q,SXY))}
  #time.full_est <- (proc.time()-ptm5)[3]
  
  #full.aic
  X<-c() #matrix(0,sum(nn),p)
  Y<-c() #matrix(0,sum(nn),1)
  #mu<-c() #matrix(0,sum(nn),1)
  ptm5 <- proc.time()
  if (vaR[vari]==123){#M个站点变化方差,vari==123,sig<-sqrt(0.5+(1:M)/M)
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,sigma[j]^(-0.5)*XX[[j]])
      Y<-rbind(Y,sigma[j]^(-0.5)*YY[[j]])
    }
    for (i in 1:p){#对第m个站点的数据,估计相应的估计量(nested)
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      #emi <- ols$residuals
      temp.infof <- extractAIC.lm(olsf)
      GCVf[1,i] <- temp.infof[2]#AIC准则
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestAICf <- which.min(GCVf[1,])
    betahatfAIC<-betahatfms[,bestAICf]
  }  else  { 
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }
    for (i in 1:p){#对第m个站点的数据,估计相应的估计量(nested)
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      #emi <- ols$residuals
      temp.infof <- extractAIC.lm(olsf)
      GCVf[1,i] <- temp.infof[2]#AIC准则
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestAICf <- which.min(GCVf[1,])
    betahatfAIC<-betahatfms[,bestAICf]
  }
  time.full_est_aic <- (proc.time()-ptm5)[3]
  
  #full.bic
  X<-c() #matrix(0,sum(nn),p)
  Y<-c() #matrix(0,sum(nn),1)
  #mu<-c() #matrix(0,sum(nn),1)
  ptm5 <- proc.time()
  if (vaR[vari]==123){#M个站点变化方差,vari==123,sig<-sqrt(0.5+(1:M)/M)
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,sigma[j]^(-0.5)*XX[[j]])
      Y<-rbind(Y,sigma[j]^(-0.5)*YY[[j]])
    }
    for (i in 1:p){#对第m个站点的数据,估计相应的估计量(nested)
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      #emi <- ols$residuals
      temp.infof2 <- extractAIC.lm(fit = olsf,k=log(length(Y)))
      GCVf[2,i] <- temp.infof2[2]#BIC准则
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestBICf <- which.min(GCVf[2,])
    betahatfBIC<-betahatfms[,bestBICf]
  }  else  { 
    for (j in 1:M){#不用for循环,数组变矩阵是否会加快速度?答:会的
      X<-rbind(X,XX[[j]])
      Y<-rbind(Y,YY[[j]])
    }
    for (i in 1:p){#对第m个站点的数据,估计相应的估计量(nested)
      olsf <- RcppEigen::fastLmPure(X=as.matrix(X[,1:i]), y=Y)#without intercept already
      #emi <- ols$residuals
      temp.infof2 <- extractAIC.lm(fit = olsf,k=log(length(Y)))
      GCVf[2,i] <- temp.infof2[2]#BIC准则
      betahatfms[,i][1:i] <- olsf$coefficients
    }
    bestBICf <- which.min(GCVf[2,])
    betahatfBIC<-betahatfms[,bestBICf]
  }
  time.full_est_bic <- (proc.time()-ptm5)[3]
  
  
  ptm6 <- proc.time()
  #########
  gradient <- function(beta,Y,X){
    (-2*crossprod(X,Y)+2*crossprod(X)%*%beta)/length(Y)
  }
  
  Ltilde = function(beta,gradc1=gradc,gradg1=gradg){crossprod(YY[[Centersite.index]]-XX[[Centersite.index]]%*%beta)/length(YY[[Centersite.index]]) - crossprod(gradc1-gradg1,beta)}
  #,center1=Center,YY1=Y,XX1=X
  #########
  time.CSL.func <- (proc.time()-ptm6)[3]
  
  ptm6 <- proc.time()
  Centersite.index <- M-1
  betainitial2 <- betahat[,p,Centersite.index]#centersite.index列会是0列,属于方案2
  time.CSL_betainit2 <- (proc.time()-ptm6)[3]
  
  ptm6 <- proc.time()
  betainitial <- apply(X=betahat[,p,],MARGIN = 1,FUN = mean)
  time.CSL_betainit <- (proc.time()-ptm6)[3]
  
  Iterative <- function(n,betainit,method=rep(1,n),Center=Centersite.index,X=XX,Y=YY){
    time.CSL_iter <- matrix(nrow=2,ncol=n)#时间max/ave;迭代次数
    betatildeN <- matrix(nrow = length(betainit),ncol = n)
    for (i in 1:n) {
      ptm7 <- proc.time()
      if(i>1){betainit=betatildeN[,(i-1)]}
      time.CSL_betainit_iter <- (proc.time()-ptm7)[3]
      
      time.gradM <- vector(length = M)
      gradM <- matrix(0,nrow = p,ncol = M)
      for (j in 1:M) {
        ptm8 <- proc.time()
        gradM[,j] <- gradient(beta=betainit,Y=Y[[j]],X=X[[j]])
        time.gradM[j] <- (proc.time()-ptm8)[3]
      } 
      time.gradM.max <- max(time.gradM)
      time.gradM.ave <- mean(time.gradM)
      
      ptm9 <- proc.time()
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
      time.CSL_center_opt <- (proc.time()-ptm9)[3]
      time.CSL_iter[1,i]=time.CSL_betainit_iter+time.gradM.max+time.CSL_center_opt#max time
      time.CSL_iter[2,i]=time.CSL_betainit_iter+time.gradM.ave+time.CSL_center_opt#ave time
      }
    return(list(betatildeN=betatildeN,time.CSL_iter=time.CSL_iter))
  }
  
  Iter <- Iterative(n=3,betainit = betainitial)#仿Jordan simulation
  betatildeN <- Iter$betatildeN
  time.CSL_est <- Iter$time.CSL_iter
  
  Iter2 <- Iterative(n=3,betainit = betainitial2)#centersite est作为初始估计
  betatildeN2 <- Iter2$betatildeN
  time.CSL_est2 <- Iter2$time.CSL_iter
  
  Iter3 <- Iterative(n=3,betainit = betainitial,method = c(2,2,2))#method 2(exact)
  betatildeN3 <- Iter3$betatildeN
  time.CSL_est3 <- Iter3$time.CSL_iter
  
  Iter4 <- Iterative(n=3,betainit = betainitial,method = c(2,1,1))#method 2,1,1(exact,appro,appro)
  betatildeN4 <- Iter4$betatildeN
  time.CSL_est4 <- Iter4$time.CSL_iter
  
  #betatildeNE <- Iterative(n=3,betainit = betainitial,method = c(1,1,2))
  #betatildeNE2 <- Iterative(n=3,betainit = betainitial,method = c(2,1,1))
  
  

  Yst<-c();Xst<-c();must<-c()
  for (j in 1:M){
    #生成预测集,与训练集同样大小,注意用了与训练集同样的变量(65-86行),训练集数据将被覆盖
    nm<-nn[j]
    X0m<-matrix(rnorm(nm*(p0),0,1),nm,p0)
    Xm<-X0m[,1:p]  
    vm<-rnorm(nm,0,sig[j])
    mum<-X0m%*%beta0
    Xst<-rbind(Xst,Xm)#所有站点的预测设计矩阵拼接成一个大矩阵
    must<-rbind(must,mum)
  }
  
  #g代表第g次重复实验,共G=100次
  ptm10 <- proc.time()
  storeloss[1] <-mean((must-Xst%*%BetaHat[,M])^2)#用最后一个站点M的模型做预测
  time.pred_p <- (proc.time()-ptm10)[3]
  
  ptm11 <- proc.time()
  storeloss[2] <-mean((must-Xst%*%betahatEMA)^2)#等权重
  time.pred_EMA <- (proc.time()-ptm11)[3]
  
  ptm12 <- proc.time()
  storeloss[3] <-mean((must-Xst%*%betahatfmax)^2)#全样本最大模型
  time.pred_fmax <- (proc.time()-ptm12)[3]
  
  ptm12 <- proc.time()
  storeloss[4] <-mean((must-Xst%*%betahatfAIC)^2)#全样本AIC
  time.pred_fAIC <- (proc.time()-ptm12)[3]
  
  ptm12 <- proc.time()
  storeloss[5] <-mean((must-Xst%*%betahatfBIC)^2)#全样本BIC
  time.pred_fBIC <- (proc.time()-ptm12)[3]
  
  ptm13 <- proc.time()
  storeloss[6] <-mean((must-Xst%*%betahatOMA)^2)#文章方法
  time.pred_OMA <- (proc.time()-ptm13)[3]

  time.pred_CSL <- vector(length = 3)#时间max/ave;迭代次数
  for (i in 1:3) {
    ptm14 <- proc.time()
    storeloss[6+i] <- mean((must-Xst%*%betatildeN[,i])^2)#CSL方法(ave init)
    time.pred_CSL[i] <- (proc.time()-ptm14)[3]
  }
  
  time.pred_CSL2 <- vector(length = 3)#时间max/ave;迭代次数
  for (i in 1:3) {
    ptm15 <- proc.time()
    storeloss[9+i] <- mean((must-Xst%*%betatildeN2[,i])^2)#CSL方法(centersite init)
    time.pred_CSL2[i] <- (proc.time()-ptm15)[3]
  }
  
  time.pred_CSL3 <- vector(length = 3)#时间max/ave;迭代次数
  for (i in 1:3) {
    ptm16 <- proc.time()
    storeloss[12+i] <- mean((must-Xst%*%betatildeN3[,i])^2)#CSL方法(ave init,method=2,2,2)
    time.pred_CSL3[i] <- (proc.time()-ptm16)[3]
  }
  
  time.pred_CSL4 <- vector(length = 3)#时间max/ave;迭代次数
  for (i in 1:3) {
    ptm17 <- proc.time()
    storeloss[15+i] <- mean((must-Xst%*%betatildeN4[,i])^2)#CSL方法(ave init,method=2,1,1)
    time.pred_CSL4[i] <- (proc.time()-ptm17)[3]
  }
  
  #计算时间汇总
  time.p.ave <- time.ols_sigma_df_est_per_site.ave+time.pred_p
  time.p.max <- time.ols_sigma_df_est_per_site.max+time.pred_p
  
  time.EMA.ave <- time.ols_sigma_df_est_per_site.ave+time.center_ema+time.pred_EMA
  time.EMA.max <- time.ols_sigma_df_est_per_site.max+time.center_ema+time.pred_EMA
  
  time.f.largest <- time.full_est_largest+time.pred_fmax
  time.f.AIC <- time.full_est_aic+time.pred_fAIC
  time.f.BIC <- time.full_est_bic+time.pred_fBIC
  
  time.OMA.ave <- time.ols_sigma_df_est_per_site.ave+time.em_per_site.ave+time.center_opt+time.pred_OMA
  time.OMA.max <- time.ols_sigma_df_est_per_site.max+time.em_per_site.max+time.center_opt+time.pred_OMA
  
  time.CSL.initave_ave.iterave <- time.ols_sigma_df_est_per_site.ave + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est[2,]) + time.pred_CSL #betainitial 是所有站点平均
  time.CSL.initave_max.itermax <- time.ols_sigma_df_est_per_site.max + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est[1,]) + time.pred_CSL
  
  time.CSL2.initcen.iterave <- time.ols_sigma_df_est_per_site[Centersite.index] + time.CSL.func + time.CSL_betainit2 + cumsum(time.CSL_est2[2,]) + time.pred_CSL2 #betainitial2 是中心站点初值
  time.CSL2.initcen.itermax <- time.ols_sigma_df_est_per_site[Centersite.index] + time.CSL.func + time.CSL_betainit2 + cumsum(time.CSL_est2[1,]) + time.pred_CSL2
  
  time.CSL3.initave_ave.iterave <- time.ols_sigma_df_est_per_site.ave + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est3[2,]) + time.pred_CSL3 #betainitial 是所有站点平均
  time.CSL3.initave_max.itermax <- time.ols_sigma_df_est_per_site.max + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est3[1,]) + time.pred_CSL3
  
  time.CSL4.initave_ave.iterave <- time.ols_sigma_df_est_per_site.ave + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est4[2,]) + time.pred_CSL4 #betainitial 是所有站点平均
  time.CSL4.initave_max.itermax <- time.ols_sigma_df_est_per_site.max + time.CSL.func + time.CSL_betainit + cumsum(time.CSL_est4[1,]) + time.pred_CSL4
  
  #timeTable <- c(time.p.max,time.EMA.max,time.f.AIC,time.OMA.max,
  #               time.CSL.initave_max.itermax,time.CSL2.initcen.itermax,
  #               time.CSL3.initave_max.itermax,time.CSL4.initave_max.itermax)
  
  
  #storelossb[g,1] <-mean((beta0-betahat[,1])^2)
  #storelossb[g,2] <-mean((beta0-betahatEMA)^2)
  #storelossb[g,3] <-mean((beta0-betahatf)^2)      
  #storelossb[g,4] <-mean((beta0-betahatMA)^2)      
  
  #sfCat(paste("Iteration ", g), sep = "\n")#实时记录sf并行状况
  #注意上述函数需sfLibrary("snowfall", character.only=TRUE)
  #且最好使用slaveOutfile参数:sfInit(parallel=TRUE, cpus=4, slaveOutfile="test.txt")
  
  return(list(storeloss=storeloss,
              time.p.ave=time.p.ave,time.p.max=time.p.max,
              time.EMA.ave=time.EMA.ave,time.EMA.max=time.EMA.max,
              time.f.largest=time.f.largest,time.f.AIC=time.f.AIC,time.f.BIC=time.f.BIC,
              time.OMA.ave=time.OMA.ave,time.OMA.max=time.OMA.max,
              time.CSL.initave_ave.iterave=time.CSL.initave_ave.iterave,time.CSL.initave_max.itermax=time.CSL.initave_max.itermax,
              time.CSL2.initcen.iterave=time.CSL2.initcen.iterave,time.CSL2.initcen.itermax=time.CSL2.initcen.itermax,
              time.CSL3.initave_ave.iterave=time.CSL3.initave_ave.iterave,time.CSL3.initave_max.itermax=time.CSL3.initave_max.itermax,
              time.CSL4.initave_ave.iterave=time.CSL4.initave_ave.iterave,time.CSL4.initave_max.itermax=time.CSL4.initave_max.itermax))
              #betahatOMA=betahatOMA,betahatEMA=betahatEMA,betahatfmax=betahatfmax,betahatfAIC=betahatfAIC,betahatfBIC=betahatfBIC,
              #betatilde=betatildeN[,1],betatilde2=betatildeN[,3],betatilde3=betatildeN2[,1],betatilde3=betatildeN2[,3]))
}




for (alphaind in 1:3) {
  alpha <- alphaM[alphaind]
  for (pind in 1:3) {
    p <- pM[pind]
    #系数大小
    beta0<-sqrt(2*alpha)*c((1:p0)^(-alpha-0.5))#(runif(p0))

    
    sfInit(parallel = T,cpus = numcpu)
    
    sfLibrary(Rsolnp)
    
    sfExport('md')
    sfExport('MM')
    sfExport('vaR')
    sfExport('ns')
    sfExport('mseM')
    sfExport('by')
    sfExport('p0')
    sfExport('p')
    sfExport('alpha')
    sfExport('beta0')
    sfExport('G')
    sfExport('extractAIC.lm')
    sfExport('deviance.lm')
    sfExport('gauss.R')
    
    
    
    for (vari in 1:length(vaR)) { #从方差等于0.5开始
      for (nsi in 1:length(ns)) {
        n1=ns[nsi] #从站点样本量等于50开始
        for (mi in 1:length(MM)) {
          M<-MM[mi] #从10个站点开始
          nn<-seq(n1, by=by,length.out=M) #为每个站点分配样本数
          
          
          betahat<-array(0,c(p,p,M))#betahat系数长度,站点内nested模型个数,站点数
          BetaHat<-array(0,c(p,M))#用AIC准则挑选的最优的nested模型的betahat系数,站点数
          betahatfms<-matrix(0,nrow=p,ncol=p)#betahat系数长度,站点内nested模型个数
          #em<-matrix(0,n,M)
          storeloss <-vector("numeric",length=md) #预测损失:方法
          #storelossb <-matrix(0,G,md)
          GCV<-matrix(0,M,p) #准则值:站点个数,on hand变量个数
          GCVf<-matrix(0,2,p)#准则值:AIC及BIC两种准则,on hand变量个数
          if (vaR[vari]==0.5){sig<-sqrt(0.5)*rep(1,M)} else #M个站点常值方差0.5
          if(vaR[vari]==1){sig<-1*rep(1,M)} else #M个站点常值方差1
          {sig<-sqrt(0.5+(1:M)/M)}#M个站点变化方差,vari==123
          sigma<-matrix(0,M,1)#M个站点的标准差向量
          df<-matrix(0,M,1)#M个站点的自由度向量
          
          sfExport('vari')
          sfExport('M')
          sfExport('nn')
          sfExport('betahat')
          sfExport('BetaHat')
          sfExport('betahatfms')
          sfExport('storeloss')
          sfExport('GCV')
          sfExport('GCVf')
          sfExport('sig')
          sfExport('sigma')
          sfExport('df')
          
          ptm <- proc.time()
          Store <- sfSapply(1:100,simfunc)
          timecost <- proc.time()-ptm
          
          StoreLoss <- matrix(data = unlist(Store[1,]),ncol = md,byrow = T)
          time.record <- matrix(data = unlist(Store[2:(1+2*3+3+4*2),]),ncol = length(timename),byrow = T)
          #set.seed(20220825)
          #storetable <- sfLapply(1:100,simfunc)
          #storeloss <- matrix(unlist(storetable),ncol=4,byrow = T)
          
          mseM[mi,,nsi,vari]<-apply(StoreLoss,2,mean)#MSE:站点数,方法,站点样本量,方差
          timeM[mi,,nsi,vari]<-apply(time.record,2,mean)
          #apply(storeloss,2,mean)
          #msebM[mi,]<-apply(storelossb,2,mean)
          
          print(paste0("sitesamplesize=",n1," sitenumbers=",M," variance=",vari))
          print(timecost)
          #防机子炸了!
          save.image(file = paste("quick_ver_1_","linear_mse_G_",G,"_p0_",p0,"_p_",p,"_by_",by,"_alpha_",alpha,"_cpu_",numcpu,".Rdata",sep=""))    }
      }
    }
    #rownames(mseM)<-MM
    sfStop()
    
    mse<-c()#归纳整合表格
    tim<-c()
    for (vari in 1:length(vaR)) {
      for (nsi in 1:length(ns)) {
        mse<-rbind(mse,mseM[,,nsi,vari])
        tim<-rbind(tim,timeM[,,nsi,vari])
      }
    }
    colnames(mse)<-c('partial-data','equal-weight','full-data.large','full-data.AIC','full-data.BIC','optimal-weight',paste0("CSL-ave-a3",1:3),paste0("CSL-cen-a3",1:3),paste0("CSL-ave-e3",1:3),paste0("CSL-ave-eaa",1:3))
    colnames(tim)<-timename
    
    cvaR<-rep(vaR,each=length(MM)*length(ns))
    cnm<-rep(rep(ns,each=length(MM)),length(vaR))
    cM<-rep(MM,length(vaR)*length(ns))
    
    #mse<-cbind(cvaR,cnm,cM,mse,timecost,timecost.DGP.sitefit,timecost.Ql,timecost.w,timecost.fullbeta.est,timecost.pred)
    mse<-cbind(cvaR,cnm,cM,mse)
    tim<-cbind(cvaR,cnm,cM,tim)
    
    
    
    write.csv(mse,paste("quick_ver_1_","linear_mse_G_",G,"_p0_",p0,"_p_",p,"_by_",by,"_alpha_",alpha,"_cpu_",numcpu,".csv",sep=""),row.names=T)
    write.csv(tim,paste("quick_ver_1_","linear_time_G_",G,"_p0_",p0,"_p_",p,"_by_",by,"_alpha_",alpha,"_cpu_",numcpu,".csv",sep=""),row.names=T)
    #write.csv(msebM,paste("msebM_G_",G,"_p0_",p0,"_p_",p,"_N_",N,".csv",sep=""),row.names=T)
    
  }
}


View(StoreLoss)
View(GCVf)



