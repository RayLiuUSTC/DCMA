rm( list = ls ( all = TRUE))

library(tidyverse)
library(data.table)
library(wesanderson)

# the values of lambda
calpha <- c(0.5,1,1.5)
cp <- c(15,30)
vaR <- c(0.5,1,123)
ns <- c(50,100,1000)
MM <- c(10,100)
pair5 <- expand.grid(p=cp,n=ns,M=MM,vaR=vaR,alpha=calpha)
pair5dim1 <- nrow(pair5)
pair5dim2 <- ncol(pair5)
#lambdaname <- sapply(1:pair5dim1,function(i){paste(names(pair5),pair5[i,],sep = "=",collapse = ",")})

lambda_vals <- function(p,n,M,vaR,alpha){
C <- seq(1e-04,1,length.out=100)
s <- p
pairCs <- expand.grid(C=C,s=s)
N <- n * M
lambdaproduce <-function(C,s){2*C*sqrt(log(p)/N)+2*C^2*s*(log(p)/sqrt(n*N)+log(p)/n)} 
res <- do.call(what = lambdaproduce,args = as.list(pairCs))
return(res)
}

lambdaTable <- matrix(0,nrow = 100*pair5dim1,ncol = pair5dim2+1)
for (i in 1:pair5dim1) {
  lambdaTable[((i-1)*100+1):(i*100),1] <- do.call(lambda_vals,as.list(pair5[i,]))
  lambdaTable[((i-1)*100+1):(i*100),2:(pair5dim2+1)] <- matrix(rep(unlist(pair5[i,]),100),nrow = 100,ncol = pair5dim2,byrow = TRUE)
}
lambdaTable <- as.data.frame(lambdaTable)
colnames(lambdaTable) <- c('lambda',names(pair5))


# import data
p15a15 <- read.csv(file = "lasso_mse_G_100_p0_50_p_15_alpha_1.5_cpu_25.csv")


p15a15 <- p15a15[,-8]#delete solnp
colnames(p15a15)[5:9] <- c("PD","DC-ew","Full","DC-opt","CSL-cv")#change the names
p15a15$cvaR <- factor(p15a15$cvaR,levels = c(0.5,1,123),labels = c("0.5","1","0.5+m/M"))#change name for var(especially change "123" to "0.5+m/M")

p15a15CSL <- data.table(p15a15)
CSLname <- do.call(paste0,args = list("CSL.cen.lambda",1:100))
p15a15CSL.longformat <- melt(p15a15CSL,measure.vars = CSLname,id.vars = c("cvaR","cnm","cM"),
            variable.factor = TRUE,value.name = "mse")
colnames(p15a15CSL.longformat) <- c("vaR","n","M","lambda","CSL")
p15a15CSL.lambda <- subset(lambdaTable,p==15&alpha==1.5) %>% arrange(n,M,vaR)
p15a15CSL.longformat1 <- p15a15CSL.longformat %>% select(n,M,vaR,lambda,CSL)%>% arrange(n,M,vaR)
p15a15CSL.longformat1$lambda <- p15a15CSL.lambda$lambda
p15a15CSL.longformat2 <- matrix(0,nrow = nrow(p15a15CSL.longformat1),ncol = 5)
p15a15 <- p15a15 %>% arrange(cnm,cM,cvaR)

for (i in 1:nrow(p15a15)) {
  p15a15CSL.longformat2[((i-1)*100+1):(i*100),1] <- rep(p15a15$'PD'[i],100)
  p15a15CSL.longformat2[((i-1)*100+1):(i*100),2] <- rep(p15a15$'DC-ew'[i],100)
  p15a15CSL.longformat2[((i-1)*100+1):(i*100),3] <- rep(p15a15$'Full'[i],100)
  p15a15CSL.longformat2[((i-1)*100+1):(i*100),4] <- rep(p15a15$'DC-opt'[i],100)
  p15a15CSL.longformat2[((i-1)*100+1):(i*100),5] <- rep(p15a15$'CSL-cv'[i],100)
}
p15a15CSL.longformat2 <- as.data.frame(p15a15CSL.longformat2)
colnames(p15a15CSL.longformat2) <- names(select(p15a15,'PD':'CSL-cv'))

p15a15CSL.longformat3 <- cbind(p15a15CSL.longformat1,p15a15CSL.longformat2)
p15a15CSL.longformat4 <- p15a15CSL.longformat3 %>% filter(M==10)
p15a15CSL.longformat5 <- melt(data=p15a15CSL.longformat4,id.vars = 1:4,variable.name = "method",value.name = "mse")
p15a15CSL.longformat5$method <- factor(x=p15a15CSL.longformat5$method,levels = c('DC-opt','Full','CSL','CSL-cv','DC-ew','PD'))

base <- ggplot(p15a15CSL.longformat5,aes(x=lambda,y=mse,color=method,lty=method))+geom_line(linewidth=0.75)
label_false <- function(labels, multi_line = FALSE, sep ="=") {
  label_both(labels = labels, multi_line = multi_line, sep = sep)
}
base

labelname <- c(bquote(DC[opt]),bquote(Full),bquote(CSL['1,o']^(c)),bquote(CSL['1,cv']^(c)),bquote(DC[ew]),bquote(PD))
p15a15M10 <- base + facet_grid(vaR~n,scales = "free",labeller = label_false)+scale_x_log10()+scale_y_log10()+
  ylab("SPEs")+xlab(expression(lambda))+theme(strip.text = element_text(face = "bold",size = rel(0.8)),
                                              legend.position = c(1,1),legend.justification = c(1,1))+
  guides(col = guide_legend(nrow = 3, byrow = T))+
  scale_color_manual(values = wes_palette("Darjeeling1",6,type = "continuous"),name="method",labels=labelname)+
  scale_linetype_manual(values =c("solid", "dashed","longdash","dotdash","twodash","dotted") ,name="method",labels=labelname)+
  labs(title = expression(paste("SPEs versus ",lambda," in the penalized surrogate likelihood function of CSL " ,"for the penalized estimation situaion")),
       subtitle = expression(paste("When M=10, q=15 and ",alpha,"= 1.5")))+coord_cartesian(ylim = c(0.001,1))
p15a15M10

#ggsave("p15a15M10.png",p15a15M10,device = "png",height = 10,width = 13)
ggsave("p15a15M10(2).png",p15a15M10,device = "png",height = 10,width = 10)
ggsave("p15a15M10(2).pdf",p15a15M10,device = "pdf",height = 10,width = 10)









































