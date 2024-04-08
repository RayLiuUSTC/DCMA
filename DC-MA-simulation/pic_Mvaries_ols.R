rm( list = ls ( all = TRUE))

# import data
load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_15_alpha_0.5_cpu_50.Rdata")
#load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_15_alpha_1.5_cpu_34.Rdata")
#load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_15_alpha_1_cpu_34.Rdata")
#load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_30_alpha_0.5_cpu_12.Rdata")
#load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_30_alpha_1_cpu_34.Rdata")
#load("C:\\Users\\86182\\OneDrive\\desktop\\pic_M_varies_CSL123_linear_G_100_p0_50_p_30_alpha_1.5_cpu_34.Rdata")



library(tidyverse)
library(data.table)
library(wesanderson)
#tablename <- c('PD','DC-ew','Full','DC-opt-solnp','DC-opt-osqp','CSLa1','CSLa2','CSLa3','CSLc1','CSLc2','CSLc3','n','vaR','M')
tablename <- c('PD','DC-ew','Full','DC-opt-osqp','CSLa1','n','vaR','M')
#tablename <- c(expression(PD),expression(DC[ew]),expression(Full),expression(DC[opt]),expression(CSL[1]^{(a)}),'n','vaR','M')
mseMo <- mseM
mseM <- mseM[,c(1:3,5,6),,]
mseSo <- mseS
mseS <- mseS[,c(1:3,5,6),,]

mseM.table <- matrix(data = NA,nrow = length(MM)*length(ns)*length(vaR),ncol = dim(mseM)[2]+3)
for (i in 1:length(ns)) {
  for (j in 1:length(vaR)) {
    mseM.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),1:dim(mseM)[2]] <- mseM[,,i,j]
    mseM.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseM)[2]+1] <- ns[i]
    mseM.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseM)[2]+2] <- vaR[j]
    mseM.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseM)[2]+3] <- MM
  }
}
colnames(mseM.table) <- tablename

mseS.table <- matrix(data = NA,nrow = length(MM)*length(ns)*length(vaR),ncol = dim(mseM)[2]+3)
for (i in 1:length(ns)) {
  for (j in 1:length(vaR)) {
    mseS.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),1:dim(mseS)[2]] <- mseS[,,i,j]
    mseS.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseS)[2]+1] <- ns[i]
    mseS.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseS)[2]+2] <- vaR[j]
    mseS.table[(((i-1)*length(ns)+j-1)*length(MM)+1):(((i-1)*length(ns)+j)*length(MM)),dim(mseS)[2]+3] <- MM
  }
}
colnames(mseS.table) <- tablename


#wide to long form
mseM.table <- data.table(mseM.table)
mseS.table <- data.table(mseS.table)
mseM.table$vaR <- factor(mseM.table$vaR,levels = c(0.5,1,123),labels = c("0.5","1","0.5+m/M"))#change name for var(especially change "123" to "0.5+m/M")
mseS.table$vaR <- factor(mseS.table$vaR,levels = c(0.5,1,123),labels = c("0.5","1","0.5+m/M"))#change name for var(especially change "123" to "0.5+m/M")

mseM.long <- melt(mseM.table,measure.vars = 1:dim(mseM)[2],id.vars = c("n","vaR",'M'),
                  variable.factor = TRUE,value.name = "mseM",variable.name = "method")
mseS.long <- melt(mseS.table,measure.vars = 1:dim(mseM)[2],id.vars = c("n","vaR",'M'),
                  variable.factor = TRUE,value.name = "mseS",variable.name = "method")
mse.long <- data.table(mseM.long,mseS=mseS.long$mseS)
mse.long <- mutate(mse.long,lower=if_else(mseM-mseS<0,mseM,mseM-mseS),upper=mseM+mseS) 
mse.long$method <- factor(x=mse.long$method,levels = c('DC-opt-osqp','Full','CSLa1','DC-ew','PD'))

#pic

#cols<-c("#982b2b","#0074b3",
#        "#f47720","#459943","#e5ce81")
##update_geom_defaults("point", list(shape = c(1,2,3,5,6)))
#base <- ggplot(mse.long,aes(x=M,y=mseM,color=method,lty=method,shape=method))
#label_false <- function(labels, multi_line = FALSE, sep ="=") {
#  label_both(labels = labels, multi_line = multi_line, sep = sep)
#}
#Mvaries <- base + facet_grid(n~vaR,scales = "free")+scale_x_log10()+scale_y_log10()+
#  ylab("SPEs")+geom_rug(sides = "b",color="grey",length = unit(0.8, "mm"),linewidth=1)+
#  geom_errorbar(aes(ymin=lower,ymax=upper),lty=1,alpha=0.7)+
#  labs(title = expression(paste("SPEs versus number of sites M for the least squares estimation situation")),
#       subtitle = expression(paste("when q=20 and ",alpha,"=1.5")))+
#  geom_line(aes(x=M,y=mseM),linewidth=0.75,alpha=0.75)+geom_point()+
#  scale_discrete_manual(aesthetics = c("color","shape","linetype"),name="method",breaks=c('DC-opt-osqp','Full','CSLa1','DC-ew','PD'),
#                        labels=c(expression(DC[opt]),"Full",expression(CSL[1]^{(a)}),expression(DC[ew]),"PD"),
#                        values = list(color = c("#982b2b","#0074b3","#f47720","#459943","#e5ce81"), linetype = c("solid", "dashed","dotdash","longdash","twodash"), shape = c(16,15,17,1,2)))
#Mvaries


#cols<-c("#982b2b","#0074b3",
#        "#f47720","#459943","#e5ce81")
##update_geom_defaults("point", list(shape = c(1,2,3,5,6)))
#base <- ggplot(mse.long,aes(x=M,y=mseM,color=method,lty=method,shape=method))+
#  geom_line(aes(x=M,y=mseM),linewidth=0.75,alpha=0.75)+geom_point()
#label_false <- function(labels, multi_line = FALSE, sep ="=") {
#  label_both(labels = labels, multi_line = multi_line, sep = sep)
#}
#Mvaries <- base + facet_grid(n~vaR,scales = "free")+scale_x_log10()+scale_y_log10()+
#  ylab("SPEs")+geom_rug(sides = "b",color="grey",length = unit(0.8, "mm"),linewidth=1)+
#  geom_errorbar(aes(ymin=lower,ymax=upper),lty=1,alpha=0.7)+
#  labs(title = expression(paste("SPEs versus number of sites M for the least squares estimation situation")),
#       subtitle = expression(paste("when q=20 and ",alpha,"=1.5")))+
#  scale_color_manual(values=c("#982b2b","#0074b3","#f47720","#459943","#e5ce81"),name="method",labels=c(expression(DC[opt]),"Full",expression(CSL[1]^{(a)}),expression(DC[ew]),"PD"))+
#  scale_shape_manual(values=c(16,15,17,1,2),name="method",labels=c(expression(DC[opt]),"Full",expression(CSL[1]^{(a)}),expression(DC[ew]),"PD"))+
#  scale_linetype_manual(values =c("solid", "dashed","dotdash","longdash","twodash") ,name="method",labels=c(expression(DC[opt]),"Full",expression(CSL[1]^{(a)}),expression(DC[ew]),"PD"))
#  
#Mvaries




cols<-c("#982b2b","#0074b3",
        "#f47720","#459943","#e5ce81")
#update_geom_defaults("point", list(shape = c(1,2,3,5,6)))
base <- ggplot(mse.long,aes(x=M,y=mseM,color=method,lty=method,shape=method))+
  geom_line(aes(x=M,y=mseM),linewidth=0.75,alpha=0.75)+geom_point()
label_false <- function(labels, multi_line = FALSE, sep ="=") {
  label_both(labels = labels, multi_line = multi_line, sep = sep)
}
#labelname <- c(expression(DC[opt]),"Full",expression(CSL[1]^(a)),expression(DC[ew]),"PD")
labelname <- c(bquote(DC[opt]),bquote(Full[A]),bquote(CSL[1]^(a)),bquote(DC[ew]),bquote(PD))
Mvaries <- base + facet_grid(n~vaR,scales = "free")+scale_x_log10()+scale_y_log10()+
  ylab("SPEs")+geom_rug(sides = "b",color="grey",length = unit(0.8, "mm"),linewidth=1)+
  geom_errorbar(aes(ymin=lower,ymax=upper),lty=1,alpha=0.7)+
  labs(title = expression(paste("SPEs versus number of sites M for the least squares estimation situation")),
       subtitle = expression(paste("when q=30 and ",alpha,"=0.5")))+
  scale_color_manual(values=c("#982b2b","#0074b3","#f47720","#459943","#e5ce81"),name="method",labels=labelname)+
  scale_shape_manual(values=c(16,15,17,1,2),name="method",labels=labelname)+
  scale_linetype_manual(values =c("solid", "dashed","dotdash","longdash","twodash") ,name="method",labels=labelname)

Mvaries









ggsave("ols_p30a05Mvaries.png",Mvaries,device = "png",height = 10,width = 13)
ggsave("ols_p30a05Mvaries.pdf",Mvaries,device = "pdf",height = 10,width = 13)

ggsave("ols_p20a15Mvaries(2).png",Mvaries,device = "png",height = 10,width = 10)
#ggsave("ols_p20a15Mvaries2.png",Mvaries,device = "png",height = 13,width = 13)
#ggsave("p20a15Mvaries.tiff",Mvaries,device = "tiff",height = 10,width = 13)
#ggsave("ols_p15a1Mvaries.pdf",Mvaries,device = "pdf",height = 10,width = 13)
#ggsave("p15a15M10.jpeg",p15a15M10,device = "jpeg",height = 10,width = 10)


#figname <- "p20a15Mvaries.eps"
#savePlot(file=figname,type="eps",dev.cur())
#postscript(file = figname)
#Mvaries
#dev.off()


