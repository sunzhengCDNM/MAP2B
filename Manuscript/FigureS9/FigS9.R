#================================= FigureS9 =================================
library(ggplot2)
library(viridis)
setwd("C:\\Users\\luck\\Desktop\\aDist")
abd_cur <- read.table("CAMI_50_abundance_change_in_ground_truth_nominmax.xls",header = T, sep="\t")
abd_cur <- abd_cur[1:2850,]
#Recall
ggplot(abd_cur,aes(x=log10(Abundance),y=Recall,fill=Datasets))+
  geom_point(aes(color=Datasets))+
  geom_smooth(aes(color=Datasets), size=1,method="loess")+
  theme_bw()+scale_y_continuous(limits = c(0.7,1))+
  scale_fill_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))+
  scale_color_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))
#Precision
ggplot(abd_cur,aes(x=log10(Abundance),y=Precision,fill=Datasets))+
  geom_point(aes(color=Datasets))+
  geom_smooth(aes(color=Datasets), size=1,method="loess")+
  theme_bw()+scale_y_continuous(limits = c(0.7,1))+
  scale_fill_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))+
  scale_color_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))
#F1
ggplot(abd_cur,aes(x=log10(Abundance),y=(Recall*Precision),fill=Datasets))+
  geom_point(aes(color=Datasets))+
  geom_smooth(aes(color=Datasets), size=1,method="loess")+
  theme_bw()+scale_y_continuous(limits = c(0.7,1))+
  scale_fill_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))+
  scale_color_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))
#Accuracy
ggplot(abd_cur,aes(x=log10(Abundance),y=Accuracy,fill=Datasets))+
  geom_point(aes(color=Datasets))+
  geom_smooth(aes(color=Datasets), size=1,method="loess")+
  theme_bw()+scale_y_continuous(limits = c(0.95,1))+
  scale_fill_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))+
  scale_color_manual(values=alpha(c("#440154FF","#B8DE29FF","#FDE725FF")))




#============================AUROC==================================================
#install.packages("PRROC")
library(PRROC)

roc_data <- read.table("CAMI_50_abundance_change_in_ground_truth_100W.txt",header = T, sep="\t")
roc_data <- read.table("CAMI_50_abundance_change_in_ground_truth_10W.txt",header = T, sep="\t")
roc_data <- read.table("CAMI_50_abundance_change_in_ground_truth_1W.txt",header = T, sep="\t")

fg <- roc_data[which(roc_data$Class=="strmgCAMI2"),][roc_data[which(roc_data$Class=="strmgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="strmgCAMI2"),][roc_data[which(roc_data$Class=="strmgCAMI2"),]$Tag == 0,]
roc <- roc.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(roc,main="ROC (all ground truth)",col="#FDE725FF")
#text(0.5,0.1,paste("AUROC=",round(roc$auc,3)),col="#FDE725FF")

fg <- roc_data[which(roc_data$Class=="rhimgCAMI2"),][roc_data[which(roc_data$Class=="rhimgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="rhimgCAMI2"),][roc_data[which(roc_data$Class=="rhimgCAMI2"),]$Tag == 0,]
roc <- roc.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(roc,col="#B8DE29FF",add=T)
#text(0.5,0.2,paste("AUROC=",round(roc$auc,3)),col="#B8DE29FF")

fg <- roc_data[which(roc_data$Class=="marmgCAMI2"),][roc_data[which(roc_data$Class=="marmgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="marmgCAMI2"),][roc_data[which(roc_data$Class=="marmgCAMI2"),]$Tag == 0,]
roc <- roc.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(roc,col="#33638DFF",add=T)
#text(0.5,0.3,paste("AUROC=",round(roc$auc,3)),col="#440154FF")#33638DFF

#============================AUPRC==================================================
fg <- roc_data[which(roc_data$Class=="strmgCAMI2"),][roc_data[which(roc_data$Class=="strmgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="strmgCAMI2"),][roc_data[which(roc_data$Class=="strmgCAMI2"),]$Tag == 0,]
pr <- pr.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(pr,main="PRC (all ground truth)",col="#FDE725FF")
#text(0.5,0.1,paste("AUPRC=",round(pr$auc.integral,3)),col="#FDE725FF")

fg <- roc_data[which(roc_data$Class=="rhimgCAMI2"),][roc_data[which(roc_data$Class=="rhimgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="rhimgCAMI2"),][roc_data[which(roc_data$Class=="rhimgCAMI2"),]$Tag == 0,]
pr <- pr.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(pr,col="#B8DE29FF",add=T)
#text(0.5,0.2,paste("AUPRC=",round(pr$auc.integral,3)),col="#B8DE29FF")

fg <- roc_data[which(roc_data$Class=="marmgCAMI2"),][roc_data[which(roc_data$Class=="marmgCAMI2"),]$Tag == 1,]
bg <- roc_data[which(roc_data$Class=="marmgCAMI2"),][roc_data[which(roc_data$Class=="marmgCAMI2"),]$Tag == 0,]
pr <- pr.curve(scores.class0 = fg$Probability, scores.class1 = bg$Probability, curve = T)
#plot(pr,col="#440154FF",add=T)
#text(0.5,0.3,paste("AUPRC=",round(pr$auc.integral,3)),col="#440154FF")##33638DFF