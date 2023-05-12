#################################################################
# Function: Multivariate statistical analysis based on distance matrix
# Call: Rscript PM_Bdiversity.R -m map_file -d dist_file -o output
# R packages used: reshape,ggplot2,pheatmap,pROC,combinat,plyr,vegan,optparse
# Last update: 2017-03-01, Shi Huang, Xiaoquan Su, Gongchao jing
#################################################################
options(warn=-1)
#Rprof()
## install necessary libraries
p <- c("reshape","ggplot2","pheatmap","pROC","combinat","plyr","vegan","optparse")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
sourcedir <- Sys.getenv("ParallelMETA")
source(sprintf('%s/Rscript/util.R',sourcedir))
# make option list and parse command line
option_list <- list(
    make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix table [Required]"),
    make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
    make_option(c("-o", "--out_dir"), type="character", default='Beta_diversity', help="Output directory [default %default]"),
    make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
    make_option(c("-n", "--dist_name"), type="character", default='Default', help="The distance metrics name such as Meta-Storms, Jensen-Shannon, Euclidean et al. [Optional, default %default]")    
    )
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$meta_data)) stop('Please input a meta data file')
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

# create output directory if needed
#if(opts$out_dir != ".") 
dir.create(outpath1<-paste(opts$out_dir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

filename<-opts$dist_file                       
metadata.filename<-opts$meta_data               
dm_name<-opts$dist_name
prefix<-opts$prefix                         

con <- file(paste(opts$out_dir,'/',prefix,'.Beta_diversity_log.txt',sep=''))
sink(con, append=TRUE)
sink(con, append=TRUE, type='message')

#--------------------------------
dm<-read.table(filename,header=T,row.names=1)
dm<-dm[order(rownames(dm)),order(colnames(dm))]
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
if(length(allmetadata)==1){metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
                           all_group<-colnames(metadata)<-colnames(allmetadata)
                           }else{
                           metadata<-allmetadata[order(rownames(allmetadata)),which(sapply(allmetadata,class)!=0)]
                           all_group<-colnames(metadata)
                           all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]
                           all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]
                           }
cat("All the sample metadata: ",all_group, "\n\n",sep=" ")
#--------------------------------Data Check
if(any((colnames(dm)==rownames(dm))==FALSE)) 
  {cat("The column names do not exactly match the row names! Please revise!")}
if(any((rownames(metadata)==rownames(dm))==FALSE)) 
  {cat("The row names in Map file do not exactly match the row names in the distance matrix! Please revise!\n")}

#--------------------------------
# Statistical test: Adonis and Anosim
#--------------------------------
stat_summ<-matrix(NA,nrow=length(all_group),ncol=4)
rownames(stat_summ)<-all_group;colnames(stat_summ)<-c("Adonis.F","Adonis.P","Anosim.R","Anosim.P")
#--------------------------------
#suppressWarnings(
for(group in all_group) {
    #--------------------------------
    if(is.element(group,all_group_f)){
    ano<-anosim(dm,metadata[,group])
    stat_summ[group,4]<-ano.P<-ano$signif
    stat_summ[group,3]<-ano.R<-ano$statistic
    cat("ANOSIM (",group,"): \n")
    cat("--------------------------------")
    print(ano)
    }
    #--------------------------------
    ado<-adonis(dm~metadata[,group])
    stat_summ[group,2]<-ado.P<-ado$aov.tab$P[1]
    stat_summ[group,1]<-ado.F<-ado$aov.tab$F.Model[1]
    cat("ADONIS/PERMANOVA (",group,"): \n")
    cat("--------------------------------\n")
    print(ado$aov.tab)
    cat("--------------------------------\n\n")

}
#)
sink(paste(outpath1,prefix,".Beta_diversity_summ.xls",sep=""));cat("\t");write.table(stat_summ,quote=FALSE,sep='\t',row.names=TRUE);sink()

#--------------------------------
# Categorical metadata: Distance boxplot
#--------------------------------
if(length(all_group_f)>=1){

for(group in all_group_f) {
    #--------------------------------
    dir.create(outpath2<-paste(outpath1,prefix,".",group,"/",sep=""))
    d<-DistBoxplot(dm,dm_name=dm_name,group=metadata[,group],group_name=group,outpath=outpath2)
    #unlink(outpath2, recursive=TRUE)
    #--------------------------------
    # print(paste(" All_between ",group," VS All_within ",group," P value (T-test)=", d$p_t,sep=""))
    # print(paste(" All_between ",group," VS All_within ",group," P value (Wilcox-test)=", d$p_w,sep=""))
    #--------------------------------
    dm_value<-data.frame(Grouping=rep(group,nrow(d$dm_value)),d$dm_value)
    if(group==all_group_f[1]){
                       dm_v<-dm_value}else{
                       dm_v<-rbind(dm_v,dm_value)}   
}
#--------------------------------
p<-qplot(x=Grouping, y=Dist, data=dm_v, geom="boxplot", fill=DistType, position="dodge",main="", ylab=paste(dm_name," Distance",sep=""))+ coord_flip()+ theme_bw()
suppressMessages(ggsave(filename=paste(outpath1,"/",prefix,".Beta_diversity_DistBoxplot.pdf",sep=""),plot=p, limitsize=TRUE, width=5, height=ifelse(length(all_group_f)>1,length(all_group_f),2)))
}

#--------------------------------
# Continuous metadata: Scatterplot
#--------------------------------
if(length(all_group_n)>=1){
for(group in all_group_n) {
    #--------------------------------
    dir.create(outpath2<-paste(outpath1,prefix,".",group,"/",sep=""))
    dm_n<-data.matrix(dist(metadata[,group]))
    dm_n_value<-dm_n[lower.tri(dm_n, diag = FALSE)]
    dm_value<-dm[lower.tri(dm, diag = FALSE)]
    dm_data<-data.frame(dm=dm_value,group=dm_n_value)
    corr<-with(dm_data,cor.test(dm_value,dm_n_value,method="spearman"))
    p<-ggplot(data=dm_data, aes(x=group,y=dm))+geom_point(alpha=0.2)+ 
       annotate("text", x=max(dm_data$group)*0.9, y=max(dm_data$dm)*0.9, label= paste("Rho=",round(corr$estimate,2),"\n","P=",round(corr$p.value,4),"\n",sep="")) +
       ylab(paste(dm_name," Distance",sep=""))+ 
       xlab(bquote(paste(~Delta, .(group),sep=" ")))+
       theme_bw() 
    if(corr$estimate>0.4 && corr$p.value<0.01) p<-p+geom_smooth(method = "loess", se=TRUE, span=1) 
suppressMessages(ggsave(filename=paste(outpath2,"/",prefix,".Beta_diversity_",group,".Scatterplot",".pdf",sep=""),plot=p, limitsize=TRUE, width=4, height=4))
}
}


#Rprof(NULL)
#summaryRprof()

sink()
sink(type='message')
invisible(system("rm -rf Rplots.pdf"))

