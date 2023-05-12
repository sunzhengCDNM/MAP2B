library(ggplot2)
  library(ggforce)
  library(umap)
  library(ape)
  library(concaveman)
  
  metadata <- read.table("metadata.txt",header=T,row.names = 1,sep="\t")
  distmap2b <- read.table("dist_2b_test.dist",header=T,row.names = 1,sep="\t")
  distmap2b_cov <- read.table("dist_2bcov_test.dist",header=T,row.names = 1,sep="\t")
  distmotu <- read.table("dist_mOTUs_test.dist",header=T,row.names = 1,sep="\t")
  distbraken <- read.table("dist_Bracken_test.dist",header=T,row.names = 1,sep="\t")
  distkracken <- read.table("dist_kraken_test.dist",header=T,row.names = 1,sep="\t")
  distKrackenuniq <- read.table("dist_KrakenUniq_test.dist",header=T,row.names = 1,sep="\t")
  distmpa3 <- read.table("dist_mpa3_test.dist",header=T,row.names = 1,sep="\t")
  
  map2b_pcoa <- pcoa(distmap2b)
  map2b_pcoa_gg <- Mpa2b_pcoa$vectors[,1:2]
  ggplot(data.frame(map2b_pcoa_gg),aes(x=map2b_pcoa_gg[,1],y=map2b_pcoa_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_map2b_gg_pcoa.pdf", width = 9, height = 4)
  
  Mpa2b_cov_pcoa <- pcoa(distmap2b_cov)
  Map2b_cov_gg <- Mpa2b_cov_pcoa$vectors[,1:2]
  ggplot(data.frame(Map2b_cov_gg),aes(x=Map2b_cov_gg[,1],y=Map2b_cov_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_Map2b_cov_gg.pdf", width = 9, height = 4)
  
  motus_pcoa <- pcoa(distmotu)
  motus_gg <- motus_pcoa$vectors[,1:2]
  ggplot(data.frame(motus_gg),aes(x=motus_gg[,1],y=motus_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_motus_gg.pdf", width = 9, height = 4)
  
  bracken_pcoa <- pcoa(distbraken)
  bracekn_gg <- bracken_pcoa$vectors[,1:2]
  ggplot(data.frame(bracekn_gg),aes(x=bracekn_gg[,1],y=bracekn_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_bracekn_gg_BC.pdf", width = 9, height = 4)
  
  kraken_pcoa <- pcoa(distkracken)
  kraken_gg <- kraken_pcoa$vectors[,1:2]
  ggplot(data.frame(kraken_gg),aes(x=kraken_gg[,1],y=kraken_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_kraken_gg_BC.pdf", width = 9, height = 4)
  
  krakenU_pcoa <- pcoa(distKrackenuniq)
  KrakenU_gg <- krakenU_pcoa$vectors[,1:2]
  ggplot(data.frame(KrakenU_gg),aes(x=KrakenU_gg[,1],y=KrakenU_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_Krackenunique_gg_BC.pdf", width = 9, height = 4)
  
  mpa3_pcoa <- pcoa(distmpa3)
  mpa3_gg <- mpa3_pcoa$vectors[,1:2]
  ggplot(data.frame(mpa3_gg),aes(x=mpa3_gg[,1],y=mpa3_gg[,2],group=metadata$diagnosis))+geom_point(aes(color=metadata$diagnosis))+theme_bw()+facet_wrap(~metadata$batch,scales = "free")+ stat_ellipse(level = 0.95)
  ggsave("IBD_mpa3_gg_BC.pdf", width = 9, height = 4)
  
  