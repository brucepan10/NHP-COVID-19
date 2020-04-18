########################################################## Figure S3 #####################################################################
library(Seurat)
library(ggplot2)
library(viridis)

load("Pancreas.integrated.Find.RData") ## input RData of Pancreas
pdf("FigureS3_Monkey_Pancreas_ACE2.pdf",5.2,5)
p1<-FeaturePlot(Pancreas.integrated.Find,features = "ACE2",order = T)+scale_color_viridis()
p1
dev.off()
pdf("FigureS3_Monkey_Pancreas_TMPRSS2.pdf",5.2,5)
p2<-FeaturePlot(Pancreas.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p2
dev.off()

load("Parotid.integrated.Find.RData") ## input RData of Parotid
pdf("FigureS3_Monkey_Parotid_ACE2.pdf",5.2,5)
p3<-FeaturePlot(Parotid.integrated.Find,features = "ACE2",order = T)+scale_color_viridis()
p3
dev.off()
pdf("FigureS3_Monkey_Parotid_TMPRSS2.pdf",5.2,5)
p4<-FeaturePlot(Parotid.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p4
dev.off()

load("PBMC.integrated.Find.RData") ## input RData of PBMC
pdf("FigureS3_Monkey_PBMC_ACE2.pdf",5.2,5)
p5<-FeaturePlot(PBMC.integrated.Find,features = "ACE2",order = T)+scale_color_viridis()
p5
dev.off()
pdf("FigureS3_Monkey_PBMC_TMPRSS2.pdf",5.2,5)
p6<-FeaturePlot(PBMC.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p6
dev.off()

load("Thyroid.integrated.Find.RData") ## input RData of Thyroid
pdf("FigureS3_Monkey_Thyroid_ACE2.pdf",5.2,5)
p7<-FeaturePlot(Thyroid.integrated.Find,features = "ACE2",order = T)+scale_color_viridis()
p7
dev.off()
pdf("FigureS3_Monkey_Thyroid_TMPRSS2.pdf",5.2,5)
p8<-FeaturePlot(Thyroid.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p8
dev.off()

load("Aorta.integrated.Find.RData") ## input RData of Aorta
pdf("FigureS3_Monkey_Artery_ACE2.pdf",5.2,5)
p9<-FeaturePlot(Artery.integrated.Find,features = "ACE2",order = T)+scale_color_viridis()
p9
dev.off()
pdf("FigureS3_Monkey_Artery_TMPRSS2.pdf",5.2,5)
p10<-FeaturePlot(Artery.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p10
dev.off()

load("Neocortex.integrated.Find.RData") ## input RData of Neocortex
dim(Neocortex.integrated.Find@assays$RNA@data)
cluster_cor_umap = as.data.frame(Embeddings(object = Neocortex.integrated.Find,reduction = "umap"))
cluster_cor_umap$ACE2<-0

pdf("FigureS3_Monkey_Neocortex_ACE2.pdf",5.2,5)
	p11 <- ggplot(cluster_cor_umap,aes(x=cluster_cor_umap$UMAP_1,y=cluster_cor_umap$UMAP_2,colour=ACE2))+geom_point(size = 1)+
	theme_classic()+
	theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA))+
	ylab ('UMAP1')+xlab('UMAP2')+
	scale_colour_gradientn(colours = c("#440153","#440153","#440153"))+
	ggtitle("ACE2")
dev.off()
pdf("FigureS3_Monkey_Neocortex_TMPRSS2.pdf",5.2,5)
p12<-FeaturePlot(Neocortex.integrated.Find,features = "TMPRSS2",order = T)+scale_color_viridis()
p12
dev.off()

require(gridExtra)
pdf("FigureS3_Monkey_ACE2_TMPRSS2.pdf",12,32)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol=2)
dev.off()



