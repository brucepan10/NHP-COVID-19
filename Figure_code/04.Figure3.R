library(viridis)
library(ggplot2)
library(Seurat)
library(ggpubr)

########################################################## 1. Figure 3A #####################################################################
## Lung
load("Lung.integrated.Find.RData")
Lung_monkey=Lung.integrated.Find
FeaturePlot(Lung_monkey,features = c("ACE2","TMPRSS2"),order = T,cols = viridis(50))

########################################################## 2. Figure 3B #####################################################################
## kidney
load("kidney.integrated.Find.RData")
Kidney_monkey=kidney.integrated.Find
FeaturePlot(Kidney_monkey,features = c("ACE2","TMPRSS2"),order = T,cols = viridis(50))

########################################################## 3. Figure 3C #####################################################################
## Liver
load("Liver.integrated.Find.RData")
Liver_monkey=Liver.integrated.Find
FeaturePlot(Liver_monkey,features = c("ACE2","TMPRSS2"),order = T,cols = viridis(50))

########################################################## 4. Figure 3D #####################################################################
lung_human=readRDS("lung_ts.rds")
load("ALL_cells_Lung.integrated.Find.RData")
cell=sort(unique(c(unique(lung_human$celltype),unique(as.character(lung_monkey$celltype)))))
lung_human$celltype=factor(lung_human$celltype,levels = cell)
lung_monkey$celltype=factor(lung_monkey$celltype,levels = cell)
DefaultAssay(lung_monkey)="RNA"
RNA_Human_sub=as.data.frame(t(as.data.frame(lung_human@assays$RNA@data[c("ACE2","TMPRSS2"),])))
RNA_Monkey_sub=as.data.frame(t(as.data.frame(lung_monkey@assays$RNA@data[c("ACE2","TMPRSS2"),])))
RNA_Human_sub$celltype=lung_human$celltype
RNA_Monkey_sub$celltype=lung_monkey$celltype
ss_Human_ACE2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$ACE2>0)[,3]))
ss_Human_TMPRSS2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$TMPRSS2>0)[,3]))
tab_type_Human=as.data.frame(table(RNA_Human_sub$celltype))
sum_Human=cbind(ss_Human_ACE2[2],ss_Human_TMPRSS2[2],tab_type_Human)
sum_Human$Species="Human"
sum_Human$exp_ACE2=unlist(tapply(RNA_Human_sub$ACE2,RNA_Human_sub$celltype,mean))
sum_Human$exp_TMPRSS2=unlist(tapply(RNA_Human_sub$TMPRSS2,RNA_Human_sub$celltype,mean))
ss_Monkey_ACE2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$ACE2>0)[,3]))
ss_Monkey_TMPRSS2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$TMPRSS2>0)[,3]))
tab_type_Monkey=as.data.frame(table(RNA_Monkey_sub$celltype))
sum_Monkey=cbind(ss_Monkey_ACE2[2],ss_Monkey_TMPRSS2[2],tab_type_Monkey)
sum_Monkey$Species="Monkey"
sum_Monkey$exp_ACE2=unlist(tapply(RNA_Monkey_sub$ACE2,RNA_Monkey_sub$celltype,mean))
sum_Monkey$exp_TMPRSS2=unlist(tapply(RNA_Monkey_sub$TMPRSS2,RNA_Monkey_sub$celltype,mean))
all=rbind(sum_Human,sum_Monkey)
all$Ratio_ACE2=all[,1]/all[,4]
all$Ratio_TMPRSS2=all[,2]/all[,4]
all[is.na(all)]=0
all1=subset(all,all[,4]>0)
tab=as.data.frame(table(all1$Var1))
tab_sub=subset(tab,tab$Freq==2)
all2=merge(all1,tab_sub[1],by="Var1")
###Bubble plot ACE2
all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order=all2_Monkey[order(all2_Monkey$Ratio_ACE2,decreasing = T),]
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = all2_Monkey_order$Var1)
b1=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = levels(all2_Monkey_order$Var1))
b2=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
###Bubble plot TMPRSS2
all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = levels(all2_Monkey_order$Var1))
b3=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = levels(all2_Monkey_order$Var1))
b4=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
pdf("Figure 3D.pdf")
ggarrange(b1,b3,b2,b4,ncol=2,nrow=2)
dev.off()

########################################################## 5. Figure 3E #####################################################################
## get Bubble plot
## input cell type information and ACE2 TMPRSS2 expression matrix of human and monkey kidney
RNA_Human_sub=read.table("Kidney_human_ACE2_TMPRSS2_celltype_info.txt",header = T,sep = "\t")
RNA_Monkey_sub=read.table("Kidney_monkey_ACE2_TMPRSS2_celltype_info.txt",header = T,sep = "\t")

ss_Human_ACE2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$ACE2>0)[,3]))
ss_Human_TMPRSS2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$TMPRSS2>0)[,3]))
tab_type_Human=as.data.frame(table(RNA_Human_sub$celltype))
sum_Human=cbind(ss_Human_ACE2[2],ss_Human_TMPRSS2[2],tab_type_Human)
sum_Human$Species="Human"
sum_Human$exp_ACE2=unlist(tapply(RNA_Human_sub$ACE2,RNA_Human_sub$celltype,mean))
sum_Human$exp_TMPRSS2=unlist(tapply(RNA_Human_sub$TMPRSS2,RNA_Human_sub$celltype,mean))

ss_Monkey_ACE2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$ACE2>0)[,3]))
ss_Monkey_TMPRSS2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$TMPRSS2>0)[,3]))
tab_type_Monkey=as.data.frame(table(RNA_Monkey_sub$celltype))
sum_Monkey=cbind(ss_Monkey_ACE2[2],ss_Monkey_TMPRSS2[2],tab_type_Monkey)
sum_Monkey$Species="Monkey"
sum_Monkey$exp_ACE2=unlist(tapply(RNA_Monkey_sub$ACE2,RNA_Monkey_sub$celltype,mean))
sum_Monkey$exp_TMPRSS2=unlist(tapply(RNA_Monkey_sub$TMPRSS2,RNA_Monkey_sub$celltype,mean))

all=rbind(sum_Human,sum_Monkey)

all$Ratio_ACE2=all[,1]/all[,4]
all$Ratio_TMPRSS2=all[,2]/all[,4]
all[is.na(all)]=0

all1=subset(all,all[,4]>0)
tab=as.data.frame(table(all1$Var1))
tab_sub=subset(tab,tab$Freq==2)
all2=merge(all1,tab_sub[1],by="Var1")

#####ACE2

all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order=all2_Monkey[order(all2_Monkey$Ratio_ACE2,decreasing = T),]
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = all2_Monkey_order$Var1)
b1=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = levels(all2_Monkey_order$Var1))
b2=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

#####TMPRSS2

all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order=all2_Monkey[order(all2_Monkey$Ratio_TMPRSS2,decreasing = T),]
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = all2_Monkey_order$Var1)
b3=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = all2_Monkey_order$Var1)
b4=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

pdf("Figure 3E.pdf")
ggarrange(b1,b3,b2,b4,ncol=2,nrow=2)
dev.off()
########################################################## 6. Figure 3F #####################################################################
load("ALL_cells_LiverL.integrated.Find.RData")
Liver_monkey=LiverL.integrated.Find
load("Liver_human.RData")
cell=sort(unique(c(unique(Liver_human$celltype),unique(Liver_monkey$celltype))))
Liver_human$celltype=factor(Liver_human$celltype,levels = cell)
Liver_monkey$celltype=factor(Liver_monkey$celltype,levels = cell)
DefaultAssay(Liver_monkey)="RNA"
RNA_Human_sub=as.data.frame(t(as.data.frame(Liver_human@assays$RNA@data[c("ACE2","TMPRSS2"),])))
RNA_Monkey_sub=as.data.frame(t(as.data.frame(Liver_monkey@assays$RNA@data[c("ACE2","TMPRSS2"),])))
RNA_Human_sub$celltype=Liver_human$celltype
RNA_Monkey_sub$celltype=Liver_monkey$celltype
ss_Human_ACE2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$ACE2>0)[,3]))
ss_Human_TMPRSS2=as.data.frame(table(subset(RNA_Human_sub,RNA_Human_sub$TMPRSS2>0)[,3]))
tab_type_Human=as.data.frame(table(RNA_Human_sub$celltype))
sum_Human=cbind(ss_Human_ACE2[2],ss_Human_TMPRSS2[2],tab_type_Human)
sum_Human$Species="Human"
sum_Human$exp_ACE2=unlist(tapply(RNA_Human_sub$ACE2,RNA_Human_sub$celltype,mean))
sum_Human$exp_TMPRSS2=unlist(tapply(RNA_Human_sub$TMPRSS2,RNA_Human_sub$celltype,mean))
ss_Monkey_ACE2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$ACE2>0)[,3]))
ss_Monkey_TMPRSS2=as.data.frame(table(subset(RNA_Monkey_sub,RNA_Monkey_sub$TMPRSS2>0)[,3]))
tab_type_Monkey=as.data.frame(table(RNA_Monkey_sub$celltype))
sum_Monkey=cbind(ss_Monkey_ACE2[2],ss_Monkey_TMPRSS2[2],tab_type_Monkey)
sum_Monkey$Species="Monkey"
sum_Monkey$exp_ACE2=unlist(tapply(RNA_Monkey_sub$ACE2,RNA_Monkey_sub$celltype,mean))
sum_Monkey$exp_TMPRSS2=unlist(tapply(RNA_Monkey_sub$TMPRSS2,RNA_Monkey_sub$celltype,mean))
all=rbind(sum_Human,sum_Monkey)
all$Ratio_ACE2=all[,1]/all[,4]
all$Ratio_TMPRSS2=all[,2]/all[,4]
all[is.na(all)]=0
all1=subset(all,all[,4]>0)
tab=as.data.frame(table(all1$Var1))
tab_sub=subset(tab,tab$Freq==2)
all2=merge(all1,tab_sub[1],by="Var1")
###Bubble plot ACE2
all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order=all2_Monkey[order(all2_Monkey$Ratio_ACE2,decreasing = T),]
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = all2_Monkey_order$Var1)
b1=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = levels(all2_Monkey_order$Var1))
b2=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_ACE2))+
  geom_point() +
  aes_string(color="exp_ACE2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
###Bubble plot TMPRSS2
all2_Monkey=subset(all2,all2$Species=="Monkey")
all2_Monkey_order=all2_Monkey[order(all2_Monkey$Ratio_TMPRSS2,decreasing = T),]
all2_Monkey_order$Var1=factor(all2_Monkey_order$Var1,levels = all2_Monkey_order$Var1)
b3=ggplot(all2_Monkey_order,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))

all2_Human=subset(all2,all2$Species=="Human")
all2_Human$Var1=factor(all2_Human$Var1,levels = all2_Monkey_order$Var1)
b4=ggplot(all2_Human,aes_(x = ~Var1, y = ~Species, size = ~Ratio_TMPRSS2))+
  geom_point() +
  aes_string(color="exp_TMPRSS2") +
  theme_bw()+
  theme(axis.text.x = element_text(size=8,angle=90, hjust=1, vjust=1))+
  scale_color_gradientn(colours = viridis(50,direction = -1,option = "C"))
pdf("Figure 3F.pdf")
ggarrange(b1,b3,b2,b4,ncol=2,nrow=2)
dev.off()
