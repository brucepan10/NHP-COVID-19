# input cluster, umap coordinates and other information of all cells
scanpy <- read.csv(file="Scanpy_cell_info.csv", header=TRUE, sep=",")
head(scanpy)
scanpy_umap <- read.csv(file="Scanpy_X_pca_umap.csv", header=TRUE, sep=",")
head(scanpy_umap)
scanpy$X_umap1 <- scanpy_umap$X_umap1
scanpy$X_umap2 <- scanpy_umap$X_umap2

colnames(scanpy)[1]="ID"
scanpy$organ<-unlist(lapply(X = as.character(scanpy$ID), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))

library(RColorBrewer)
library(ggplot2)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
a=getPalette(9)

########################################################## 1. Figure 1B #####################################################################
## get umap_organ graph
pdf("umap_organ9.pdf",width=11,height=10)
ggplot(data=scanpy,aes(x=X_umap1, y =X_umap2))+
  geom_point(alpha=1, size=0.5,aes(colour=organ))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_color_manual(values = a)+
  theme_classic()+
  labs(x='UMAP1',y='UMAP2',title = "organ")+
  theme(plot.title = element_text(size=rel(2),hjust = 0.5),legend.position = 'right')+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  theme(legend.title = element_text(colour="black", size=16))+
  theme(legend.text = element_text(colour="black", size = 10))
dev.off()

## get umap_organ graph legend
sub_scanpy<-scanpy[!duplicated(scanpy$organ),]
pdf("umap_organ9_legend.pdf",width=11,height=10)
ggplot(data=sub_scanpy,aes(x=X_umap1, y =X_umap2))+
  geom_point(alpha=1, size=0.5,aes(colour=organ))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_color_manual(values = a)+
  theme_classic()+
  labs(x='UMAP1',y='UMAP2',title = "organ")+
  theme(plot.title = element_text(size=rel(2),hjust = 0.5),legend.position = 'right')+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  theme(legend.title = element_text(colour="black", size=16))+
  theme(legend.text = element_text(colour="black", size = 10))
dev.off()

## get organ cellnum bar graph
organ_cellnum<-as.data.frame(table(scanpy$organ))
organ_cellnum$Var1=factor(organ_cellnum$Var1,levels = rev(as.character(organ_cellnum$Var1)))
a_rev=rev(a)

ggplot(organ_cellnum,aes(Var1,log10(Freq),fill=Var1))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values = a_rev)+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

########################################################## 2. Figure 1C #####################################################################
### celltype annotation
## split cluster 67 to distinguish between Oligodendrocyte and Oligodendrocyte precursor cell
head(scanpy)
scanpy$Batch<-"Others"
pos1<-which(scanpy$louvain==67)
scanpy$Batch[pos1]<-67
pos2<-which(scanpy$Batch=="Others")
ggplot()+
  geom_point(data=scanpy[pos2,],aes(x=X_umap1, y =X_umap2),alpha=1, size=0.001,colour="gainsboro")+
  geom_point(data=scanpy[pos1,],aes(x=X_umap1, y =X_umap2),alpha=1, size=0.001,colour="red")+
  labs(x='',y='',title = 67)+
  theme(plot.title = element_text(size=rel(2),hjust = 0.5))+
  theme(legend.position="none")

scanpy$celltype<-0
scanpy$celltype[which(scanpy$louvain==0|scanpy$louvain==1|scanpy$louvain==2)]<-"Parotid-acinar cell"
scanpy$celltype[which(scanpy$louvain==3|scanpy$louvain==28)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==4)]<-"Follicular cell"
scanpy$celltype[which(scanpy$louvain==5|scanpy$louvain==10|scanpy$louvain==13|scanpy$louvain==16|scanpy$louvain==17|scanpy$louvain==27|scanpy$louvain==29|scanpy$louvain==32|scanpy$louvain==37|scanpy$louvain==61)]<-"Hepatocyte"
scanpy$celltype[which(scanpy$louvain==6)]<-"CD4+ T cell"
scanpy$celltype[which(scanpy$louvain==7|scanpy$louvain==30)]<-"Proximal tubule cell 1"
scanpy$celltype[which(scanpy$louvain==8)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==9)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==11)]<-"B cell"
scanpy$celltype[which(scanpy$louvain==12)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==14|scanpy$louvain==22|scanpy$louvain==31)]<-"Thick ascending limb cell"
scanpy$celltype[which(scanpy$louvain==15)]<-"Principal cell 1"
scanpy$celltype[which(scanpy$louvain==18)]<-"CD4+ T cell"
scanpy$celltype[which(scanpy$louvain==19)]<-"Principal cell 2"
scanpy$celltype[which(scanpy$louvain==20)]<-"Ciliated cell"
scanpy$celltype[which(scanpy$louvain==21)]<-"CD8+ T cell"
scanpy$celltype[which(scanpy$louvain==23)]<-"Connecting/Distal convoluted tubule cell"
scanpy$celltype[which(scanpy$louvain==24)]<-"Parotid-acinar cell"
scanpy$celltype[which(scanpy$louvain==25)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==26|scanpy$louvain==57)]<-"Excitatory neuron"
scanpy$celltype[which(scanpy$louvain==33)]<-"Pulmonary alveolar type 2 cell"
scanpy$celltype[which(scanpy$louvain==34)]<-"Liver-endothelial cell"
scanpy$celltype[which(scanpy$louvain==35)]<-"Endothelial cell"
scanpy$celltype[which(scanpy$louvain==36)]<-"Proximal tubule cell 2"
scanpy$celltype[which(scanpy$louvain==38)]<-"Thin limb cell"
scanpy$celltype[which(scanpy$louvain==39)]<-"Striated duct cell"
scanpy$celltype[which(scanpy$louvain==40)]<-"Proximal tubule cell 2"
scanpy$celltype[which(scanpy$louvain==41)]<-"Lung-macrophage"
scanpy$celltype[which(scanpy$louvain==42)]<-"Pulmonary alveolar type 1 cell"
scanpy$celltype[which(scanpy$louvain==43)]<-"Vascular smooth muscle cell"
scanpy$celltype[which(scanpy$louvain==44)]<-"Mucous acinar cell"
scanpy$celltype[which(scanpy$louvain==45|scanpy$louvain==73|scanpy$louvain==75)]<-"Parotid-acinar cell"
scanpy$celltype[which(scanpy$louvain==46)]<-"Pancreas-islet cell"
scanpy$celltype[which(scanpy$louvain==47|scanpy$louvain==60)]<-"Stromal cell"
scanpy$celltype[which(scanpy$louvain==48|scanpy$louvain==50)]<-"Intercalated cell 1"
scanpy$celltype[which(scanpy$louvain==49)]<-"Intercalated cell 2"
scanpy$celltype[which(scanpy$louvain==51)]<-"Macrophage"
scanpy$celltype[which(scanpy$louvain==52)]<-"Podocyte"
scanpy$celltype[which(scanpy$louvain==53)]<-"Pancreas-duct cell"
scanpy$celltype[which(scanpy$louvain==54)]<-"Pancreas-stromal cell"
scanpy$celltype[which(scanpy$louvain==55)]<-"Lung-macrophage"  
scanpy$celltype[which(scanpy$louvain==56)]<-"Inhibitory neuron"
scanpy$celltype[which(scanpy$louvain==58)]<-"Club cell"
scanpy$celltype[which(scanpy$louvain==59)]<-"Lung-smooth muscle cell"
scanpy$celltype[which(scanpy$louvain==62)]<-"Astrocyte"
scanpy$celltype[which(scanpy$louvain==63)]<-"B cell"
scanpy$celltype[which(scanpy$louvain==64)]<-"Adipocyte"
scanpy$celltype[which(scanpy$louvain==65)]<-"Myoepithelial cell"
scanpy$celltype[which(scanpy$louvain==66|scanpy$louvain==69|scanpy$louvain==70)]<-"Pancreas-acinar cell"
scanpy$celltype[which(scanpy$louvain==68)]<-"Lung-stromal cell"
scanpy$celltype[which(scanpy$louvain==71)]<-"Cholangiocyte"
scanpy$celltype[which(scanpy$louvain==72)]<-"Thyroid-endothelial cell"
scanpy$celltype[which(scanpy$louvain==74)]<-"Thyroid-smooth muscle cell"
scanpy$celltype[which(scanpy$X_umap2>0&scanpy$louvain==67)]<-"Oligodendrocyte"
scanpy$celltype[which(scanpy$X_umap2<0&scanpy$louvain==67)]<-"Oligodendrocyte precursor cell"

### color palette
best_color<- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
               "#006FA6","#A30059","#FFE4E1","#0000A6","#B79762",
               "#004D43","#8FB0FF","#997D87","#5A0007","#809693",
               "#1B4400","#4FC601","#3B5DFF","#BA0900","#FF2F80",
               "#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA",
               "#DDEFFF","#7B4F4B","#A1C299","#0AA6D8","#00A087FF",
               "#4DBBD5FF","#E64B35FF","#3C5488FF","#0067A5","#63FFAC",
               "#F38400","#A1CAF1", "#C2B280","#848482","#E68FAC",  
               "#F99379", "#604E97","#F6A600", "#B3446C","#DCD300",
               "#882D17", "#8DB600","#654522", "#E25822", "#2B3D26",
               "#191970","#000080",
               "#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493",
               "#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F",
               "#00FF00","#76EE00","#FFF68F","Yellow1","Gold1",
               "DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24",
               "#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347",
               "#FF4500","#FF1493","#FF6EB4","#EE30A7","#8B008B")

## get umap_celltype graph
pdf("umap_celltype_best_color.pdf",width=15,height=10)
p1<-ggplot(data=scanpy,aes(x=X_umap1, y =X_umap2))+
  geom_point(alpha=1, size=0.5,aes(colour=celltype))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values = best_color)+
  theme_classic()+
  labs(x='UMAP1',y='UMAP2',title = "celltype")+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 8))
p1
dev.off()

## get umap_celltype graph legend
sub_scanpy<-scanpy[!duplicated(scanpy$celltype),]
pdf("umap_celltype_best_color_legend.pdf",width=15,height=10)
p2<-ggplot(data=sub_scanpy,aes(x=X_umap1, y =X_umap2))+
  geom_point(alpha=1, size=0.5,aes(colour=celltype))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values = best_color)+
  theme_classic()+
  labs(x='UMAP1',y='UMAP2',title = "celltype_legend")+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 8))
p2
dev.off()

## get celltype cellnum bar graph 
celltype_cellnum<-as.data.frame(table(scanpy$celltype))
celltype_cellnum$Var1=factor(celltype_cellnum$Var1,levels = rev(as.character(celltype_cellnum$Var1)))
b_rev=rev(best_color[1:44])

ggplot(celltype_cellnum,aes(Var1,log10(Freq),fill=Var1))+
  geom_bar(stat="identity")+coord_flip()+
  scale_fill_manual(values = b_rev)+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

########################################################## 3. Figure S1 #####################################################################
## Figure S1A: get violin plot of UMI count
ggplot(scanpy,aes(x=factor(organ),y=log10(n_counts),fill=organ,color=organ))+
  geom_violin(alpha=0.95,width=1)+
  labs(x="",y="log10(UMI Count)")+
  theme_classic()+
  scale_y_continuous(limits = c(2,5))+
  scale_fill_manual(values = a)+
  scale_color_manual(values = a)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  guides(fill=FALSE)

## Figure S1B: get violin plot of gene count
ggplot(scanpy,aes(x=factor(organ),y=log10(n_genes),fill=organ,color=organ))+
  geom_violin(alpha=0.95,width=1)+
  labs(x="",y="log10(Gene Count)")+
  theme_classic()+
  scale_y_continuous(limits = c(2,4))+
  scale_fill_manual(values = a)+
  scale_color_manual(values = a)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
  guides(fill=FALSE)

## Figure S1C: get celltype marker expression heatmap
# input celltype expression data frame normalized
bulk_gene_exp<-read.table("organ9_bulk_gene_rlogTransformation.txt",sep="\t",header=T)
bulk_gene_exp<-bulk_gene_exp[order(rownames(bulk_gene_exp)),]

# input marker list of 44 celltypes of 9 organs
marker<-read.table("organ9_celltype44_marker.txt",sep="\t")
marker1<-unique(marker)

library(pheatmap)
library("viridis")

pheatmap(bulk_gene_exp[,as.vector(marker1$V1[marker1$V1%in%colnames(bulk_gene_exp)])],scale="column",cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("dark blue","blue", "royal blue","sky blue","white","salmon","red","maroon","black"))(100))
