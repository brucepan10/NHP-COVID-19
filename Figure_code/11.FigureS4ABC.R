library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)

########################################################## 1.Figure S4A #####################################################################
load("ALL_cells_kidney.integrated.Find.RData")
Kidney_monkey=kidney.integrated.Find
meta_monkey_kidney=Kidney_monkey@meta.data
meta_monkey_kidney$celltype="NNN"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==0]="Proximal tuble S1"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==2]="Proximal tuble cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==4]="Proximal tuble S3"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==1]="Thick ascending limb cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==3]="Thick ascending limb cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==7]="Thick ascending limb cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==18]="Thick ascending limb cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==12]="Thin limb cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==6]="Principle cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==5]="Principle cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==8]="Principle cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==14]="Principle cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==9]="Distal convoluted tubule cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==19]="Distal convoluted tubule cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==10]="Connecting tubule cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==11]="Endothelial cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==13]="Intercalated cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==15]="Intercalated cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==17]="Intercalated cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==16]="Podocyte"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==20]="Mesangial cell"
Kidney_monkey$celltype=meta_monkey_kidney$celltype
Idents(Kidney_monkey)=Kidney_monkey$celltype
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  "black",
  "#FFFF00","#1CE6FF","#FF34FF",
  "#D16100"
)
pdf("FigureS4A.pdf")
DimPlot(Kidney_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("kidney")
dev.off()
########################################################## 2.Figure S4B #####################################################################
DefaultAssay(Kidney_monkey)="RNA"
deg=FindMarkers(Kidney_monkey,ident.1="Proximal tuble cell3",
                ident.2="Proximal tuble cell1",
                min.pct = 0,p_val_adj=1,logfc.threshold = 0)
deg$gene=rownames(deg)
deg$FC=exp(1)^deg$avg_logFC
deg$Log2Foldchange=log2(deg$FC)
deg_Sig_UP=subset(deg,deg$p_val_adj<0.01 & deg$FC>1.5)
deg_Sig_DOWN=subset(deg,deg$p_val_adj<0.01 & deg$FC < 2/3)
deg_Sig_NONE=subset(deg,deg$p_val_adj>0.01 | abs(deg$FC)<1.5)
deg_Sig_UP$Tag="Sig_UP"
deg_Sig_DOWN$Tag="Sig_DOWN"
deg_Sig_NONE$Tag="Sig_NONE"
total=rbind(deg_Sig_UP,deg_Sig_DOWN,deg_Sig_NONE)
total$log10_Padj=log10(total$p_val_adj)
total$log10_Padj[total$log10_Padj < -300]= -300
la=total[match(c("SLC7A13","SLC5A2","ACE2"),total$gene),]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
c=getPalette(9)
pdf("FigureS4B.pdf")
ggplot(total,aes(log2(FC),-log10_Padj))+ 
  geom_point(aes(color = Tag),size = 1)+
  geom_text_repel(data=la,aes(x=log2(FC),y=-log10_Padj,label=gene),size=3,color="green")+
  scale_color_manual(values = c(c[2], "lightgrey",c[1])) +
  geom_point(data=la,aes(log2(FC),-log10_Padj),size=3)+
  geom_hline(yintercept = 2,linetype = 2) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2)+
  theme(legend.title=element_blank())+theme(legend.background = element_blank())+theme(legend.key = element_blank())+
  theme(legend.text = element_text(size=20))+theme_bw()+theme(axis.text.x = element_text(size = 20))+theme(axis.text.y = element_text(size = 20))+
  ylim(0,310)+xlim(-3,3)
dev.off()
########################################################## 3.Figure S4C #####################################################################
pdf("Figure S4C_1.pdf")
FeaturePlot(kidney.integrated.Find,features = "SLC5A2",
            order = T)+
            scale_color_viridis()
            ggtitle("Kidney SLC5A2")
dev.off()

pdf("Figure S4C_2.pdf")
FeaturePlot(kidney.integrated.Find,features = "SLC7A13",
            order = T)+
            scale_color_viridis()
            ggtitle("Kidney SLC7A13")
dev.off()
