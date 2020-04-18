library(data.table)
library(Matrix)
library(Signac)
library(Seurat)
library(reticulate)
library(ggplot2)
library(RColorBrewer)
library(JASPAR2018)
library(viridis)
set.seed(1234)

########################################################## Cluster ATAC & Integrate with RNA data #####################################################################
raw_count=as.data.frame(fread("Kidney.count_peak_Monkey.matrix"))
raw_meta1=read.csv("Kidney1.QCstats.csv")
raw_meta2=read.csv("Kidney2.QCstats.csv")
raw_meta=rbind(raw_meta1,raw_meta2)
peak=as.data.frame(raw_count$peak)
peak$tag1=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",1))
peak$tag2=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",2))
peak$tag3=unlist(lapply(strsplit(as.character(peak[,1]),"_"),"[",3))
peak$Tag1=paste(peak$tag1,peak$tag2,sep = ":")
peak$Tag=paste(peak$Tag1,peak$tag3,sep = "-")
rownames(raw_count)=peak$Tag
raw_count=raw_count[,-1]
my_counts=as(as.matrix(raw_count), "dgCMatrix")
rownames(raw_meta)=raw_meta[,1]
raw_meta=raw_meta[,-1]
my_meta=raw_meta
kidney <- CreateSeuratObject(
  counts = my_counts,
  assay = 'peaks',
  project = 'Monkey',
  min.cells = 10,
  meta.data = my_meta
)
kidney@meta.data$uniqueNuclearFrags=log10(kidney@meta.data$uniqueNuclearFrags)
########################################################## 1. Figure S5AB #####################################################################
pdf("FigureS5A.pdf")
VlnPlot(
  object = kidney,
  features = c('uniqueNuclearFrags'),
  pt.size = 0.1) + NoLegend()+
  scale_color_brewer(palette = "Set1")
dev.off()
pdf("FigureS5B.pdf")
VlnPlot(
  object = kidney,
  features = c('tssProportion'),
  pt.size = 0.1) + NoLegend()+
  scale_color_brewer(palette = "Set1")
dev.off()
kidney <- subset(kidney, subset = uniqueNuclearFrags > 3 & tssProportion > 0.1)
save(kidney,file = "Signac_Kidney.RDS")
###ATAC data
load("Signac_Kidney.RDS")
kidney.atac=kidney
input=as.data.frame(fread("Kidney.count_promoter_Monkey.matrix"))
tss=read.table("promoter_4k_Macaca_fascicularis.bed",header = F)
tss$Name=unlist(lapply(strsplit(as.character(tss$V4),":"),"[",1))
tss$Name=make.unique(tss$Name)
rownames(input)=tss$Name
input=input[,-1]
cell=colnames(kidney.atac@assays$peaks)
sub_input=input[,match(cell,colnames(input))]
kidney.atac[["ACTIVITY"]] <- CreateAssayObject(counts = sub_input)
DefaultAssay(kidney.atac) <- "ACTIVITY"
kidney.atac <- FindVariableFeatures(kidney.atac)
kidney.atac <- NormalizeData(kidney.atac)
kidney.atac <- ScaleData(kidney.atac)
DefaultAssay(kidney.atac) <- "peaks"
VariableFeatures(kidney.atac) <- names(which(Matrix::rowSums(kidney.atac) > 100))
kidney.atac <- RunLSI(kidney.atac, n = 50, scale.max = NULL)
kidney.atac <- RunUMAP(kidney.atac, reduction = "lsi", dims = 1:50)
###RNA data
load("ALL_cells_kidney.integrated.Find.RData")
kidney.rna=kidney.integrated.Find
###transfer anchors
transfer.anchors <- FindTransferAnchors(reference = kidney.rna, query = kidney.atac, 
                                        features = VariableFeatures(object = kidney.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", 
                                        reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = kidney.rna$celltype, 
                                     weight.reduction = kidney.atac[["lsi"]])
kidney.atac <- AddMetaData(kidney.atac, metadata = celltype.predictions)
kidney.atac.filtered <- subset(kidney.atac, subset = prediction.score.max > 0)
kidney.atac.filtered$predicted.id <- factor(kidney.atac.filtered$predicted.id, levels = sort(unique(kidney.rna$celltype)))  # to make the colors match
###Co-embedding
genes.use <- VariableFeatures(kidney.rna)
refdata <- GetAssayData(kidney.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = kidney.atac[["lsi"]])
kidney.atac[["integrated"]] <- imputation
coembed <- merge(x = kidney.rna, y = kidney.atac)
DefaultAssay(coembed) <- "integrated"
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)
meta_data=coembed@meta.data
meta_data$tech="N"
meta_data$name=rownames(meta_data)
meta_data$tech[grepl("BC",meta_data$name)]<-"DNBelab C4 scATAC"
meta_data$tech[!(grepl("BC",meta_data$name))]<-"DNBelab C4 scRNA"
coembed@meta.data$tech=meta_data$tech
########################################################## 2. Figure 5B #####################################################################
pdf("Figure5B.pdf")
DimPlot(coembed,group.by = "tech")+
  scale_color_manual(values = c("#35B779FF","#440154FF"))
dev.off()
########################################################## 3. Figure 5C #####################################################################
coembed_ATAC=subset(coembed,subset=tech=="DNBelab C4 scATAC",invert=F)
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
pdf("Figure5C.pdf")
DimPlot(coembed_ATAC,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)
dev.off()
########################################################## 4. Figure 5D #####################################################################
species <- 9606
opts <- list()
opts["species"] <- species
out  <- TFBSTools::getMatrixSet(JASPAR2018, opts)
names(out) <- paste(names(out), TFBSTools::name(out),sep = "_")
motif=out
seq1=read.table("seq1.txt",sep = "\t")
seq2=read.table("seq2.txt",sep = "\t")
seq3=read.table("seq3.txt",sep = "\t")
seq4=read.table("seq4.txt",sep = "\t")
motif_ix_1 <- matchMotifs(motif, as.character(seq1[2,1])) 
motif_ix_2 <- matchMotifs(motif, as.character(seq2[2,1])) 
motif_ix_3 <- matchMotifs(motif, as.character(seq3[2,1])) 
motif_ix_4 <- matchMotifs(motif, as.character(seq4[2,1])) 
r1=t(as.matrix(motifMatches(motif_ix_1)))
r2=t(as.matrix(motifMatches(motif_ix_2)))
r3=t(as.matrix(motifMatches(motif_ix_3)))
r4=t(as.matrix(motifMatches(motif_ix_4)))
out=as.data.frame(cbind(r1,r2,r3,r4))
colnames(out)=c("Peak_139736","Peak_139735","Peak_139734","Peak_139733")
out$TF=rownames(out)
out$count=apply(out,1,function(x) sum(x=="TRUE"))
sub_out=subset(out,out$count>0)
write.table(sub_out,"ACE2_peak_Motif.xls",sep = "\t",quote = FALSE,row.names = FALSE)
########################################################## 5. Figure 5G #####################################################################
load("ALL_cells_kidney.integrated.Find.RData")
pdf("Figure 5G_1.pdf")
FeaturePlot(kidney.integrated.Find,features = "IL6R",
            cols = viridis(50,option = "D"),order = T,max.cutoff = 'q95',
            min.cutoff = 'q5')+
  ggtitle("Kidney IL6R")
dev.off()

kidney_RNA=as.data.frame(kidney.integrated.Find@assays$RNA@data)
kidney_RNA_sub=as.data.frame(t(kidney_RNA[match(c("ACE2","IL6R"),rownames(kidney_RNA)),]))
kidney_coor=as.data.frame(Embeddings(object = kidney.integrated.Find,reduction = "umap"))
kidney_res=cbind(kidney_coor,kidney_RNA_sub)
kidney_res$Positive=kidney_res$ACE2>0 & kidney_res$IL6R>0
sum_both=subset(kidney_res,kidney_res$Positive=="TRUE")

pdf("Figure 5G_2.pdf")
ggplot()+
  geom_point(data=kidney_res,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=0.0001,colour="#0D0887FF")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=25),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_point(data=sum_both,aes(x=UMAP_1, y =UMAP_2),alpha=1, size=1,colour="#F0F921FF")+
  ggtitle("Kidney ACE2/IL6R+")+
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank()) 
dev.off()