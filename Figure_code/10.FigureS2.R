########################################################## Figure S2 #####################################################################
########################################################## 1. Figure S2A #####################################################################
### Clustering of Lung
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Lung_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
kidney = CreateSeuratObject(counts = count, meta.data = meta)
kidney <- subset(kidney, subset = nFeature_RNA > 200)
kidney.list <- SplitObject(object = kidney, split.by = "Batch")
for (i in 1:length(x = kidney.list)) {
    kidney.list[[i]] <- NormalizeData(object = kidney.list[[i]], verbose = FALSE)
    kidney.list[[i]] <- FindVariableFeatures(object = kidney.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- kidney.list
kidney.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
kidney.integrated <- IntegrateData(anchorset = kidney.anchors, dims = 1:30)
DefaultAssay(object = kidney.integrated) <- "integrated"
kidney.integrated <- ScaleData(object = kidney.integrated, verbose = FALSE)
kidney.integrated <- RunPCA(object = kidney.integrated, npcs = 30, verbose = FALSE)
kidney.integrated <- RunUMAP(object = kidney.integrated, reduction = "pca",
    dims = 1:30)
kidney.integrated = FindNeighbors(object = kidney.integrated,k.param=40,dims = 1:30)
kidney.integrated.Find <- FindClusters(object = kidney.integrated,resolution = 0.5)
kidney.integrated.Find <- RunTSNE(object = kidney.integrated.Find,dims = 1:30,seed.use = 2)
save(kidney.integrated.Find,file="ALL_cells_kidney.integrated.Find.RData")
inputmarkers <- FindAllMarkers(kidney.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
lung_monkey=Lung.integrated.Find
meta_monkey=lung_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="Ciliated cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="Ciliated cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==16]="Ciliated cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="Macrophage"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="Macrophage"
meta_monkey$celltype[meta_monkey$seurat_clusters==14]="Macrophage"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="Pulmonary alveolar Type1 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==10]="Pulmonary alveolar Type1 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="Pulmonary alveolar Type2 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==8]="Pulmonary alveolar Type2 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==9]="Pulmonary alveolar Type2 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==12]="Pulmonary alveolar Type2 cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==11]="Endothelial cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==17]="Pericyte"
meta_monkey$celltype[meta_monkey$seurat_clusters==13]="Fibroblast"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="Muscle cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==18]="Muscle cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="Club cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==15]="Cycling macrophage"
lung_monkey$celltype=meta_monkey$celltype

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

pdf("Cluster_Lung.pdf")
DimPlot(lung_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("lung")
dev.off()

########################################################## 2. Figure S2B #####################################################################
### Clustering of Kidney
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Kidney_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
kidney = CreateSeuratObject(counts = count, meta.data = meta)
kidney <- subset(kidney, subset = nFeature_RNA > 200)
kidney.list <- SplitObject(object = kidney, split.by = "Batch")
for (i in 1:length(x = kidney.list)) {
    kidney.list[[i]] <- NormalizeData(object = kidney.list[[i]], verbose = FALSE)
    kidney.list[[i]] <- FindVariableFeatures(object = kidney.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- kidney.list
kidney.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
kidney.integrated <- IntegrateData(anchorset = kidney.anchors, dims = 1:30)
DefaultAssay(object = kidney.integrated) <- "integrated"
kidney.integrated <- ScaleData(object = kidney.integrated, verbose = FALSE)
kidney.integrated <- RunPCA(object = kidney.integrated, npcs = 30, verbose = FALSE)
kidney.integrated <- RunUMAP(object = kidney.integrated, reduction = "pca",
    dims = 1:30)
kidney.integrated = FindNeighbors(object = kidney.integrated,k.param=40,dims = 1:30)
kidney.integrated.Find <- FindClusters(object = kidney.integrated,resolution = 0.5)
kidney.integrated.Find <- RunTSNE(object = kidney.integrated.Find,dims = 1:30,seed.use = 2)
save(kidney.integrated.Find,file="ALL_cells_kidney.integrated.Find.RData")
inputmarkers <- FindAllMarkers(kidney.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Kidney_monkey=kidney.integrated.Find
meta_monkey_kidney=Kidney_monkey@meta.data
meta_monkey_kidney$celltype="NNN"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==0]="Proximal tuble cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==2]="Proximal tuble cell"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==4]="Proximal tuble cell"
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
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==13]="Intercalated cell 1"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==15]="Intercalated cell 1"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==17]="Intercalated cell 2"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==16]="Podocyte"
meta_monkey_kidney$celltype[meta_monkey_kidney$seurat_clusters==20]="Stormal cell"
Kidney_monkey$celltype=meta_monkey_kidney$celltype

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

pdf("Cluster_Kidney.pdf")
DimPlot(Kidney_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("kidney")
dev.off()

########################################################## 3. Figure S2C #####################################################################
### Clustering of Liver
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Liver_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
liver = CreateSeuratObject(counts = count, meta.data = meta)
liver <- subset(liver, subset = nFeature_RNA > 200)
liver.list <- SplitObject(object = liver, split.by = "Batch")
for (i in 1:length(x = liver.list)) {
    liver.list[[i]] <- NormalizeData(object = liver.list[[i]], verbose = FALSE)
    liver.list[[i]] <- FindVariableFeatures(object = liver.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- liver.list
liver.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
liver.integrated <- IntegrateData(anchorset = liver.anchors, dims = 1:30)
DefaultAssay(object = liver.integrated) <- "integrated"
liver.integrated <- ScaleData(object = liver.integrated, verbose = FALSE)
liver.integrated <- RunPCA(object = liver.integrated, npcs = 30, verbose = FALSE)
liver.integrated <- RunUMAP(object = liver.integrated, reduction = "pca",
    dims = 1:30)
liver.integrated = FindNeighbors(object = liver.integrated,k.param=40,dims = 1:30)
liver.integrated.Find <- FindClusters(object = liver.integrated,resolution = 0.5)
liver.integrated.Find <- RunTSNE(object = liver.integrated.Find,dims = 1:30,seed.use = 2)
save(liver.integrated.Find,file="ALL_cells_liver.integrated.Find.RData")
inputmarkers <- FindAllMarkers(liver.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Liver_monkey=LiverL.integrated.Find
meta_monkey_Liver=Liver_monkey@meta.data
meta_monkey_Liver$celltype="NNN"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==0]="Hepatocyte"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==1]="Hepatocyte"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==2]="Hepatocyte"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==3]="Hepatocyte"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==4]="Endothelial cell"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==9]="Endothelial cell"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==5]="NK T"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==6]="Kupffer cell"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==7]="Heptic stellate cell"
meta_monkey_Liver$celltype[meta_monkey_Liver$seurat_clusters==8]="Cholangiocyte"
Liver_monkey$celltype=meta_monkey_Liver$celltype

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

pdf("Cluster_Liver.pdf")
DimPlot(Liver_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25[-c(3,7)])+
  ggtitle("liver")
dev.off()

########################################################## 4. Figure S2D #####################################################################
### Clustering of PBMC
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("PBMC_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
liver = CreateSeuratObject(counts = count, meta.data = meta)
liver <- subset(liver, subset = nFeature_RNA > 200)
liver.list <- SplitObject(object = liver, split.by = "Batch")
for (i in 1:length(x = liver.list)) {
    liver.list[[i]] <- NormalizeData(object = liver.list[[i]], verbose = FALSE)
    liver.list[[i]] <- FindVariableFeatures(object = liver.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- liver.list
liver.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
liver.integrated <- IntegrateData(anchorset = liver.anchors, dims = 1:30)
DefaultAssay(object = liver.integrated) <- "integrated"
liver.integrated <- ScaleData(object = liver.integrated, verbose = FALSE)
liver.integrated <- RunPCA(object = liver.integrated, npcs = 30, verbose = FALSE)
liver.integrated <- RunUMAP(object = liver.integrated, reduction = "pca",
    dims = 1:30)
liver.integrated = FindNeighbors(object = liver.integrated,k.param=40,dims = 1:30)
liver.integrated.Find <- FindClusters(object = liver.integrated,resolution = 0.5)
liver.integrated.Find <- RunTSNE(object = liver.integrated.Find,dims = 1:30,seed.use = 2)
PBMC.integrated.Find <- liver.integrated.Find
save(PBMC.integrated.Find,file="ALL_cells_PBMC.integrated.Find.RData")
inputmarkers <- FindAllMarkers(liver.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
PBMC_monkey=PBMC.integrated.Find
meta_monkey=PBMC_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="B cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="CD4 T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="CD4 T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="CD8 Naive T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="CD4 T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="CD8+ memory T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="B cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="CD8 Naive T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==8]="CD4 T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==9]="NK"
meta_monkey$celltype[meta_monkey$seurat_clusters==10]="Dentritc cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==11]="CD4 T cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==12]="Platelet /CD14+ Monocyte"
meta_monkey$celltype[meta_monkey$seurat_clusters==13]="CD16+ Monocyte"
PBMC_monkey$celltype=meta_monkey$celltype

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

pdf("Cluster_PBMC.pdf")
DimPlot(PBMC_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("PBMC")
dev.off()

########################################################## 5. Figure S2E #####################################################################
### Clustering of Neocortex
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Neocortex_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
Neocortex = CreateSeuratObject(counts = count, meta.data = meta)
Neocortex <- subset(Neocortex, subset = nFeature_RNA > 200)
Neocortex.list <- SplitObject(object = Neocortex, split.by = "Batch")
for (i in 1:length(x = Neocortex.list)) {
    Neocortex.list[[i]] <- NormalizeData(object = Neocortex.list[[i]], verbose = FALSE)
    Neocortex.list[[i]] <- FindVariableFeatures(object = Neocortex.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- Neocortex.list
Neocortex.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Neocortex.integrated <- IntegrateData(anchorset = Neocortex.anchors, dims = 1:30)
DefaultAssay(object = Neocortex.integrated) <- "integrated"
Neocortex.integrated <- ScaleData(object = Neocortex.integrated, verbose = FALSE)
Neocortex.integrated <- RunPCA(object = Neocortex.integrated, npcs = 30, verbose = FALSE)
Neocortex.integrated <- RunUMAP(object = Neocortex.integrated, reduction = "pca",
    dims = 1:30)
Neocortex.integrated = FindNeighbors(object = Neocortex.integrated,k.param=40,dims = 1:30)
Neocortex.integrated.Find <- FindClusters(object = Neocortex.integrated,resolution = 0.5)
Neocortex.integrated.Find <- RunTSNE(object = Neocortex.integrated.Find,dims = 1:30,seed.use = 2)
save(Neocortex.integrated.Find,file="ALL_cells_Neocortex.integrated.Find.RData")
inputmarkers <- FindAllMarkers(Neocortex.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Neocortex_monkey=Neocortex.integrated.Find
meta_monkey=Neocortex_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==15]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==12]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==14]="Excitatory neuron"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="Astrocyte"
meta_monkey$celltype[meta_monkey$seurat_clusters==8]="Oligodendrocyte"
meta_monkey$celltype[meta_monkey$seurat_clusters==16]="Microglia"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="SST"
meta_monkey$celltype[meta_monkey$seurat_clusters==11]="Oligodendrocyte precursor cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==13]="SV2C"
meta_monkey$celltype[meta_monkey$seurat_clusters==10]="VIP"
meta_monkey$celltype[meta_monkey$seurat_clusters==9]="PVALB"
Neocortex_monkey$celltype=meta_monkey$celltype

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

pdf("Cluster_Neocortex.pdf")
DimPlot(Neocortex_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("neocortex")
dev.off()

########################################################## 6. Figure S2F #####################################################################
### Clustering of Parotid
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Parotid_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
liver = CreateSeuratObject(counts = count, meta.data = meta)
liver <- subset(liver, subset = nFeature_RNA > 200)
liver.list <- SplitObject(object = liver, split.by = "Batch")
for (i in 1:length(x = liver.list)) {
    liver.list[[i]] <- NormalizeData(object = liver.list[[i]], verbose = FALSE)
    liver.list[[i]] <- FindVariableFeatures(object = liver.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- liver.list
liver.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
liver.integrated <- IntegrateData(anchorset = liver.anchors, dims = 1:30)
DefaultAssay(object = liver.integrated) <- "integrated"
liver.integrated <- ScaleData(object = liver.integrated, verbose = FALSE)
liver.integrated <- RunPCA(object = liver.integrated, npcs = 30, verbose = FALSE)
liver.integrated <- RunUMAP(object = liver.integrated, reduction = "pca",
    dims = 1:30)
liver.integrated = FindNeighbors(object = liver.integrated,k.param=40,dims = 1:30)
liver.integrated.Find <- FindClusters(object = liver.integrated,resolution = 0.5)
liver.integrated.Find <- RunTSNE(object = liver.integrated.Find,dims = 1:30,seed.use = 2)
Parotid.integrated.Find <- liver.integrated.Find
save(Parotid.integrated.Find,file="ALL_cells_Parotid.integrated.Find.RData")
inputmarkers <- FindAllMarkers(liver.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Parotid_monkey=Parotid.integrated.Find
meta_monkey=Parotid_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="Serous acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==12]="Serous acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="Serous acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==11]="Serous acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="Striated duct cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==9]="Striated duct cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==10]="Striated duct cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="Intercalated duct cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="Mucous acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="Myoepithelial cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="Myoepithelial cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==8]="Macrophage"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="Stromal cell"
Parotid_monkey$celltype=meta_monkey$celltype

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

pdf("Cluster_Parotid.pdf")
DimPlot(Parotid_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("parotid")
dev.off()

########################################################## 7. Figure S2G #####################################################################
### Clustering of Aorta
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Aorta_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]

Aorta = CreateSeuratObject(counts = count, meta.data = meta)
Aorta <- subset(Aorta, subset = nFeature_RNA > 200)
Aorta.list <- SplitObject(object = Aorta, split.by = "Batch")
for (i in 1:length(x = Aorta.list)) {
    Aorta.list[[i]] <- NormalizeData(object = Aorta.list[[i]], verbose = FALSE)
    Aorta.list[[i]] <- FindVariableFeatures(object = Aorta.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- Aorta.list
Aorta.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Aorta.integrated <- IntegrateData(anchorset = Aorta.anchors, dims = 1:30)
DefaultAssay(object = Aorta.integrated) <- "integrated"
Aorta.integrated <- ScaleData(object = Aorta.integrated, verbose = FALSE)
Aorta.integrated <- RunPCA(object = Aorta.integrated, npcs = 30, verbose = FALSE)
Aorta.integrated <- RunUMAP(object = Aorta.integrated, reduction = "pca",
    dims = 1:30)
Aorta.integrated = FindNeighbors(object = Aorta.integrated,k.param=40,dims = 1:30)
Aorta.integrated.Find <- FindClusters(object = Aorta.integrated,resolution = 0.5)
Aorta.integrated.Find <- RunTSNE(object = Aorta.integrated.Find,dims = 1:30,seed.use = 2)
save(Aorta.integrated.Find,file="ALL_cells_Aorta.integrated.Find.RData")
inputmarkers <- FindAllMarkers(Aorta.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Aorta_monkey=Aorta.integrated.Find
meta_monkey_Aorta=Aorta_monkey@meta.data
meta_monkey_Aorta$celltype="NNN"
meta_monkey_Aorta$celltype[meta_monkey_Aorta$seurat_clusters==0]="Smooth muscle cell"
meta_monkey_Aorta$celltype[meta_monkey_Aorta$seurat_clusters==3]="Smooth muscle cell"
meta_monkey_Aorta$celltype[meta_monkey_Aorta$seurat_clusters==1]="Adipocyte"
meta_monkey_Aorta$celltype[meta_monkey_Aorta$seurat_clusters==4]="Endothelial cell"
meta_monkey_Aorta$celltype[meta_monkey_Aorta$seurat_clusters==2]="Myofibroblast"
Aorta_monkey$celltype=meta_monkey_Aorta$celltype

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

pdf("Cluster_Aorta.pdf")
DimPlot(Aorta_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("Aorta")
dev.off()


########################################################## 8. Figure S2H #####################################################################
### Clustering of Thyroid
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Thyroid_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
Thyroid = CreateSeuratObject(counts = count, meta.data = meta)
Thyroid <- subset(Thyroid, subset = nFeature_RNA > 200)
Thyroid.list <- SplitObject(object = Thyroid, split.by = "Batch")
for (i in 1:length(x = Thyroid.list)) {
    Thyroid.list[[i]] <- NormalizeData(object = Thyroid.list[[i]], verbose = FALSE)
    Thyroid.list[[i]] <- FindVariableFeatures(object = Thyroid.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- Thyroid.list
Thyroid.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Thyroid.integrated <- IntegrateData(anchorset = Thyroid.anchors, dims = 1:30)
DefaultAssay(object = Thyroid.integrated) <- "integrated"
Thyroid.integrated <- ScaleData(object = Thyroid.integrated, verbose = FALSE)
Thyroid.integrated <- RunPCA(object = Thyroid.integrated, npcs = 30, verbose = FALSE)
Thyroid.integrated <- RunUMAP(object = Thyroid.integrated, reduction = "pca",
    dims = 1:30)
Thyroid.integrated = FindNeighbors(object = Thyroid.integrated,k.param=40,dims = 1:30)
Thyroid.integrated.Find <- FindClusters(object = Thyroid.integrated,resolution = 0.5)
Thyroid.integrated.Find <- RunTSNE(object = Thyroid.integrated.Find,dims = 1:30,seed.use = 2)
save(Thyroid.integrated.Find,file="ALL_cells_Thyroid.integrated.Find.RData")
inputmarkers <- FindAllMarkers(Thyroid.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Thyroid_monkey=Thyroid.integrated.Find
meta_monkey=Thyroid_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="Follicular cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="Follicular cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="Follicular cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="Follicular cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="Endothelial cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="Stromal cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="Muscle cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="Adipocyte"
Thyroid_monkey$celltype=meta_monkey$celltype

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

pdf("Cluster_Thyroid.pdf")
DimPlot(Thyroid_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("thyroid")
dev.off()

########################################################## 9. Figure S2I #####################################################################
### Clustering of Pancreas
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)

count=as.data.frame(fread("Pancreas_rawcount.txt"))
rownames(count)=count[,1]
count=count[,-1]
meta=as.data.frame(colnames(count))
meta$Batch=unlist(lapply(strsplit(as.character(meta$`colnames(count)`),":"),"[",1))
rownames(meta)=meta[,1]
meta=meta[2]
Pancreas = CreateSeuratObject(counts = count, meta.data = meta)
Pancreas <- subset(Pancreas, subset = nFeature_RNA > 200)
Pancreas.list <- SplitObject(object = Pancreas, split.by = "Batch")
for (i in 1:length(x = Pancreas.list)) {
    Pancreas.list[[i]] <- NormalizeData(object = Pancreas.list[[i]], verbose = FALSE)
    Pancreas.list[[i]] <- FindVariableFeatures(object = Pancreas.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- Pancreas.list
Pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Pancreas.integrated <- IntegrateData(anchorset = Pancreas.anchors, dims = 1:30)
DefaultAssay(object = Pancreas.integrated) <- "integrated"
Pancreas.integrated <- ScaleData(object = Pancreas.integrated, verbose = FALSE)
Pancreas.integrated <- RunPCA(object = Pancreas.integrated, npcs = 30, verbose = FALSE)
library("reticulate")
use_condaenv(condaenv="r-reticulate", conda="/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/Anaconda2/bin/conda")
Pancreas.integrated <- RunUMAP(object = Pancreas.integrated, reduction = "pca",
    dims = 1:30)
Pancreas.integrated = FindNeighbors(object = Pancreas.integrated,k.param=40,dims = 1:30)
Pancreas.integrated.Find <- FindClusters(object = Pancreas.integrated,resolution = 0.5)
Pancreas.integrated.Find <- RunTSNE(object = Pancreas.integrated.Find,dims = 1:30,seed.use = 2)
save(Pancreas.integrated.Find,file="ALL_cells_Pancreas.integrated.Find.RData")
inputmarkers <- FindAllMarkers(Pancreas.integrated.Find)
write.table(inputmarkers,"ALL_cells_clustermarkers.txt",sep="\t",quote = FALSE)

### Cell type annotation
Pancreas_monkey=Pancreas.integrated.Find
meta_monkey=Pancreas_monkey@meta.data
meta_monkey$celltype="NNN"
meta_monkey$celltype[meta_monkey$seurat_clusters==0]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==1]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==2]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==4]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==5]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==3]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==6]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==11]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==13]="Acinar cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==8]="unknown"
meta_monkey$celltype[meta_monkey$seurat_clusters==7]="Ductal cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==10]="Beta cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==9]="Stromal cell"
meta_monkey$celltype[meta_monkey$seurat_clusters==12]="Alpha cell"
Pancreas_monkey$celltype=meta_monkey$celltype
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

pdf("Cluster_Pancreas.pdf")
DimPlot(Pancreas_monkey,group.by = "celltype",label = T)+NoLegend()+
  scale_color_manual(values = c25)+
  ggtitle("pancreas")
dev.off()

