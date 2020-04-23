# NHP-COVID-19
Source code of manuscript: [Single-cell atlas of a non-human primate reveals new pathogenic mechanisms of COVID-19](https://www.biorxiv.org/content/10.1101/2020.04.10.022103v1)


**Requirements**  
**1.** These scripts have been tested on Windows or Linux, but should be supported by Mac OSX.  
**2.** Users should have R version 3.5 or higher, python version 3.6, and several packages as indicated in the scripts.


**Descriptions of the scripts**  
01.Clustering_organ9_Figure2AB_4D.py  
   Clustering analysis of all single cells within our dataset;
   UMAP projection of ACE2, TMPRSS2, TMEM27, IDO2, DNAJC12 and ANPEP expression in all single cells.
   
02.Figure1BC_FigureS1.R  
   UMAP visualization of all single cells from the dataset colored by tissue/organ and 44 clusters in the first round of clustering; 
   Quality control of the single-cell RNA-seq libraries.  
   
03.Figure2CDE.R  
   UMAP projection of ACE2<sup>+</sup>/TMPRSS2<sup>+</sup> cells; 
   Bubble plots showing the level of expression of TMPRSS2 and ACE2 genes and the ratio of expressing cells in the indicated cell types; 
   Barplot indicating the percentage of ACE2 and TMPRSS2 expressing cells within each cell cluster.  
   
04.Figure3.R  
   Comparative analysis of ACE2 and TMPRSS2 expression between monkey and human.
   
05.Figure4ABC.R  
   Co-expression analysis of ACE2 in monkey tissues.
   
06.Figure5BCDG_FigureS5AB.R  
   Chromatin accessibility analysis reveals epigenetic regulation of ACE2 in kidney.
   
07.Figure5E.R  
   Ratio of ACE2<sup>+</sup> cells in each cell type of kidney.
   
08.Figure5H_left_FigureS5D.py  
   UMAP projection of IL6R and ACE2 expression in all kidney single cells in human.  
09.Figure5H_right.R  
   UMAP projection of cells with IL6R<sup>+</sup>/ACE2<sup>+</sup> cells in all kidney single cells in human. 
   
10.FigureS2.R  
   Clustering analysis of cells from each organ.
   
11.FigureS3.R  
   ACE2 and TMPRSS2 expression in each tissue/organ. 
   
12.FigureS4ABC.R  
   Specific accessible chromatin in each cluster of kidney.  


**Data availability**  
All raw sequencing data will be made available upon request at the peer review stage, and accession codes will be available before publication.



