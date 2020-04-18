########################################################## 1. Clustering of 9 organs #####################################################################
import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys

os.makedirs(sys.argv[1]) #create the output directory
os.chdir(sys.argv[1]) 

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/organ9_concatenate.h5ad'  # the file that will store the analysis results
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300, frameon=False)  # low dpi (dots per inch) yields small inline figures

### input gene expression matrix of each organ
#Kidney
data1 = sc.read_csv('Kidney_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata1 = data1.T
adata1
adata1.X.shape
#Liver
data2 = sc.read_csv('Liver_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata2 = data2.T
adata2
adata2.X.shape
#Lung
data3 = sc.read_csv('Lung_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata3 = data3.T
adata3
adata3.X.shape
#Pancreas
data4 = sc.read_csv('Pancreas_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata4 = data4.T
adata4
adata4.X.shape
#Parotid
data5 = sc.read_csv('Parotid_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata5 = data5.T
adata5
adata5.X.shape
#PBMC
data6 = sc.read_csv('PBMC_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata6 = data6.T
adata6
adata6.X.shape
#Aorta
data7 = sc.read_csv('Aorta_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata7 = data7.T
adata7
adata7.X.shape
#Thyroid
data8 = sc.read_csv('Thyroid_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata8 = data8.T
adata8
adata8.X.shape
#Neocortex
data9 = sc.read_csv('Neocortex_rawcount.txt',delimiter='\t', first_column_names=None, dtype='float32')
adata9 = data9.T
adata9
adata9.X.shape

adata = adata1.concatenate(adata2,adata3,adata4,adata5,adata6,adata7,adata8,adata9, join='outer')
adata.X = np.array(pd.DataFrame(adata.X).fillna(0))

adata

### save raw count matrix of 9 organs
pd.DataFrame(data=adata.X.T, index=adata.var_names,columns=adata.obs_names).to_csv('raw_mat.csv',sep="\t",float_format='%.0f')

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

### Preprocessing
#Show those genes that yield the highest fraction of counts in each single cells, across all cells.
sc.pl.highest_expr_genes(adata, n_top=20)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = ['ND5','COX3','ATP8','COX1','ND6','ND3','ND4L','COX2','ND1','CYTB','ATP6','ND4','ND2']
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)

sc.settings.autosave = False
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
plt.savefig('./figures/03.cm_scatter.pdf')
sc.pl.scatter(adata, x='n_counts', y='n_genes')
plt.savefig('./figures/04.cg_scatter.pdf')
sc.settings.autosave = True
adata

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
### save normalize matrix of 9 organs
pd.DataFrame(data=adata.X.T, index=adata.var_names,columns=adata.obs_names).to_csv('normalize_mat.csv',sep="\t")
sc.pp.log1p(adata)
### save normalize and ln(x+1) matrix of 9 organs
pd.DataFrame(data=adata.X.T, index=adata.var_names,columns=adata.obs_names).to_csv('normalize_log_mat.csv',sep="\t")
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var['highly_variable']]

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)
### save scale matrix of 9 organs
pd.DataFrame(data=adata.X.T, index=adata.var_names,columns=adata.obs_names).to_csv('scale_mat.csv',sep="\t")

### Principal component analysis
sc.tl.pca(adata, n_comps=50, random_state=int(9), svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

adata.write(results_file)
adata

### Computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=int(9))

### Embedding the neighborhood graph
sc.tl.umap(adata, random_state=int(9))

### plot classical biomarkers
sc.settings.autosave = False
#Neocortex
marker = ['SLC17A7','GAD1','MOG','PDGFRA','SLC1A3','AIF1','CLDN5','KCNJ8','PTPRC','ADGRE1','C1QB','FLT1','PECAM1','PDGFRB','PVALB','VIP','SST']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Neocortex_marker.pdf')

marker = ['SLC17A7','GAD1','MOG','PDGFRA','SLC1A3','AIF1','CLDN5','KCNJ8','PTPRC','ADGRE1','C1QB','FLT1','PECAM1','PDGFRB','PVALB','VIP','SST']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Neocortex_marker.pdf')

#Kidney
marker = ['PLA2R1', 'EGF', 'LRP2', 'FXYD4', 'EMCN', 'SLC8A1', 'SLC12A3', 'FOXI1','FLT1', 'VWF']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Kidney_marker.pdf')

marker = ['PLA2R1', 'EGF', 'LRP2', 'FXYD4', 'EMCN', 'SLC8A1', 'SLC12A3', 'FOXI1','FLT1', 'VWF']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Kidney_marker.pdf')

#Lung
marker = ['FOXJ1', 'SCGB1A1', 'AGER', 'SFTPB', 'SFTPC','SCGB3A2', 'VWF','DCN', 'PDGFRB','PTPRC', 'MK167','FLT1', 'EMCN']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Lung_marker.pdf')

marker = ['FOXJ1', 'SCGB1A1', 'AGER', 'SFTPB', 'SFTPC','SCGB3A2', 'VWF','DCN', 'PDGFRB','PTPRC', 'MK167','FLT1', 'EMCN']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Lung_marker.pdf')

#Liver
marker = ['SEC16B', 'KRT19', 'PDGFRB', 'DCN','FLT1', 'FCN2','VWF','PECAM1','PTPRC','CD163','GZMB','EMCN']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Liver_marker.pdf')

marker = ['SEC16B', 'KRT19', 'PDGFRB', 'DCN','FLT1', 'FCN2','VWF','PECAM1','PTPRC','CD163','GZMB','EMCN']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Liver_marker.pdf')

#Pancreas
marker = ['GCG', 'INS', 'KRT19', 'DCN','FLT1', 'VWF', 'EMCN', 'PTPRC']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Pancreas_marker.pdf')

marker = ['GCG', 'INS', 'KRT19', 'DCN','FLT1', 'VWF', 'EMCN', 'PTPRC']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Pancreas_marker.pdf')

#Parotid
marker = ['PRR27', 'SLC17A4', 'ATP6V1B1', 'DCN', 'MYLK', 'FLT1', 'VWF', 'EMCN', 'PTPRC']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_Parotid_marker.pdf')

marker = ['PRR27', 'SLC17A4', 'ATP6V1B1', 'DCN', 'MYLK', 'FLT1', 'VWF', 'EMCN', 'PTPRC']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_Parotid_marker.pdf')

#PBMC
marker = ['MS4A1', 'GZMA', 'GZMK', 'GZMB', 'LEF1', 'PF4', 'FCGR3A']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker)
plt.savefig('./figures/raw_gene_expression_PBMC_marker.pdf')

marker = ['MS4A1', 'GZMA', 'GZMK', 'GZMB', 'LEF1', 'PF4', 'FCGR3A']
marker = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker, use_raw=False)
plt.savefig('./figures/scaled_gene_expression_PBMC_marker.pdf')

sc.settings.autosave = True

### Clustering the neighborhood graph
sc.tl.louvain(adata,resolution=float(4),random_state=int(9)) # set resolution
sc.pl.umap(adata, color=['louvain'])
adata.write(results_file)

### Finding marker genes
sc.settings.verbosity = 2  # reduce the verbosity

sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save='.pdf')

### Get a table with the scores and groups.
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
## top25_ranked_genes_per_cluster 
trgpc = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores', 'logfoldchanges','pvals_adj']}).head(25)
print(trgpc)
trgpc.to_csv(path_or_buf="./write/top25_ranked_genes_per_cluster.csv")

### Get cluster, umap coordinates and other information of all cells
adata.obs[['n_genes','percent_mito','n_counts','louvain']].to_csv('./write/Scanpy_cell_info.csv')
adata.obsm.to_df()[['X_pca1','X_pca2','X_umap1','X_umap2']].to_csv('./write/Scanpy_X_pca_umap.csv')

adata
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

########################################################## 2. Figure 2A-B 4D #####################################################################
### figure2AB
sc.settings.autosave = False
marker = ['ACE2', 'TMPRSS2']
marker = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker,sort_order=True,size=3)
plt.savefig('Figure2A-B.pdf')

### figure4D
marker = ['TMEM27','IDO2','DNAJC12','ANPEP']
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1,sort_order=True,size=3.5)
plt.savefig('Figure4D.pdf')
sc.settings.autosave = True