import numpy as np
import pandas as pd
import scanpy as sc
import os
import sys

sc.settings.verbosity = 3

sc.settings.set_figure_params(dpi=500)

adata = sc.read('Human_kidney.h5ad') # input adata of human kidney

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

sc.settings.autosave = False
########################################################## 1. Figure 5H left #####################################################################
marker = ['IL6R']
marker1 = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1,sort_order=True, size=3.5,use_raw=False)
plt.savefig('Human_kidney_ACE2_feature.pdf')

########################################################## 2. Figure S5D #####################################################################
marker = ['ACE2']
marker2 = [gene for gene in marker if (gene in adata.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker2,sort_order=True, size=3.5,use_raw=False)
plt.savefig('Human_kidney_IL6R_feature.pdf')
sc.settings.autosave = True
