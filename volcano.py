#pip 3 install bioinfokit
#pip3 install -U scikit-learn
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp 
import seaborn

from bioinfokit import analys, visuz
#df = analys.get_data('diffstats_allregionsold.txt').data
df = pd.read_csv('diffstats_allregionsold_volcano.txt', sep = "\t")

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
#VOLCANO PLOT, remember to adjust log2FC threshold 
visuz.gene_exp.volcano(df=df, lfc='log2FC', pv = 'FDR', geneid = 'Gene name', lfc_thr=0.2, pv_thr=0.05, \
    plotlegend=True, legendpos='upper right', legendanchor=(1.46,1), color=("#E10600FF", "grey", "#00239CFF"), \
    sign_line=True, genenames=('RGS1', 'CD163', 'LINC01094', 'HLA-DRA', 'ADAMTS2', 'SST', 'BDNF', 'PCDH8', 'MIR7-3HG', 'CALB1'), \
        gfont=4, dim=(5,5), xlm=(-1.0,1.0,0.2), figtype='pdf', dotsize=4, axylabel = '-log10(adjusted p-value)'.translate(SUB))

#y-axis is -log10(FDR), for example y=7 means that FDR = 10^-7 = 0.0000001, because -log(0.0000001)=7 