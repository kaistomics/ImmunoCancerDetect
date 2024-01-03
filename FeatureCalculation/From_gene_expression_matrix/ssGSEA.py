import os
import csv
import sys
import numpy as np
import pandas as pd
from datetime import datetime
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.plot import heatmap

inputmatrix = sys.argv[1] 
genes = sys.argv[2] # kegg_immune, kegg_cancer, go_bp_immune, go_bp_cancer

from gseapy import Biomart
bm = Biomart()

start = datetime.now()

## 01. open gene expression matrix

input = pd.read_csv(inputmatrix, sep='\t')
input.set_index('Symbol', inplace=True)
print("INPUT format ...")
print(input.iloc[:2, :2])

## 02. download enrichr library, prepare geneset dictionary

names = gp.get_library_name()
kegg_2021 = gp.get_library(name='KEGG_2021_Human', organism='Human')
go_bp_2021 = gp.get_library(name='GO_Biological_Process_2021', organism='Human')

kegg_immune = {}
common_immune = ['Antigen processing and presentation',
 'B cell receptor signaling pathway',
 'Chemokine signaling pathway',
 'Complement and coagulation cascades',
 'Cytosolic DNA-sensing pathway',
 'Fc epsilon RI signaling pathway',
 'Fc gamma R-mediated phagocytosis',
 'Hematopoietic cell lineage',
 'IL-17 signaling pathway',
 'Intestinal immune network for IgA production',
 'Leukocyte transendothelial migration',
 'NOD-like receptor signaling pathway',
 'Natural killer cell mediated cytotoxicity',
 'Neutrophil extracellular trap formation',
 'Platelet activation',
 'RIG-I-like receptor signaling pathway',
 'T cell receptor signaling pathway',
 'Th1 and Th2 cell differentiation',
 'Toll-like receptor signaling pathway']
for term in common_immune:
    kegg_immune[term] = kegg_2021[term]

kegg_cancer = {}
common_cancer = ['Central carbon metabolism in cancer', 'Endometrial cancer', 'PD-L1 expression and PD-1 checkpoint pathway in cancer',
        'Transcriptional misregulation in cancer']
for term in common_cancer:
    kegg_cancer[term] = kegg_2021[term]

go_bp_immune = {}
go_bp_immune_terms = ['T cell activation involved in immune response (GO:0002286)',
 'innate immune response (GO:0045087)',
 'innate immune response activating cell surface receptor signaling pathway (GO:0002220)',
 'lymphocyte activation involved in immune response (GO:0002285)',
 'mucosal immune response (GO:0002385)',
 'natural killer cell activation involved in immune response (GO:0002323)',
 'negative regulation of immune response (GO:0050777)',
 'negative regulation of innate immune response (GO:0045824)',
 'neutrophil activation involved in immune response (GO:0002283)',
 'positive regulation of cytokine production involved in immune response (GO:0002720)',
 'positive regulation of immune effector process (GO:0002699)',
 'positive regulation of immune response (GO:0050778)',
 'positive regulation of production of molecular mediator of immune response (GO:0002702)',
 'regulation of adaptive immune response (GO:0002819)',
 'regulation of cytokine production involved in immune response (GO:0002718)',
 'regulation of humoral immune response (GO:0002920)',
 'regulation of immune effector process (GO:0002697)',
 'regulation of immune response (GO:0050776)',
 'regulation of innate immune response (GO:0045088)']
for term in go_bp_immune_terms:
    go_bp_immune[term] = go_bp_2021[term]

go_bp_cancer = {}
go_cancer = ['cellular response to tumor necrosis factor (GO:0071356)',
 'negative regulation of tumor necrosis factor production (GO:0032720)',
 'negative regulation of tumor necrosis factor superfamily cytokine production (GO:1903556)',
 'negative regulation of tumor necrosis factor-mediated signaling pathway (GO:0010804)',
 'positive regulation of tumor necrosis factor production (GO:0032760)',
 'positive regulation of tumor necrosis factor superfamily cytokine production (GO:1903557)',
 'regulation of tumor necrosis factor production (GO:0032680)',
 'regulation of tumor necrosis factor-mediated signaling pathway (GO:0010803)',
 'response to tumor necrosis factor (GO:0034612)',
 'tumor necrosis factor-mediated signaling pathway (GO:0033209)']
for term in go_cancer:
    go_bp_cancer[term] = go_bp_2021[term]

if genes == 'kegg_immune':
    geneset = kegg_immune
elif genes == 'kegg_cancer':
    geneset = kegg_cancer
elif genes == 'go_bp_immune':
    geneset = go_bp_immune
elif genes == 'go_bp_cancer':
    geneset = go_bp_cancer   
    
print('Geneset Prepared')

## 03. run ssGSEA (rough mode, permutation mode)

## 03-1. Rough mode
# ss = gp.ssgsea(data=input,
#                gene_sets=geneset, #go_bp_cancer,
#                outdir=None,
#                sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
#                no_plot=True)

#ss.res2d.sort_values('Name')
# nes = ss.res2d.pivot(index='Term', columns='Name', values='NES')
# nes.head(10)

## 03-2. Permutation mode
ss_permut = gp.ssgsea(data=input,
               gene_sets=geneset,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' for your custom metric
               permutation_num=20, # set permutation_num > 0, it will act like prerank tool
               no_plot=True, # skip plotting, because you don't need these figures
               processes=4, seed=9)

ss_permut_res = ss_permut.res2d.sort_values('FWER p-val')
ss_permut_res.to_csv(genes+"_result.txt", sep='\t')

#ss_permut_fwer_filt = ss_permut.res2d[ss_permut.res2d['FWER p-val'] < 0.05]

nes_permut = ss_permut.res2d.pivot(index='Term', columns='Name', values='NES')
nes_permut.to_csv(genes+"_NES.txt", sep='\t')

#from gseapy.plot import heatmap
#heatmap(nes_permut_gtex)

end = datetime.now()
print("Running time: {}".format(str(end-start)[:10]))
