#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from itertools import cycle, islice
import random 
from pylab import *

df = pd.read_csv('your_vcf_file_here.vcf', sep="\t")
df.columns = ["GENE", 'POS', 'REF', 'ALT', 'DP', 'AF', "INFO"]
df['INFO'] = df['INFO'].apply(str)
df['GENE'] = df['GENE'].apply(str)

#aa change
df_aa = df["INFO"]
df_aa = pd.DataFrame(df_aa)
df_sub = df_aa.INFO.str.split("|", expand=True)

df_aa_subset = df_sub[10]
df_aa_subset = pd.DataFrame(df_aa_subset)


#first column
df_gene = df["GENE"]
df_gene = pd.DataFrame(df_gene)
df_sub_gene = df_gene.GENE.str.split("_", expand=True)

df_gene_subset = df_sub_gene[1]
df_gene_subset = pd.DataFrame(df_gene_subset)


df['AA_change'] = df_aa_subset
df['Gene_Segment'] = df_gene_subset
df = df[['Gene_Segment', 'POS', 'REF', 'ALT', 'DP', 'AF', 'AA_change']]

#cleaning the protein column
df['AA_change'] = df['AA_change'].apply(str)
df_protein = df["AA_change"]
df_protein = pd.DataFrame(df_protein)
df_protein = df_protein.AA_change.str.split(".", expand=True)
df_protein = df_protein[1]
df_protein = pd.DataFrame(df_protein)



df['SNV'] = df_protein
df = df[['Gene_Segment', 'POS', 'REF', 'ALT', 'DP', 'AF', 'SNV']]


df.to_csv('clean_vcf.csv', header=False, index=None)


#Up to this point, we will have clean csv files with the information we need for making our plots.

###################################################################################################################
###################################################################################################################
####################### FILTERING BY FREQUENCY ####################################################################
###################################################################################################################
###################################################################################################################



Nhe = pd.read_csv('variant_frequencies_all.csv', sep=',')
Nhe = Nhe[~Nhe.Reference.str.contains("N")] # drop "SNPs" called in the barcode region
Nhe = Nhe[['Frequency', 'Gene', 'AA_change', 'Position', 'Reference', 'Change']] # drop columns we don't need
Nhe = Nhe[Nhe['Frequency'] > 3]  # drop any SNPs < 3% frequency  / Change to 0.03 if your frequencies have not been converted to percents

print("\n SNV's with frequencies above 1% are: ")
Nhe.to_csv('variant_frequencies_above_3percent.csv', sep='\t', header=True, index=None)

#end
