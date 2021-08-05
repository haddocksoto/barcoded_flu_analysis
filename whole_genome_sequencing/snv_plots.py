#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches


# Import your vcf files

# here, Nhe and Pst refers to the registration amrks that were added to the HA gene segment. 
#Each sample was prepped to sequence all 8 gene segments, but only a few of those had variatns exceeding out 3% frequency treshold


Nhe = pd.read_csv('Nhe_2019_wgs/variant_frequencies_above_3percent.csv', sep='\t')
Pst = pd.read_csv('Pst_2019_wgs/variant_frequencies_above_3percent.csv', sep='\t')

#selecting data for each gene segment in each stock (here, only the ones that showed variants above 3%):

Nhe_HA = Nhe[Nhe['Gene'] == 'HA']
Nhe_PB1 = Nhe[Nhe['Gene'] == 'PB1']
Nhe_PB2 = Nhe[Nhe['Gene'] == 'PB2']


Pst_HA = Pst[Pst['Gene'] == 'HA']
Pst_PB1 = Pst[Pst['Gene'] == 'PB1']
Pst_PB2 = Pst[Pst['Gene'] == 'PB2']


#this code could be simplified, but that will be done in the future

# plot frequency of the SNP against nucleotide position 

figure(figsize=(8,4), dpi=300)

mpl.rcParams['axes.spines.left'] = True
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = True

axes = plt.gca()
axes.set_xlim([1, 13501]) #our IAV genome has barcode, registration marks and Nluc additions. That's why is larger that wild type IAV
axes.set_ylim([0,100])

PB1 = [0, 2280]

#adding box for each gene segment
axes.axvspan(0, 2280, facecolor='#C5C5C9', alpha=0.3, zorder=1)
axes.axvspan(4554, 7354, facecolor='#C5C5C9', alpha=0.3, zorder=1)
axes.axvspan(6900, 7354, facecolor='#3BC6C8', alpha=0.3, zorder=1) # this one is for the NLuc insertion
axes.axvspan(9157, 10672, facecolor='#C5C5C9', alpha=0.3, zorder=1)
axes.axvspan(12082, 12841, facecolor='#C5C5C9', alpha=0.3, zorder=1)

#Adding labels for each gene segment (and NLuc)

plt.text(900,103,'PB2', fontweight='bold', size=12)
plt.text(3050,103,'PB1', fontweight='bold', size=12)
plt.text(6800,102,'NLuc', fontweight='bold', size=8)
plt.text(5700,103,'PA', fontweight='bold', size=12)
plt.text(8050,103,'HA', fontweight='bold', size=12)
plt.text(9700,103,'NP', fontweight='bold', size=12)
plt.text(11150,103,'NA', fontweight='bold', size=12)
plt.text(12300,103,'M', fontweight='bold', size=12)
plt.text(13000,103,'NS', fontweight='bold', size=12)


#Assign colors to each data set in the legend
Nhe_patch = mpatches.Patch(color='blue', label='NheI')
Pst_patch = mpatches.Patch(color='gold', label='PstI')
plt.legend(handles=[Nhe_patch, Pst_patch], loc='center left')

#add albels
plt.ylabel("SNV Frequency (%)")
plt.xlabel("\n IAV Gene Position \n \n ")

#plotting the data
plt.scatter(Nhe_HA.Position+7354, Nhe_HA.Frequency, alpha=0.7, facecolors='blue', edgecolors='blue', zorder=10)
plt.scatter(Nhe_PB1.Position, Nhe_PB1.Frequency, alpha=0.7, facecolors='blue', edgecolors='blue', zorder=10)
plt.scatter(Nhe_PB2.Position+2280, Nhe_PB2.Frequency, alpha=0.7, facecolors='blue', edgecolors='blue', zorder=10)
plt.scatter(Pst_HA.Position+7354, Pst_HA.Frequency, alpha=0.7, facecolors='#FAD401', edgecolors='#FAD401', zorder=10)
plt.scatter(Pst_PB1.Position, Pst_PB1.Frequency, alpha=0.7, facecolors='#FAD401', edgecolors='#FAD401', zorder=10)
plt.scatter(Pst_PB2.Position+2280, Pst_PB2.Frequency, alpha=0.7, facecolors='#FAD401', edgecolors='#FAD401', zorder=10)

#we labeled some variants above a frequency of 10%. here is the code we used for that

plt.text(508+7354,49.6164,'  K153E' ,horizontalalignment='left', color='#B33040', size=8, fontweight='bold')
plt.text(508+7354,68.7833,'  K154E' ,horizontalalignment='left', color='#B33040', size=8, fontweight='bold')
plt.text(515+7354,15.0136,'  G155E' ,horizontalalignment='left', color='black', size=8, fontweight='bold')

#added a nice title
plt.title('\n YOUR_TITLE_HERE \n\n\n')

#we set up a 10% line to distinguish tose non-synonymopus variants that exceeded 10%
plt.axhline(y = 10, color = 'black', linestyle = 'dashed', alpha=0.1)

#saving the figures
plt.savefig('SNV_plot_skewed.png', dpi=300, bbox_inches='tight')
plt.savefig('SNV_plot_skewed.pdf', dpi=300, bbox_inches='tight') 
