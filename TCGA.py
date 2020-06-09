# Load python packages
import pandas as pd
import numpy as np
import xlrd
import matplotlib.pyplot as plt

# ENV settings
datadir = './data/'

# Load L-Score KRAS Signature Genes
l_score = pd.read_csv(datadir+'l_score_genes.txt', header=None, sep='\t')
l_score = l_score[l_score.columns[:-1]]
l_score.columns = ['ENSG', 'EntrezID', 'Symbol', 'Direction']
l_up = l_score[l_score['Direction']=='up']
l_dn = l_score[l_score['Direction']=='dn']

l_up.insert(1, 'ENSG_pref', [ensg.split('.')[0] for ensg in l_up['ENSG']])
l_dn.insert(1, 'ENSG_pref', [ensg.split('.')[0] for ensg in l_dn['ENSG']])

# Load S-Score KRAS Signature Genes
s_score = pd.read_csv(datadir+'s_score_genes.txt', header=None, sep='\t')
s_score.columns = ['ENSG', 'EntrezID', 'Symbol', 'Direction']
s_up = s_score[s_score['Direction']=='UP']
s_dn = s_score[s_score['Direction']=='DN']

s_up.insert(1, 'ENSG_pref', [ensg.split('.')[0] for ensg in s_up['ENSG']])
s_dn.insert(1, 'ENSG_pref', [ensg.split('.')[0] for ensg in s_dn['ENSG']])



# Load TCGA PAAD Sample annotation
tcga_paad_annotation = pd.read_excel(datadir+'tcga.paad.sample_data.xlsx', sheet_name=None)
tcga_paad_annotation_table = tcga_paad_annotation['FreezeSamples'].T.reset_index(drop=True).set_index(0).T

# Select only "high-purity" samples from TCGA PAAD cohort
tcga_paad_purity_table = tcga_paad_annotation_table[['Tumor Sample ID', 'ABSOLUTE Purity']]

# Determine whether sample is high or low purity
tcga_paad_purity_table.insert(tcga_paad_purity_table.shape[1], 'Purity Group', 
                              ['High' if tcga_paad_purity_table.loc[i]['ABSOLUTE Purity'] >= 0.33 
                               else 'Low' if tcga_paad_purity_table.loc[i]['ABSOLUTE Purity'] < 0.33 else None 
                               for i in tcga_paad_purity_table.index])



# Load TCGA-PAAD expression profiles
tcga_paad_exp = pd.read_csv(datadir+'tcga.paad.hiseqv2', sep='\t', index_col=0)

# Keep only high-purity samples with expression data from UCSC Xena
tcga_paad_exp_high_purity = tcga_paad_exp[[pat_id[:-1] for pat_id in list(tcga_paad_purity_table[tcga_paad_purity_table['Purity Group']=='High']['Tumor Sample ID'])]]

# Get raw TPM expression values from TCGA data (don't remove TPM pseudocount of 0.001)
tcga_paad_tpm_high_purity = 2**tcga_paad_exp_high_purity

# log10 transform TPM matrix
tcga_paad_log10tpm_high_purity = np.log10(tcga_paad_tpm_high_purity)

# Calculate average log10(TPM) for each gene
tcga_paad_log10tpm_high_purity_mean = tcga_paad_log10tpm_high_purity.mean(axis=1)

# Calculate gene expression relative to mean for every gene (subtract sample log10(TPM) from average log10(TPM))
tcga_paad_log10tpm_high_purity_ratio = tcga_paad_log10tpm_high_purity.subtract(tcga_paad_log10tpm_high_purity_mean, axis=0)

# Save raw expression (TPM) values
df = pd.DataFrame(data=tcga_paad_log10tpm_high_purity_ratio).T
df.to_csv('TCGA_Raw_TPMs.csv')


# Membership
print('Genes in L-Score Up gene set not found in current expression table:', set(l_up['Symbol'])-set(tcga_paad_log10tpm_high_purity_ratio.index))
print('Genes in L-Score Down gene set not found in current expression table:', set(l_dn['Symbol'])-set(tcga_paad_log10tpm_high_purity_ratio.index))
print('Genes in S-Score Up gene set not found in current expression table:', set(s_up['Symbol'])-set(tcga_paad_log10tpm_high_purity_ratio.index))
print('Genes in S-Score Down gene set not found in current expression table:', set(s_dn['Symbol'])-set(tcga_paad_log10tpm_high_purity_ratio.index))

# L-Score up-genes: {'CXCL8', 'MYDGF'}
l_up = l_up.replace({'CXCL8':'IL8', 'MYDGF':'C19orf10'})

# L-Score down-genes: {'CCSAP', 'CEP57L1', 'SUGP1'}
l_dn = l_dn.replace({'CCSAP':'C1orf96', 'CEP57L1':'C6orf182', 'SUGP1':'SF4'})

# S-Score down-genes: {'TMEM237'}
s_dn = s_dn.replace({'TMEM237':'ALS2CR4'})



# L-Score
tcga_paad_high_purity_lscore = tcga_paad_log10tpm_high_purity_ratio.loc[l_up['Symbol']].mean().subtract(tcga_paad_log10tpm_high_purity_ratio.loc[l_dn['Symbol']].mean())

lsdf = pd.DataFrame(data=tcga_paad_high_purity_lscore).T
lsdf.to_csv('LScores.csv')

# Plot bargraph of L-Score for each sample (samples sorted by L-Score)
tcga_paad_high_purity_lscore.sort_values().plot.bar(width=0.8, figsize=(16, 6))
plt.ylabel('L-Score', fontsize=20)
plt.tick_params(axis='y', labelsize=16)



# S-Score
tcga_paad_high_purity_sscore = tcga_paad_log10tpm_high_purity_ratio.loc[s_up['Symbol']].mean().subtract(tcga_paad_log10tpm_high_purity_ratio.loc[s_dn['Symbol']].mean())

ssdf = pd.DataFrame(data=tcga_paad_high_purity_sscore).T
ssdf.to_csv('SScores.csv')

# Plot bar-chart of L-Score for each sample (samples sorted by L-Score)
tcga_paad_high_purity_sscore.sort_values().plot.bar(width=0.8, figsize=(16, 6))
plt.ylabel('S-Score', fontsize=20)
plt.tick_params(axis='y', labelsize=16)

# Plot scatter plot comparing L-Score (x-axis) to S-Score (y-axis)
score_align = pd.concat([tcga_paad_high_purity_lscore.rename('L-Score'), 
                         tcga_paad_high_purity_sscore.rename('S-Score')],
                        axis=1, sort=False)
plt.figure(figsize=(6,6))
plt.scatter(score_align['L-Score'], score_align['S-Score'], s=50)
plt.xlabel('L-Score', fontsize=20)
plt.ylabel('S-Score', fontsize=20)
plt.tick_params(axis='both', labelsize=16)

