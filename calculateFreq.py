#!/usr/bin/python
"Calculate frequency of mutations withing inverted repeats stratified by various parameters"
import os
import pandas as pd

HOME = '/home/fa_bula/PCAWG/'
REPEATS_DIR = os.path.join(HOME, 'nonBDNA/inverted_repeats/gff/')
RESULTS_DIR = os.path.join(HOME, 'results')
MUTATIONS_DIR = os.path.join(HOME, 'mutations_PCAWG')
cancer_types  = ['BLCA', 'BRCA', 'CESC', 'HNSC', 'LUAD', 'LUSC']

byvars_dict = {'targets': ['loc']}
for c in cancer_types:
    byvars_dict[c] = byvars_dict['targets'] + [c + '_sample']

mutation_files = [os.path.join(RESULTS_DIR, 'chr' + str(i) + '.csv') for i in range(1,23)]
df_list = []

dtypes = {'chr':'string', 'pos':'int64'}
for c in cancer_types:
    dtypes[c + '_mutation'] = 'string'
    dtypes[c + '_sample'] = 'string'

enrichment = pd.read_csv(os.path.join(HOME, "PCAWG_enrichment.txt"), sep="\t")
enrichment["Sample"] = enrichment["Sample"].apply(lambda x: x.strip()) # Remove all leading/triling spaces from Sample name

for f in mutation_files:
    df0 = pd.read_csv(f, sep='\t', dtype=dtypes)
    df1 = df0.loc[(((df0['ref'] == 'C') & (df0['rel_pos_inv'] == 1)) | ((df0['ref'] == 'G') & (df0['rel_pos'] == 1)))]
    df1 = df0.loc[((df0['ref'] == 'C')| (df0['ref'] == 'G'))]
    del df0
    df = df1.loc[df1['loc'] == 'L']
    del df1
    df_list += [df]
    del df

df_all = pd.DataFrame().append(df_list)
del df_list

targets_freq = df_all.drop_duplicates(subset=['chr', 'pos'] + byvars_dict['targets']).groupby(byvars_dict['targets']).size().to_frame(name="targets").reset_index()

for key in byvars_dict:
    if key != 'targets':
        freq0 = df_all.drop_duplicates(subset=['chr', 'pos'] + byvars_dict[key]).groupby(byvars_dict[key]).size().to_frame(name="mutations").reset_index()
        freq = pd.merge(freq0, targets_freq, how='left', on=byvars_dict['targets'])
        freq = freq.merge(enrichment[['Sample', 'APOBEC_enrichment']], how='left', left_on=[key + '_sample'], right_on=['Sample'])
        freq = freq[byvars_dict[key] + ['mutations', 'targets', "APOBEC_enrichment"]]
        freq.to_csv(os.path.join(RESULTS_DIR, key + '_freq.csv'), sep='\t', header=True, index=False)
