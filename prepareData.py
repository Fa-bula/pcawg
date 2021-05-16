#!/usr/bin/python
"Prepare all information about inverted repeats positions in a single file"
import os
import pandas as pd
from Bio import SeqIO


HOME = '/home/fa_bula/PCAWG/'
REPEATS_DIR = os.path.join(HOME, 'nonBDNA/inverted_repeats/gff/')
RESULTS_DIR = os.path.join(HOME, 'results')
GENOME_DIR = os.path.join(HOME, 'genome')
MUTATIONS_DIR = os.path.join(HOME, 'mutations_PCAWG')
CHR =  os.environ['SGE_TASK_ID']
if CHR == '23':
    CHR = 'X'

if CHR == '24':
    CHR = 'Y'

cancer_types  = ['BLCA', 'BRCA', 'CESC', 'HNSC', 'LUAD', 'LUSC']

# Reading genome
genome_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(GENOME_DIR, "hg19.fa"), "fasta"))
genome_seq = genome_dict['chr' + CHR].seq.upper()

# Reading file with inverted repeats
repeats_file = os.path.join(REPEATS_DIR, 'chr' + CHR + '_IR.gff')
repeats = pd.read_csv(repeats_file, sep='\t|;', comment='#', engine='python',
                 names=["chr", "ABCC", "type", "begin", "end", "foo", "bar", "baz", "id", "spacer",
                        "repeat", "perms", "minloop", "composition", "sequence", "subset"],
                 usecols=["chr", "begin", "end", "spacer", "repeat", "sequence", "id", "perms"])

repeats['chr'] = repeats['chr'].replace('chr', '', regex=True)
repeats['id'] = repeats['id'].replace('ID=', '', regex=True)
repeats['spacer'] = repeats['spacer'].replace('spacer=', '', regex=True)
repeats['repeat'] = repeats['repeat'].replace('repeat=', '', regex=True)
repeats['sequence'] = repeats['sequence'].replace('sequence=', '', regex=True)
repeats = repeats.loc[repeats['perms'] == "perms=1"]

outdata = []                       # List to be converted to Data Frame
for index, row in repeats.iterrows():
    for i,c in enumerate(row['sequence']):
        d = {}
        spacer = int(row['spacer'])
        repeat = int(row['repeat'])
        d['chr'] = row['chr']
        d['id'] = row['id']
        d['pos'] = row['begin'] + i
        d['spacer'] = spacer
        d['repeat'] = repeat
        d['ref'] = c.upper()
        # Previous nucleotide
        if i == 0:
            d['prev'] = genome_seq[d['pos'] - 1 - 1] # -1 because we want previous and -1 because characters in genome_seq enumerated from 0
        else:
            d['prev'] = row['sequence'][i-1].upper()
        # Next Nucleotide
        if i + 1 == len(row['sequence']):
            d['next'] = genome_seq[d['pos'] - 1 + 1] # +1 because we want next and -1 because characters in genome_seq enumerated from 0  
        else:
            d['next'] = row['sequence'][i+1].upper()
        # Position location inside inverted repeat
        if i < repeat:
            d['loc'] = 'S1'     # Stem 1
            d['rel_pos'] = i + 1 # Relative position (1 is first) in Stem 1
            d['rel_pos_inv'] = repeat - i # Relative inverted position (1 is last) in Stem 1
        else:
            if i < repeat + spacer:
                d['loc'] = 'L'  # Loop
                d['rel_pos'] = i + 1 - repeat # Relative position in Loop
                d['rel_pos_inv'] = repeat + spacer - i # Relative inverted position (1 is last) in Loop
            else:
                d['loc'] = 'S2' # Stem 2
                d['rel_pos'] = i + 1 - repeat - spacer # Relative position in Stem 2
                d['rel_pos_inv'] = 2 * repeat + spacer - i # Relative inverted position (1 is last) in Stem 2
            
        outdata.append(d)

# Data frame with all position in inverted repeats
out = pd.DataFrame(outdata)[['chr', 'pos', 'prev', 'ref', 'next', 'loc', 'rel_pos', 'rel_pos_inv', 'spacer', 'repeat', 'id']]
del outdata, repeats

# Merge mutation informations into out Data Frame
for cancer_type in cancer_types:
    m = pd.read_csv(os.path.join(MUTATIONS_DIR, "snv_mnv_" + cancer_type + "-US.tsv"),
                    sep='\t', low_memory=False, usecols=['CHROM', 'POS', 'ALT', 'sample'])
    m.columns = ['chr', 'pos', cancer_type + '_mutation', cancer_type + '_sample']
    # Filter mutations by current chromosome
    mf = m.loc[m['chr'] == CHR]
    out = out.merge(mf, how='left', on = ['chr', 'pos'])

out.to_csv(os.path.join(RESULTS_DIR, 'chr' + CHR + '.csv'), sep='\t', index=False)
