import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

if __name__ == '__main__':
    dvf_table = pd.read_csv(sys.argv[1], sep='\t', dtype={'pvalue': np.float64}, header=0)
    is_significant = dvf_table['pvalue'] < 0.05
    significant_rows = dvf_table.loc[is_significant]
    
    sorted_rows = significant_rows.sort_values(by=['pvalue'])
    significant_ids = list()
    for index, row in sorted_rows.iterrows():
        split_row = row[0].split()
        significant_ids.append(split_row[0])

    filtered_contigs = open(sys.argv[3], "w")

    with open(sys.argv[2]) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in significant_ids:
                SeqIO.write(record, filtered_contigs, "fasta")
