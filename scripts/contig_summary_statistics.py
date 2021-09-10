import sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    prefix = sys.argv[1]
    record = list(SeqIO.parse(sys.argv[2], "fasta"))
    num_seq = len(record)
    sizes = [len(rec) for rec in record]

    # calculate a 5-number summary
    quartiles = np.percentile(sizes, [25, 50, 75])
    data_min = min(sizes)
    data_max = max(sizes)
    data_q1 = int(quartiles[0])
    data_median = int(quartiles[1])
    data_q3 = int(quartiles[2])

    # write summary to output file
    output_file = open(prefix + ".contig.summary.txt", "w")
    output_file.write("Number of assembled contigs: " + str(num_seq)+'\n')
    output_file.write("Min length: " + str(data_min) + "bp"+'\n')
    output_file.write("Q1 length: " + str(data_q1) + "bp"+'\n')
    output_file.write("Median length: " + str(data_median) + "bp"+'\n')
    output_file.write("Q3 length: " + str(data_q3) + "bp"+'\n')
    output_file.write("Max length: " + str(data_max) + "bp"+'\n')
    output_file.close()
    
    # export contig length histogram
    plt.hist(sizes, bins='auto')
    plt.title(
    "%i Assembled contigs\nLengths %i to %i" % (num_seq, data_min, data_max)
)
    plt.xlabel("Sequence length (bp)")
    plt.ylabel("Count")
    plt.savefig(prefix + ".contigs.length_histogram.png")
