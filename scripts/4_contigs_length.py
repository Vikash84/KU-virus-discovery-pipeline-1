import sys
from Bio import SeqIO

if __name__ == '__main__':
    prefix = sys.argv[1]
    record = list(SeqIO.parse(sys.argv[2], "fasta"))
    sizes = [len(rec) for rec in record]

    # write each sequence's length to output file
    with open(prefix + ".contigs.len", "w") as output_file:
        for rec in record:
            output_file.write(str(len(rec)) + '\n')