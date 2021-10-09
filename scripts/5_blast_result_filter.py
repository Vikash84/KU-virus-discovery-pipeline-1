import argparse
import pandas as pd
import os

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for filtering Blast result in format 6')

parser.add_argument('--input', '-i', metavar='blast.txt',
        help='Input blast result file',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='blast_filtered.txt',
        help='Onput file name where filtered result will be written')
    
parser.add_argument('--exclude', '-e', metavar='excluding_keyword.txt',
        help='List of keywords to be excluded from blast result',
        type=lambda x: is_valid_file(parser, x)) 

parser.add_argument('--aln_len', '-al', type=int, default=0,
        help='Hit below this align length will be filtered')

parser.add_argument('--aln_contig_len_prop', type=float, default=0,
        help='Hit below this align length / contig length will be filtered')

args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=0)

aln_len_filter = df['ALN_LEN'] > args.aln_len
aln_contig_len_prop_filter = ( df['ALN_LEN'] / df['QUERY_LEN'] ) > args.aln_contig_len_prop

composite_filter = aln_len_filter & aln_contig_len_prop_filter

if args.exclude:
    with open(args.exclude) as f:
        content = f.readlines()

        excluding_keywords = [x.strip() for x in content]
        exclude_filter = df['REF_TITLE'].str.contains('(?i)' + '|'.join(excluding_keywords)) == False

        composite_filter = composite_filter & exclude_filter

filtered_df = df[composite_filter]
evalue_sorted_df = filtered_df.sort_values(by='EVALUE')
sorted_filtered_df = evalue_sorted_df.sort_values(by='BITSCORE', ascending=False)

sorted_filtered_df.to_csv(args.output, index=False, sep='\t')
