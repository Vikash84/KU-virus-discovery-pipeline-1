import argparse
import pandas as pd
import os

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for filtering blast result by taxonomy')

parser.add_argument('--input', '-i', metavar='blast_result.txt',
        help='Input blast result file',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='reference_summary_filtered.txt',
        help='Output file name where filtered result will be written')

parser.add_argument('--level', metavar='e.g.Superkingdom',
        help='Taxonmy level to be filtered be filtered')

parser.add_argument('--name', metavar='e.g.Viruses',
        help='Taxonmy name to be filtered be filtered')

args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=0)

taxon_filter = df[args.level] == args.name
filtered_df = df[taxon_filter]
filtered_df.to_csv(args.output, index=False, sep='\t')

