import argparse
import pandas as pd
import os

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for filtering reference mapping result')

parser.add_argument('--input', '-i', metavar='reference_summary.txt',
        help='Input reference mapping result file',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='reference_summary_filtered.txt',
        help='Output file name where filtered result will be written')

parser.add_argument('--min_avg_cov', type=float, default=1.0,
        help='Mapping below this average coverage will be filtered')

args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t', header=0, usecols=lambda x: x != "START")
df = df.rename(columns={"END" : "LEN"})

min_avg_cov_filter = df['AVG_COV'] > args.min_avg_cov
composite_filter = min_avg_cov_filter
filtered_df = df[composite_filter]
sorted_filtered_df = filtered_df.sort_values(by='PERCENT_COVERED', ascending=False)
sorted_filtered_df.to_csv(args.output, index=False, sep='\t')
