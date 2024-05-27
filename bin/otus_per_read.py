#!/usr/bin/env python3
import gzip
import argparse
import csv
import logging
import math

import pandas as pd
import numpy as np
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord
from enum import Enum
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='')
    # def file_choices(choices,fname):
    #     if any([fname.endswith(ext) for ext in choices]):
    #         return fname
    #     else:
    #         parser.error("file doesn't end with one of {}".format(choices))
        
    # parser.add_argument('-c', '--cluster', required = True)
    # parser.add_argument('-s', '--sample', required = True)
                        # type=lambda s:file_choices(("fastq.gz","fq.gz", "fq", "fastq"), s),
                        # help='uc file which we will extract regions from')
    parser.add_argument('-o', '--output', nargs='?')
    parser.add_argument('file', nargs='+')
    args = parser.parse_args()

    print(args.file)

    file_by_region = []
    for file in args.file:
        region = Path(file).stem
        # file_by_region[region] = pd.read_csv(file, sep='\t')
        df = pd.read_csv(file, sep='\t')
        file_by_region.append((region, df))

    (region, merged_df) = file_by_region[0]
    # print("first region:", region)
    merged_df = merged_df.add_suffix(f'_{region}')

    # print([r for (r,_) in file_by_region])
    # print([r for (r,_) in file_by_region[1:]])

    

    for (next_region, next_df) in file_by_region[1:]:
        next_df = next_df.add_suffix(f'_{next_region}')
        merged_df = pd.merge(
            merged_df,
            next_df,
            left_on=f'Column.QueryId_sample_{region}',
            right_on=f'Column.QueryId_sample_{next_region}',
            suffixes=['', f'_{next_region}'],
            how='outer',
            validate='1:1'
        )

        def use_non_na(row):
            first = row[f'Column.QueryId_sample_{region}']
            if isinstance(first, float) and math.isnan(first):
                return row[f'Column.QueryId_sample_{next_region}'] 
            else:
                return first

        merged_df[f'Column.QueryId_sample_{region}'] = merged_df[[f'Column.QueryId_sample_{region}', f'Column.QueryId_sample_{next_region}']].apply(
            use_non_na,
            axis = 'columns'
        )
    
    # print(merged_df.columns)
    cols = [f'Column.QueryId_sample_{region}'] + [f'OTU_{r}' for (r,_) in file_by_region]
    # print(cols)

    result = merged_df[cols]
    result.columns = ['ReadId'] + [r for (r,_) in file_by_region]
    print(result)
    
    if args.output:
        result.to_csv(args.output, sep='\t', index=False, na_rep='')


def open_file_based_on_ext(fpath):
    if fpath.endswith('.gz'):
        return gzip.open(fpath, "rt")
    else:
        return open(fpath)

main()


