#!/usr/bin/env python3
import gzip
import argparse
import csv
import logging

import pandas as pd
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord
from enum import Enum

class Column(Enum):
            RecordType = 'RecordType'
            ClusterNum = 'ClusterNum'
            SeqLength = 'SeqLength'
            PctIdentity = 'PctIdentity'
            Strand = 'Strand'
            NA_6 = 'NA_6'
            NA_7 = 'NA_7'
            Alignment = 'Alignment'
            QueryId = 'QueryId'
            ReferenceId = 'ReferenceId'

class AdditionalColumn(Enum):
            PlainQueryId = 'PlainQueryId'
            PlainRefId = 'PlainRefId'

#
# Run with:
# bin/recover_reads_from_otus.py -s <(cat output/full-multi-region/vsearch-derep/$region/*Sample*/*.uc) \
#           -p output/full-multi-region/vsearch-derep/$region/pooled_reads/pooled_reads.uc \
#           -c output/full-multi-region/vsearch-cluster/$region/pooled_reads/pooled_reads.uc.tsv.gz  \
#           -o scratch/otu-recovery/$region.tsv
def main():
    parser = argparse.ArgumentParser(description='')
    # def file_choices(choices,fname):
    #     if any([fname.endswith(ext) for ext in choices]):
    #         return fname
    #     else:
    #         parser.error("file doesn't end with one of {}".format(choices))
        
    parser.add_argument('-c', '--cluster', required = True)
    parser.add_argument('-s', '--sample', required = True)
                        # type=lambda s:file_choices(("fastq.gz","fq.gz", "fq", "fastq"), s),
                        # help='uc file which we will extract regions from')
    parser.add_argument('-p', '--pooled', required = True)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    sample_derep = read_uc_file(args.sample)
    pooled_derep = read_uc_file(args.pooled)
    cluster_df = read_uc_file(args.cluster)

    print("samples: ", args.sample)
    print("pooled: ", args.pooled)
    print("cluster: ", args.cluster)

    print("Pooled reads derep")

    test_id = 'ffffedeb-4ea3-43b1-928c-17dd973edbaf'
    
    
    #print(pooled_read_derep.iloc[1:5, :])

    print(pooled_derep.loc[(pooled_derep[AdditionalColumn.PlainQueryId] == test_id)
                                | (pooled_derep[AdditionalColumn.PlainRefId] == test_id)])

    print(cluster_df.loc[(cluster_df[AdditionalColumn.PlainQueryId] == test_id)])

    first_merge = pd.merge(sample_derep, pooled_derep,
                           left_on=AdditionalColumn.PlainRefId,
                           right_on=AdditionalColumn.PlainQueryId,
                           suffixes=['_sample', '_pooled'])
    
    print(first_merge)

    merged = pd.merge(first_merge, cluster_df.add_suffix('_cluster'), 
                      left_on=f'{AdditionalColumn.PlainRefId}_pooled', 
                      right_on=f'{AdditionalColumn.PlainQueryId}_cluster')

    print("Merged:")
    
    print(merged.columns)


    r = pd.DataFrame( 
         merged[
            [f'{Column.QueryId}_sample', f'{Column.QueryId}_pooled', f'{Column.QueryId}_cluster', f'{AdditionalColumn.PlainRefId}_cluster', f'{Column.ClusterNum}_cluster']
        ]
    )
    r.columns =[
        f'{Column.QueryId}_sample', f'{Column.QueryId}_pooled', f'{Column.QueryId}_cluster', f'{AdditionalColumn.PlainRefId}_cluster', 'OTU'
    ]
    r['OTU'] = r['OTU'].map(lambda otu: f'OTU_{otu+1}') 

    print(r)

    if args.output:
        r.to_csv(args.output, sep='\t')


def read_uc_file(path):
    df = pd.read_csv(path, sep='\t', header=None, names = [c for c in Column])
    df = fmt_uc_dataframe(df)
    df[AdditionalColumn.PlainQueryId] = df[Column.QueryId].map(strip_id)
    df[AdditionalColumn.PlainRefId] = df[Column.ReferenceId].map(strip_id)
    return df

def strip_id(full_id):
    return full_id.split(';', 1)[0]

def simplify_ref_id(row):
    return strip_id(replace_stars_ref_id(row))

def replace_stars_ref_id(row):
    if row[Column.ReferenceId] == '*':
        return row[Column.QueryId]
    return row[Column.ReferenceId]

def fmt_uc_dataframe(df):
    no_cs = df.loc[df[Column.RecordType] != 'C']
    no_cs.loc[:, [Column.ReferenceId]] = no_cs[[Column.QueryId, Column.ReferenceId]].apply(replace_stars_ref_id, axis="columns")

    return no_cs

def find_otus_for_query_seq(row, cluster_df):
    plain_query_id = row[Column.QueryId].split(';', 1)[0]
    plain_ref_id = row[Column.ReferenceId].split(';', 1)[0]

    id_in_cluster = plain_query_id if row[Column.RecordType] in ['C', 'S'] else plain_ref_id
    print(id_in_cluster)


def open_file_based_on_ext(fpath):
    if fpath.endswith('.gz'):
        return gzip.open(fpath, "rt")
    else:
        return open(fpath)

main()


