#!/usr/bin/env python3
import gzip
import argparse
import csv
import logging
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description='')
    def file_choices(choices,fname):
        if any([fname.endswith(ext) for ext in choices]):
            return fname
        else:
            parser.error("file doesn't end with one of {}".format(choices))
        
    parser.add_argument('-p', '--positions', required = True,
                        help='position (tsv) file output from ITSx')
    parser.add_argument('-f', '--fastq', required = True,
                        type=lambda s:file_choices(("fastq.gz","fq.gz", "fq", "fastq"), s),
                        help='fastq file which we will extract regions from')
    parser.add_argument('-o', '--outfile', required = True,
                        help='output fastq file')
    args = parser.parse_args()

    with open_file_based_on_ext(args.fastq) as handle:
        original_reads = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))

    with open(args.positions) as file:
        regions = extract_regions(original_reads, file)
        
        if args.outfile.endswith('.gz'):
            with bgzf.BgzfWriter(args.outfile, "wb") as outgz:
                SeqIO.write(regions, outgz, "fastq") 
        else:
            SeqIO.write(regions, args.outfile, "fastq") 

def open_file_based_on_ext(fpath):
    if fpath.endswith('.gz'):
        return gzip.open(fpath, "rt")
    else:
        return open(fpath)

def extract_regions(original_reads, position_file):
    positions = csv.reader(position_file, delimiter="\t")

    print(original_reads)
    regions = []
    for line in positions:
        id = line[0]
        lsu = line[6]

        if 'Not found' in lsu:
            logging.warning("LSU region was not found for position ('%s')" % id)
            continue

        [r_min_1, r_max] = [int(r) for r in lsu.replace('LSU: ', '').split('-')]
        r_min = r_min_1 - 1 # adjust for 1 based indexing

        print(id, r_min, r_max)

        if id not in original_reads:
            logging.warning("id ('%s') not found in fastq file" % id)
            continue
        else:
            o = original_reads[id]

            letter_annotations = slice_quality_scores(o, r_min, r_max)
            seq_region = o.seq[r_min:r_max]
            print(seq_region)

            if len(seq_region) == 0:
                logging.warning("region [%d:%d] not present for read with id ('%s')" % (r_min, r_max, id))
                continue

            subregion = SeqRecord(
                id = o.id, name = o.name, description = o.description, 
                annotations = o.annotations, features = o.features, dbxrefs = o.dbxrefs,
                seq = seq_region, letter_annotations = letter_annotations)
            
            yield subregion

def slice_quality_scores(seq_record, r_min, r_max):
    subscores = {}
    for key,scores in seq_record.letter_annotations.items():
        subscores[key] = scores[r_min:r_max]
    return subscores

main()


