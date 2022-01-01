#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
from collections import defaultdict

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Converting traitar data to standardized table of traits'
epi = """DESCRIPTION:
Formatting traitar table with at least the following columns:
[sample, phenotype, prediction]
... to the following:
[genome, domain, phylum, class, order, family, genus, species, trait1, ..., traitN]

A metadata file of all genomes is used to get the genome taxonomy data.

Output written as tsv file to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('traitar_output', metavar='traitar_output', type=str,
                    help='Traitar output file')
parser.add_argument('genome_metadata', metavar='genome_metadata', type=str,
                    help='Tab-delim metadata file that contains at least 2 columns: "accession" & "gtdb_taxonomy')
parser.add_argument('--version', action='version', version='0.0.1')


def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def parse_meta(infile):
    """
    Return: {genome_accession : taxonomy}
    """
    if infile is None:
        return None
    logging.info('Loading file: {}'.format(infile))
    regex = re.compile(r'[^a-zA-Z0-9_-]+')
    header = {}
    meta = {}
    msg = 'Cannot find column "{}"'
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            line = _decode(line).rstrip().split('\t')
            if i == 0:
                header = {x.lower():ii for ii,x in enumerate(line)}
                continue
            try:
                accession = regex.sub('_', str(line[header['accession']]))
            except KeyError:
                raise KeyError(msg.format('accession'))
            try:
                taxonomy = line[header['gtdb_taxonomy']].split(';')
            except KeyError:
                raise KeyError(msg.format('gtdb_taxonomy'))                
            meta[accession] = taxonomy
    logging.info('  No. of accessions: {}'.format(len(meta.keys())))
    return meta

def parse_traitar(infile, model='phypat+PGL'):
    """
    Parsing gene annotation file.
    Return: {genome : {trait : score}}
    """
    logging.info('Loading file: {}'.format(infile))
    trt = defaultdict(dict)
    header = {}
    status = {'records' : 0}
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            # line parse
            line = _decode(line).rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            if len(line) < len(header.keys()):
                msg = 'line {}: less columns than header; skipping!'
                logging.warning(msg.format(i+1))
                continue
            # genome
            try:
                genome_name = line[header['genome']]
            except KeyError:
                msg = 'Cannot find "{}" column in "{}"'
                raise KeyError(msg.format('genome', infile))
            # phenotype model
            try:
                phen_model = str(line[header['phenotype_model']])
            except KeyError:
                msg = 'Cannot find "{}" column in "{}"'
                raise KeyError(msg.format('phenotype_model', infile))
            except IndexError:
                msg = 'No model listed in line: {}'
                raise IndexError(msg.format(i+1))
            if phen_model != model:
                continue
            # phenotype
            try:
                trt_name = str(line[header['phenotype']])
            except KeyError:
                msg = 'Cannot find "{}" column in "{}"'
                raise KeyError(msg.format('phenotype', infile))
            # phenotype score
            try:
                trt_score = line[header['prediction_score']]
            except KeyError:
                msg = 'Cannot find "{}" column in "{}"'
                raise KeyError(msg.format('prediction_score', infile))
            # adding info            
            trt[genome_name][trt_name] = trt_score
            status['records'] += 1
    # status
    logging.info('  No. of records: {}'.format(status['records']))
    return trt

def get_all_trt(trt):
    """
    All phenotype names
    """
    all_trt = set()
    for genome in trt.keys():
        for trt_name in trt[genome].keys():
            all_trt.add(trt_name)
    logging.info('  No. of trait columns: {}'.format(len(all_trt)))
    return sorted(all_trt)

def write_trait_table(trt, meta):
    """
    Writing table of annotations
    """
    logging.info('Writing table to STDOUT...')
    # all annotations
    all_trt = get_all_trt(trt)
    header = ['genome', 'domain', 'phylum', 'class', 'order', 'family',
              'genus', 'species']
    print('\t'.join(header + all_trt))
    status = {'records' : 0}
    for genome in meta.keys():
        # genome taxonomy
        try:
            taxonomy = meta[genome]
        except KeyError:
            msg = 'Cannot find "{}" in metadata'
            raise KeyError(msg.format(genome))
        # counts
        trt_cnts = []
        for trt_name in all_trt:
            try:
                x = trt[genome][trt_name]
            except KeyError:
                x = 0
            trt_cnts.append(x)
        # writing line
        trt_cnts = [str(x) for x in trt_cnts]
        print('\t'.join([genome] + taxonomy + trt_cnts))
        status['records'] += 1
    logging.info('  No. of records written: {}'.format(status['records']))


def main(args):
    # parsing genome metadata
    meta = parse_meta(args.genome_metadata)
    # loading traitar data
    trt = parse_traitar(args.traitar_output)
    # writing table
    write_trait_table(trt, meta)    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
