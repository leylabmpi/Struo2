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

desc = 'Creating table of gene content per genome'
epi = """DESCRIPTION:
Create table of annoted genome content.
Example: number of genes annotated as each COG category, per genome.

The script uses genes annotated with UniRef IDs via the Struo2 pipeline.

If --genome-metadata is provided (genes per genome), then the annotation
acounts are normalized by total genes per genome. The output files will 
then be fractions (0-1) instead of integers.

Output is written to STDOUT; tab-delim file with a header.
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('annotations', metavar='annotations', type=str,
                    help='Struo2 HUMAnN3 gene annotation file')
parser.add_argument('--content', type=str,
                    choices = ['uniref', 'cog', 'cog-categories', 'pfam', 'pfam-categories'],
                    default='/ebio/abt3_projects/databases_no-backup/humann3/201901/utility_mapping/map_eggnog_uniref90.txt.gz',
                    help='UniRef <=> COG mapping file')
parser.add_argument('--cog-map', type=str,
                    default='/ebio/abt3_projects/databases_no-backup/humann3/201901/utility_mapping/map_eggnog_uniref90.txt.gz',
                    help='UniRef <=> COG mapping file')
parser.add_argument('--cog-cat', type=str,
                    default='/ebio/abt3_projects/databases_no-backup/humann3/201901/utility_mapping/hierarchy/COG.tsv',
                    help='COG <=> COG-category mapping file')
parser.add_argument('--pfam-map', type=str,
                    default='/ebio/abt3_projects/databases_no-backup/humann3/201901/utility_mapping/map_pfam_uniref90.txt.gz',
                    help='UniRef <=> Pfam mapping file')
parser.add_argument('--pfam-cat', type=str,
                    default='/ebio/abt3_projects/databases_no-backup/humann3/201901/utility_mapping/hierarchy/pfam.tsv',
                    help='Pfam <=> Pfam-category mapping file')
parser.add_argument('--genome-metadata', type=str,
                    help='Tab-delim metadata file that contains at least 2 columns: "accession" & "protein_count' + \
                         ' If provided, counts are normalized by total genes per genome')
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

def parse_map(infile):
    """
    Parse `UniRef <=> annotation` mapping file
    Return: {uniref => annot}
    """
    logging.info('Loading file: {}'.format(infile))
    idx = {}
    with _open(infile) as inF:
        for line in inF:
            line = _decode(line).rstrip().split('\t')
            for x in line[1:]:
                idx[x] = line[0]
    return idx    

def parse_cat(infile, annot_map, group_type):
    """
    parsing annotation category file
    Return: {uniref => cat}
    """
    logging.info('Loading file: {}'.format(infile))
    annot_cat = {}
    with _open(infile) as inF:
        for line in inF:
            line = _decode(line).rstrip().split('\t')
            if group_type == 'cog':
                annot_cat[line[0]] = line[2]
            elif group_type == 'pfam':
                annot_cat[line[0]] = line[1]
            else:
                msg = 'group_type "{}" not recognized'
                raise ValueError(msg.format(group_type))
    logging.info('Creating UniRef => category index...')
    annot_idx = {}
    for uniref,annot in annot_map.items():
        try:
            annot_idx[uniref] = annot_cat[annot]
        except KeyError:
            continue
    annot_map = None
    return annot_idx

def parse_annot(infile, annot_map):
    """
    Parsing gene annotation file.
    Return: {genome : annotations}
            {genome : taxonomy} 
    """
    logging.info('Loading file: {}'.format(infile))
    annot = defaultdict(dict)
    genome = {}
    header = {}
    status = {'genes' : 0, 'annot' : 0, 'no_annot' : 0}
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            #if i > 50000:    # DEBUG!!
            #    break
            status['genes'] += 1
            # line parse
            line = _decode(line).rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            # genome
            genome_name = line[header['genome_name']]
            try:
                _ = genome[genome_name]
            except KeyError:
                genome[genome_name] = [line[header['domain']],
                                       line[header['phylum']],
                                       line[header['class']],
                                       line[header['order']],
                                       line[header['family']],
                                       line[header['genus']],
                                       line[header['species']]]
            # annotations
            annotation = line[header['annotation']]
            ## change annotation to category?
            if annot_map is not None:
                try:
                    annotation = annot_map[annotation]
                    status['annot'] += 1
                except KeyError:
                    status['no_annot'] += 1
                    continue
            try:
                annot[genome_name][annotation] += 1
            except KeyError:
                annot[genome_name][annotation] = 1
    # status
    logging.info('No. of genomes: {}'.format(len(genome.keys())))
    logging.info('No. of genes: {}'.format(status['genes']))    
    if annot_map is not None:
        logging.info('No. genes w/ annotations: {}'.format(status['annot']))        
    return {'annotations' : annot, 'genomes' : genome}

def get_all_annots(annotations):
    """
    All annotations
    """
    all_annots = set()
    for genome in annotations.keys():
        for annot_name in annotations[genome].keys():
            all_annots.add(annot_name)
    logging.info('  No. of annotation columns: {}'.format(len(all_annots)))
    return sorted(all_annots)

def write_annot_table(annots, meta):
    """
    Writing table of annotations
    """
    logging.info('Writing table to STDOUT...')
    # all annotations
    all_annots = get_all_annots(annots['annotations'])
    header = ['genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    print('\t'.join(header + all_annots))
    for genome,taxonomy in annots['genomes'].items():
        # norm
        if meta is not None:
            try:
                protein_count = float(meta[genome])
            except KeyError:
                msg = 'Cannot find "{}" in metadata'
                raise KeyError(msg.format(genome))
        else:
            protein_count = None    
        # counts
        annot_cnts = []
        for annot in all_annots:
            try:
                x = annots['annotations'][genome][annot]
            except KeyError:
                x = 0
            # norm
            if protein_count is not None:
                x = round(x / protein_count, 9)
            annot_cnts.append(x)
        # writing line
        annot_cnts = [str(x) for x in annot_cnts]
        print('\t'.join([genome] + taxonomy + annot_cnts))

def parse_meta(infile):
    if infile is None:
        return None
    logging.info('Loading file: {}'.format(infile))
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
                accession = line[header['accession']]
            except KeyError:
                raise KeyError(msg.format('accession'))
            try:
                protein_count = line[header['protein_count']]
            except KeyError:
                raise KeyError(msg.format('protein_count'))
            meta[accession] = protein_count
    return meta

def main(args):
    # parsing metadata (if needed)
    meta = parse_meta(args.genome_metadata)
    # parsing mapping file (if needed)
    annot_map = None
    if args.content != 'uniref':
        # mapping file
        if args.content.startswith('cog'):
            annot_map = parse_map(args.cog_map)
            if args.content.endswith('categories'):
                annot_map = parse_cat(args.cog_cat, annot_map, 'cog')
        if args.content.startswith('pfam'):
            annot_map = parse_map(args.pfam_map)
            if args.content.endswith('categories'):
                annot_map = parse_cat(args.pfam_cat, annot_map, 'pfam')
    # loading annotations
    annots = parse_annot(args.annotations, annot_map)
    # writing table
    write_annot_table(annots, meta)
        

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
