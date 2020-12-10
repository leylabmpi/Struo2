#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import gzip
import argparse
import logging

desc = 'Filtering the cluster reps to just those that lack current annotations'
epi = """DESCRIPTION:
Filtering cluster reps to just those that are lacking any 
annotation data.
Output (filtered fasta) written to STDOUT
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('cluster_membership', metavar='cluster_membership', type=str,
                    help='mmseqs cluster membership file (format: cluster_rep<tab>cluster_member)')
parser.add_argument('query_hits', metavar='query_hits', type=str,
                    help='blast-formatted table of hits (cluster_rep <=> target_db_seqs)')
parser.add_argument('cluster_reps_fasta', metavar='cluster_reps_aa', type=str,
                    help='mmseqs cluster representatives fasta file')
parser.add_argument('--hit-columns', type=str,
                    default='qseqid,sseqid,evalue,pident,alnlen,slen',
                    help='Hit table output columns (default: %(default)s)')      
parser.add_argument('--min-pident', type=float, default=90,
                    help='Min % identity of hit (default: %(default)s)')
parser.add_argument('--min-cov', type=float, default=80,
                    help='Min % alignment coverage of subject sequence length (default: %(default)s)')
parser.add_argument('--reps-metadata', type=str, default=None,
                    help='Cluster reps metadata in order to provide more summary info about filtering (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def _open(infile):
    if infile.endswith('.gz'):
        return gzip.open(infile, 'rb')
    else:
        return open(infile)

def read_membership(infile):
    """
    Reading in cluster membership table: cluster_rep<tab>cluster_member
    Return: 
      dict: {cluster_member : cluster_id} 
    """
    logging.info('Reading file: {}'.format(infile))
    mem = {}
    with _open(infile) as inF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if len(line) < 2:
                continue
            mem[line[1]] = line[0]
    logging.info('  No. of cluster members: {}'.format(len(mem.keys())))
    n_clst = len(set(mem.values()))
    logging.info('  No. of clusters: {}'.format(n_clst))
    return mem

def read_hits(infile, mem, colnames, min_pident=0, min_cov=0):
    """ 
    Loading query hits.
    Return: 
      set(cluster rep)  # clusters with hits
    """
    logging.info('Loading hits table...')    
    clusts = []
    idx = {x:i for i,x in enumerate(colnames.split(','))}
    with _open(infile) as inF:
        for i,line in enumerate(inF):            
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            if len(line) < 2:
                msg = 'Line {}: <2 values in hits table'
                raise ValueError(msg.format(i+1))
            else:
                # filtering to just acceptable annotations
                ## percent identity
                pident = 0
                try:
                    pident = float(line[idx['pident']])
                except KeyError:
                    pass
                if pident < min_pident:
                    continue
                ## coverage of target seq
                cov = 0
                try:
                    cov = float(line[idx['slen']]) / float(line[idx['alnlen']]) * 100
                except KeyError:
                    pass
                if cov < min_cov:
                    continue
                ## adding clusterID to set of genes w/ acceptable annotation
                qseqid = line[idx['qseqid']]
                try:
                    clusts.append(mem[qseqid])
                except KeyError:
                    msg = 'Cannot find "{}" in cluster membership'
                    raise KeyError(msg.format(qseqid))
    clusts = set(clusts)
    msg = '  No of clusters with acceptable annotations: {}'
    logging.info(msg.format(len(clusts)))
    return clusts

def filter_fasta(infile, clust_w_annot, meta=None):
    """
    Filtering out sequences that are in the set(clust_w_annot)
    Return:
      None
    """
    logging.info('Filtering input fasta...')
    cnts = {'all' : 0, 'filtered' : 0, 'kept' : 0}
    genome_cnt = {}
    to_keep = False
    with _open(infile) as inF:
        for line in inF:            
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip()
            if line.startswith('>'):
                cnts['all'] += 1
                # already has annotation?
                line = line.lstrip('>').rstrip().split(' ')[0]                
                if line in clust_w_annot:
                    to_keep = False
                    cnts['filtered'] += 1
                else:
                    to_keep = True
                    cnts['kept'] += 1
                    print('>' + line)
                # metadata stats
                if meta is not None and to_keep is True:
                    try:
                        genome = meta[line]
                    except KeyError:
                        genome = 'OTHER'
                    try:
                        genome_cnt[genome] += 1
                    except KeyError:
                        genome_cnt[genome] = 1
            elif to_keep == True:
                print(line)
    # status
    logging.info('  No. of total seqs: {}'.format(cnts['all']))
    logging.info('  No. of filtered seqs: {}'.format(cnts['filtered']))
    logging.info('  No. of retained seqs: {}'.format(cnts['kept']))
    ## w/ metadata
    if meta is not None:
        msg = '    No. of retained seqs for {}: {}'
        for genome,cnt in genome_cnt.items():
            logging.info(msg.format(genome, cnt))

def read_metadata(infile):
    logging.info('Reading file: {}'.format(infile))
    header = {}
    meta = {}
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            seqid = line[header['seq_uuid']]
            genome = line[header['genome_name']]
            meta[seqid] = genome
    return meta
    
def main(args):
    # cluster metadata
    if args.reps_metadata is not None:
        meta = read_metadata(args.reps_metadata)
    else:
        meta = None
    
    # determining which clusters already have acceptable annotations
    clust_w_annot = read_hits(args.query_hits,
                              read_membership(args.cluster_membership),
                              colnames = args.hit_columns,
                              min_pident = args.min_pident,
                              min_cov = args.min_cov)
    # filtering fasta to just those lacking acceptable annotations
    filter_fasta(args.cluster_reps_fasta, clust_w_annot, meta)
    
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

