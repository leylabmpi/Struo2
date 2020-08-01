#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import logging
import argparse
from pprint import pprint

desc = 'Propagating all DIAMOND-based annotations to all genes in each mmseqs cluster'
epi = """DESCRIPTION:
Taking DIAMOND hits of gene cluster reps to a UniRef DB and applying them
to each member of the cluster.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('diamond_hits', metavar='diamond_hits', type=str,
                    help='tab-delim table of diamond hits')
parser.add_argument('genes_fasta_prot', metavar='genes_fasta_prot', type=str,
                    help='Genes in amino acid fasta format')
parser.add_argument('genes_names', metavar='genes_names', type=str,
                    help='Gene metadata such as origin genome & taxonomy')
parser.add_argument('cluster_membership', metavar='cluster_membership', type=str,
                    help='Cluster membership of each genome')
parser.add_argument('--in-nuc', type=str, default=None,
                    help='Input nucleotide fasta (default: %(default)s)')
parser.add_argument('--out-nuc', type=str, default='annotated.fna',
                    help='Output nucleotide fasta (default: %(default)s)')
parser.add_argument('--out-prot', type=str, default='annotated.faa',
                    help='Output amino acid fasta (default: %(default)s)')
parser.add_argument('--dmnd-columns', type=str, default='qseqid,sseqid,evalue,pident,alnlen,slen',
                    help='Diamond output columns (default: %(default)s)')      
parser.add_argument('--meta-columns', type=str,
                    default='uuid,gene_name,domain,phylum,class,order,family,genus,species,taxID,genomeID,genome_len',
                    help='Gene metadata table columns (default: %(default)s)')
parser.add_argument('--keep-unclassified', action='store_true', default=False,
                    help='Keep gene clusters without DIAMOND hit? (default: %(default)s)')
parser.add_argument('--min-pident', type=float, default=50,
                    help='Min % identity of hit (default: %(default)s)')
parser.add_argument('--min-cov', type=float, default=80,
                    help='Min % alignment coverage of subject sequence length (default: %(default)s)')
parser.add_argument('--threads', type=int, default=1,
                    help='Threads used for diamond (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def load_dmnd_hits(infile, colnames, min_pident=0, min_cov=0):
    """ Loading query hits
    Return: {query => cluster_rep : hit => unirefID}  
    """
    logging.info('Loading hits table...')
    dmnd = {}
    skipped = {'pident' : 0, 'cov' : 0}
    idx = {x:i for i,x in enumerate(colnames.split(','))}
    with open(infile) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            if len(line) < 2:
                raise ValueError('Line {}: <2 values in hits table'.format(i+1))
            else:
                # filtering
                ## percent identity
                pident = 0
                try:
                    pident = float(line[idx['pident']])
                except KeyError:
                    pass
                if pident < min_pident:
                    skipped['pident'] += 1
                    continue
                ## coverage of target seq
                cov = 0
                try:
                    cov = float(line[idx['slen']]) / float(line[idx['alnlen']]) * 100
                except KeyError:
                    pass
                if cov < min_cov:
                    skipped['cov'] += 1
                    continue
                # loading
                qseqid = line[idx['qseqid']]
                dmnd[qseqid] = line[idx['sseqid']]
    logging.info('  No. of excluded hits w/ < min-pident: {}'.format(skipped['pident']))
    logging.info('  No. of excluded hits w/ < min-cov: {}'.format(skipped['cov']))
    logging.info('  No. of annotated genes: {}'.format(len(dmnd.keys())))
    return dmnd
    
def load_cluster_mmshp(infile):
    """ Loading cluster membership, formatted by 'mmseqs createtsv'
    Return: {cluster_member : cluster_rep} 
    """
    logging.info('Loading cluster membership...')
    mmshp = {}    
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')  # [member,clusterID]
            if line[0] == '':
                continue
            if len(line) < 2:
                raise ValueError('<2 values in cluster membership table')
            mmshp[line[1]] = line[0]
    logging.info('  No. of clusters: {}'.format(len(mmshp.keys())))
    return mmshp

def mmshp_annotate(mmshp, dmnd):
    """ Converting {cluster_member | cluster_rep} to {cluster_member : dmnd_sseqid} 
    """
    logging.info('Propogating sseqid values to all members of each cluster...')
    # cluster representatives with diamond hit
    rep_idx = {}   # cluster_id : sseqid
    annot = 0
    uncls = 0
    for cluster_member,cluster_rep in mmshp.items():
        try:
            rep_idx[cluster_member] = dmnd[cluster_rep]
            annot += 1
        except KeyError:
            rep_idx[cluster_member] = 'unclassified'
            uncls += 1
    logging.info('  No. of members: {}'.format(len(rep_idx.keys())))
    logging.info('    No. of annotated: {}'.format(annot))
    logging.info('    No. of unclassifieds: {}'.format(uncls))
    return rep_idx

def load_gene_metadata(infile, colnames):
    """ Loading gene metadata
    Return: {geneID : [metadata]} 
    """
    logging.info('Loading gene metadata...')
    colnames = colnames.split(',')
    meta = {}            
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            if len(line) < 2:
                raise ValueError('<2 values in gene metadata table')
            if len(line) != len(colnames):
                print(line)
                msg = 'No. of columns in metadata table != --meta-columns'
                raise ValueError(msg)
            meta[line[0]] = line[1:]
    return meta

def propagate_info(infile, outfile, mmshp, meta, colnames, keep_unclassified=True, seq_len_multi=1):
    """
    sequence header: gene_family|gene_length|taxonomy
    gene_len = nucleotide length (hence, the seq_len_multi)
    """
    idx = {x.lower():i for i,x in enumerate(colnames.split(',')[1:])}
    to_skip = False
    skipped = 0
    kept = 0
    seq = {}
    seq_header = None
    with open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            if line.startswith('>'):      # sequence header
                # writing last seq
                try:
                    seq = seq[seq_header]                    
                    seq_header = seq_header.format(len(seq.rstrip('*')) * seq_len_multi)
                    outF.write('\n'.join([seq_header, seq]) + '\n')
                    seq = {}
                    seq_header = None
                except KeyError:
                    pass
                # formatting for this seq
                to_skip = False
                uuid = line.lstrip('>').rstrip()
                # getting sseqid 
                try:
                    sseqid = mmshp[uuid]
                except KeyError:
                    msg = 'Cannot find "{}" in cluster membership'
                    raise KeyError(msg.format(uuid))
                # getting metadata
                try:
                    m = meta[uuid] 
                except KeyError:
                    msg = 'Cannot find "{}" in the gene metadata'
                    raise KeyError(msg.format(uuid))
                # skip/keep sequence?
                if keep_unclassified is False and sseqid.startswith('unclassified'):
                    to_skip = True
                    skipped += 1
                    continue
                else:                    
                    genus = m[idx['genus']]
                    species = m[idx['species']]                
                    taxid = m[idx['taxid']]
                    species = '__'.join([species, taxid])
                    taxonomy = ';'.join([genus, species])
                    seq_header = '|'.join(['>' + sseqid, '{}', taxonomy])
                    seq[seq_header] = ''
                    kept += 1
            elif to_skip is False:      # sequence
                seq[seq_header] += line.rstrip()
        # last seq
        if to_skip is False and len(seq[seq_header]) > 0:
            seq = seq[seq_header]
            seq_header = seq_header.format(len(seq.rstrip('*')) * seq_len_multi)
            outF.write('\n'.join([seq_header, seq]) + '\n')        

    # status
    logging.info('File written: {}'.format(outfile))
    logging.info('  No. of seqs written: {}'.format(kept))
    if keep_unclassified is False:
        logging.info('  No. of seqs excluded due to unclassified: {}'.format(skipped))

def write_table(mmshp, meta, colnames):
    colnames = colnames.split(',')[1:]
    idx = {x:i for i,x in enumerate(colnames)}
    print('\t'.join(['gene_uuid', 'annotation'] + colnames))
    for cluster_member,sseqid in mmshp.items():
        m = meta[cluster_member]
        line = [cluster_member, sseqid] + [m[idx[x]] for x in colnames]
        print('\t'.join(line))
        
def main(args):
    # load diamond hits
    dmnd = load_dmnd_hits(args.diamond_hits, args.dmnd_columns,
                          min_pident=args.min_pident, min_cov=args.min_cov)
    # load cluster membership
    mmshp = load_cluster_mmshp(args.cluster_membership)
    mmshp = mmshp_annotate(mmshp, dmnd)
    dmnd = None
    # load genes metadata
    meta = load_gene_metadata(args.genes_names, args.meta_columns)
        
    # for each gene [fna, faa]:
    ## apply information to sequence header
    propagate_info(args.genes_fasta_prot, args.out_prot, mmshp, meta,
                   args.meta_columns, keep_unclassified=args.keep_unclassified,
                   seq_len_multi=3)
    if args.in_nuc is not None:
        propagate_info(args.in_nuc, args.out_nuc, mmshp, meta,
                       args.meta_columns, keep_unclassified=args.keep_unclassified,
                       seq_len_multi=1)

    # create dict: {gene : [clusterID, dmnd-hit, gene_metadata]}
    write_table(mmshp, meta, args.meta_columns)
    
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
