#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import logging
import argparse
from pprint import pprint
from collections import defaultdict

desc = 'Propagating all BLAST-table-formatted annotations to all genes in each mmseqs cluster'
epi = """DESCRIPTION:
Taking DIAMOND/mmseqs hits of gene cluster reps to a UniRef DB and applying them
to each member of the cluster.

Table of gene metadata written to STDOUT
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
                    help='diamond/mmseqs output columns (default: %(default)s)')      
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

def _open(infile):
    if infile.endswith('.gz'):
        return gzip.open(infile, 'rb')
    else:
        return open(infile)

def load_query_hits(infile, colnames, min_pident=0, min_cov=0):
    """ 
    Loading query hits
    Return: {query => cluster_rep : hit => unirefID}  
    """
    logging.info('Loading hits table...')    
    dmnd = {}
    skipped = {'pident' : 0, 'cov' : 0}
    idx = {x:i for i,x in enumerate(colnames.split(','))}
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
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
    """ 
    Loading cluster membership, formatted by 'mmseqs createtsv'
    Return: {cluster_member : cluster_rep} 
    """
    logging.info('Loading cluster membership...')
    mmshp = {}    
    with _open(infile) as inF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')  # [member,clusterID]
            if line[0] == '':
                continue
            if len(line) < 2:
                raise ValueError('<2 values in cluster membership table')
            mmshp[line[1]] = line[0]
    logging.info('  No. of clusters: {}'.format(len(mmshp.keys())))
    return mmshp

def mmshp_annotate(mmshp, dmnd):
    """ 
    Converting {cluster_member : cluster_rep} to {cluster_member : target_hit_seq_id} 
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

def load_gene_metadata(infile):
    """ Loading gene metadata
    Return: {geneID : 'body' : [metadata], header : header_line]
    """
    logging.info('Loading gene metadata...')
    colnames = {}
    meta = defaultdict(dict)
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            if infile.endswith('.gz'):
                line = line.decode('utf-8')                
            line = line.rstrip().split('\t')
            if i == 0:
                meta['header'] = line
                colnames = {x.lower():ii for ii,x in enumerate(line[1:])}
                continue
            if line[0] == '':
                continue
            if len(line) < 2:
                raise ValueError('<2 values in gene metadata table')
            meta['body'][line[0]] = line[1:]
    return meta, colnames

def propagate_info(infile, outfile, mmshp, meta, colnames,
                   keep_unclassified=True, seq_len_multi=1):
    """
    sequence header: gene_family|gene_length|taxonomy
    gene_len = nucleotide length (hence, the seq_len_multi)
    """
    outdir = os.path.split(outfile)[0]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    to_skip = False
    skipped = 0
    kept = 0
    seq = {}
    seq_header = None
    rename_idx = {}
    with _open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            if line.startswith('>'):      # sequence header
                # writing previous seq
                try:
                    seq = seq[seq_header]
                    seq_len = len(seq.rstrip('*')) * seq_len_multi                    
                    seq_header = seq_header.format(seq_len)
                    outF.write('\n'.join([seq_header, seq]) + '\n')
                    rename_idx[uuid] = [seq_header.lstrip('>'), seq_len]
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
                    m = meta['body'][uuid] 
                except KeyError:
                    msg = 'Cannot find "{}" in the gene metadata'
                    raise KeyError(msg.format(uuid))
                # skip/keep sequence?
                if keep_unclassified is False and sseqid.startswith('unclassified'):
                    to_skip = True
                    skipped += 1
                    continue
                else:                    
                    genus = m[colnames['genus']]
                    species = m[colnames['species']]                
                    taxid = m[colnames['taxid']]
                    species = '__'.join([species, taxid])
                    taxonomy = ';'.join([genus, species])
                    seq_header = '|'.join(['>' + sseqid, '{}', taxonomy])
                    seq[seq_header] = ''
                    kept += 1
            elif to_skip is False:       # sequence
                seq[seq_header] += line.rstrip()
        # last seq
        try:
            seq_len = len(seq[seq_header])
        except KeyError:
            seq_len = 0
        if to_skip is False and seq_len > 0:
            seq = seq[seq_header]
            seq_len = len(seq.rstrip('*')) * seq_len_multi
            seq_header = seq_header.format(seq_len)
            outF.write('\n'.join([seq_header, seq]) + '\n')        
            rename_idx[uuid] = [seq_header, seq_len]

    # status
    logging.info('File written: {}'.format(outfile))
    logging.info('  No. of seqs written: {}'.format(kept))
    if keep_unclassified is False:
        logging.info('  No. of seqs excluded due to unclassified: {}'.format(skipped))
        
    return rename_idx

# def check_meta(meta, mmshp):
#     missing = []
#     for gene_id,m in meta['body'].items():
#         try:
#             _ = mmshp[gene_id]
#         except KeyError:
#             missing.append(gene_id)
#     if len(missing) > 0:
#         msg = '-- Genes in metadata but missing from membership --'
#         logging.warning(msg)
#         logging.warning('Number: {}'.format(len(missing)))
#         logging.warning('Genes:')
#         print('\n'.join(missing))
#         sys.exit(1)

def write_table(meta, mmshp, humann_name_idx, keep_unclassified=False):
    """
    writing updated metadata table
    """
    logging.info('Writing gene metadata to stdout...')
    skipped = 0
    print('\t'.join(meta['header'] + ['annotation', 'humann_seq_id', 'seq_len_nuc']))
    for seq_id,m in meta['body'].items():
        try:
            annot = mmshp[seq_id]
        except KeyError:
            msg = 'Cannot find "{}" in cluster membership'
            raise KeyError(msg.format(seq_id))
        try:
            seq_header,seq_len = humann_name_idx[seq_id]
        except KeyError:
            if annot.startswith('unclassified'):
                seq_header = 'None'
            else:
                msg = 'Cannot find "{}" in humann name index'
                raise KeyError(msg.format(seq_id))
        if keep_unclassified is False and annot.startswith('unclassified'):
            skipped += 1
        else:
            print('\t'.join([seq_id] + m + [annot, seq_header, str(seq_len)]))
    if keep_unclassified is False:
        logging.info('  No. of genes excluded due to unclassified: {}'.format(skipped))
        
def main(args):
    # load query hits (eg., from diamond)
    dmnd = load_query_hits(args.diamond_hits, args.dmnd_columns,
                           min_pident=args.min_pident, min_cov=args.min_cov)
    # load cluster membership
    ## propagate cluster-level annotations to all members of each cluster
    mmshp = load_cluster_mmshp(args.cluster_membership)
    mmshp = mmshp_annotate(mmshp, dmnd)
    dmnd = None
    # load genes metadata
    meta,meta_columns = load_gene_metadata(args.genes_names)
    # checking that genes are present in membership
    #check_meta(meta, mmshp)
    
    # for each gene [fna, faa]:
    ## apply information to sequence header
    humann_name_idx = propagate_info(args.genes_fasta_prot, args.out_prot, mmshp, meta,
                                     meta_columns, keep_unclassified=args.keep_unclassified,
                                     seq_len_multi=3)
    if args.in_nuc is not None:
        propagate_info(args.in_nuc, args.out_nuc, mmshp, meta,
                       meta_columns, keep_unclassified=args.keep_unclassified,
                       seq_len_multi=1)
    # writing updated metadata table
    write_table(meta, mmshp, humann_name_idx, keep_unclassified=args.keep_unclassified)
    
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
