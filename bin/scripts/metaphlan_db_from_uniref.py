#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import bz2
import pickle
import copy
import logging
import argparse
import collections
import numpy as np
from pprint import pprint

desc = 'Identify species-specific genes'
epi = """DESCRIPTION:
Updating a metaphlan3 database (pkl & the nucleotide fasta).
Using genes annotated for humann3 as cluster.
Instead of using NCBI taxonomy, GTDB taxonomy will be used.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('gene_fasta', metavar='gene_fasta', type=str,
                    help='all genes in fasta format (nucleotide)')
parser.add_argument('gene_metadata', metavar='gene_metadata', type=str,
                    help='tab-delim table of gene metadata with humann annotations')
parser.add_argument('metaphlan_pkl', metavar='metaphlan_pkl', type=str,
                    help='Metaphlan database pkl file to view/update')
parser.add_argument('metaphlan_fasta', metavar='metaphlan_fasta', type=str,
                    help='Metaphlan database fasta file')
parser.add_argument('gtdb_taxdump_names', metavar='gtdb_taxdump_names', type=str,
                    help='GTDB taxdump names.dmp file')
parser.add_argument('--pkl-n', type=int, default=0,
                    help='Number of entries to view in the pkl database file. If <1, no viewing (default: %(default)s)')
parser.add_argument('--pkl-viruses', action='store_true', default=False,
                    help='Include viruses when viewing the pkl database file? (default: %(default)s)')
parser.add_argument('--ext-abs-cutoff', type=int, default=10,
                    help='Max number of external genomes (default: %(default)s)')
parser.add_argument('--ext-rel-cutoff', type=float, default=0.1,
                    help='Max fraction of external genomes to total genomes (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def _open(infile):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, 'rb')
    elif infile.endswith('.gz'):
        return gzip.open(infile, 'rb')
    else:
        return open(infile)

def _decode(line, infile):
    """
    Decoding input, depending on the file extension
    """
    if infile.endswith('.gz') or infile.endswith('.bz2'):
        line = line.decode('utf-8')
    return line
    
def view_mphn_pkl(infile, N=10, viruses=False):
    """
    Reading in metaphlan database pkl file and printing
    # Sequence formatting: (NCBI_taxid)(UniRef90_cluster)(CDS_name)
    #   Note that the UniRef90 clusterID does NOT include "Uniref90_"
    """
    logging.info('Reading in: {}'.format(infile))
    db = pickle.load(bz2.open(infile, 'r'))
    m = '-- taxonomy of geach genome: [taxonomy, tax-ids_per-tax-level, genome_length] --'
    logging.info(m)
    cnt = 0
    for tax in db['taxonomy'].keys():
        if viruses is False and 'k__Viruses' in tax or 'k__Viroids' in tax:
            continue
        else:
            cnt += 1
        if cnt > N:
            break
        print([tax, db['taxonomy'][tax]])

    m = '-- marker gene: [gene_id : {ext: [], score: float, clade: str, len: gene_len, taxon: str}] --'
    logging.info(m)
    logging.info('  "gene_id" matches the sequence headers of the marker fasta')
    logging.info('  ext: GCA-IDs of non-species genomes that have markers')
    cnt = 0
    for mkr in db['markers'].keys():
        if viruses is False and 'k__Viruses' in db['markers'][mkr]['taxon']:            
            continue
        else:
            cnt += 1
        if cnt > N:
            break
        print([mkr, db['markers'][mkr]])
    exit(0)

def get_tax_ids(tax, idx):
    """
    getting tax_ids based on taxon names
    """
    tax_ids = []
    for x in tax:
        try:
            tax_ids.append(idx[x])
        except KeyError:
            msg = 'WARNING: {} not found in taxdump names file'
            logging.info(msg.format(x))
            tax_ids.append('')
    if len(tax) != len(tax_ids):
        tax = ','.join(tax)
        raise ValueError('Problem finding taxids for {}'.format(tax))
    return tax_ids
    
def create_metadata_indices(meta_file, tax_id_idx):
    """
    Loading Struo gene metadata info and creating an index.
    Args:
      meta_file : struo2 genes metadata file
      marker_len_idx : the length of the gene sequence      
    Return:
      return: {cluster_id : {(taxonomy) : { genome_id : {gene_ids : [gene_ids], genome_len : len}}}}
    """
    logging.info('Loading gene metadata...')
    # genes
    header = {}
    idx = collections.defaultdict(dict)
    tax_levs = ['domain','phylum','class','order','family','genus','species']
    with _open(meta_file) as inF:
        for i,line in enumerate(inF):
            line = _decode(line, meta_file)
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x.lower():i for i,x in enumerate(line)}
                continue
            # genes
            cluster_id = line[header['annotation']]
            gene_uuid = line[header['seq_uuid']]
            genome_id = line[header['genome_name']]
            genome_len = line[header['genome_length_bp']]
            tax = [line[header[x]] for x in tax_levs]
            tax_ids = get_tax_ids(tax, tax_id_idx)                
            genome_tax_id = line[header['taxid']]
            humann_id = line[header['humann_seq_id']]
            gene_len = line[header['seq_len_nuc']]
            ## adding to data structure
            try:
                # adding genes to cluster
                idx[cluster_id][genome_id]['gene_ids'].append([gene_uuid,
                                                               humann_id,
                                                               gene_len])
            except KeyError:
                # genome metadata + first gene of cluster
                idx[cluster_id] = {genome_id : {'gene_ids' : [[gene_uuid,
                                                               humann_id,
                                                               gene_len]],
                                                'genome_len' : genome_len,
                                                'taxonomy' : tax,
                                                'taxids' : tax_ids,
                                                'genome_tax_id' : genome_tax_id}}
    #metadata_summary(idx)
    return idx


def tax_gene_cnt(d):
    for tax,v1 in d.items():
        for genome_id,v2 in v1.items():
            pass
    

def cluster_tax_breadth(clust_idx, max_ext=10):
    """
    Determining the taxonomic breadth of each cluster.
    Filter multi-copy genes 
    Which are species specific?
    """
    logging.info('Determining taxonomic breadth of each cluster...')
    for clust_id,v1 in clust_idx.items():
        tax_cnt = {}
        for genome_id,v2 in v1.items():
            pass
            #n_genes = 

def filter_multicopy(clust_idx):
    for clust_id,v1 in clust_idx.items():
        genome_cnts = []
        for genome_id,v2 in v1.items():
            genome_cnts.append(len(v2['gene_ids']))
        print(genome_cnts)
    sys.exit()
    
            
def sum_stats(x, label):
    """
    Print stats of a distribution (vector)
    """
    y = '  {}: min={}, q1={}, mean={}, q3={}, max={}'
    y = y.format(label,
                 round(np.min(x),2),
                 round(np.quantile(x,0.25),2),
                 round(np.mean(x),2),
                 round(np.quantile(x,0.75),2),
                 round(np.max(x),2))
    logging.info(y)

def metadata_summary(idx):
    """
    Summarizing the gene metadata
    """
    tax_per_cluster = []
    genomes_per_tax = []
    genes_per_genome = []
    for cluster_id,v in idx.items():
        tax_per_cluster.append(len(v.keys()))
        for tax,vv in v.items():
            genomes_per_tax.append(len(vv.keys()))
            for genomeID,gene_ids in vv.items():
                genes_per_genome.append(len(set(gene_ids)))
    sum_stats(tax_per_cluster, 'Clades per cluster')
    sum_stats(genomes_per_tax, 'Gemomes per clade')
    sum_stats(genes_per_genome, 'Genes per genome')  

def get_full_tax(idx):
    """
    Getting full taxonomy for all genes.
    return: {(tax_lev1, tax_lev2, ...) : genome_ID : genome_len}
    """
    logging.info('Compiling the taxonomy for all genomes...')
    tax_idx = collections.defaultdict(dict)
    for cluster_id,v in idx.items():
        for tax,vv in v.items():
            for genome_id,x in vv.items():
                tax_idx[tax][genome_id] = x['genome_len']
    n_genomes = 0
    for tax,v in tax_idx.items():
        n_genomes += len(v.keys())
    logging.info('  Total number of genomes: {}'.format(n_genomes))
    # return
    return tax_idx

def filter_multi_copy_clusters(idx):
    """
    {cluster_id : {taxonomy : {genomeid : [gene_uuid,...]}}}
    """
    logging.info('Filtering out multi-copy genes...')
    clust_cnt = collections.defaultdict(dict)
    to_remove = []
    for cluster_id,v in idx.items():
        per_genome_copy = {}
        for tax,vv in v.items():
            for genome_id,x in vv.items():
                per_genome_copy[genome_id] = len(set(x['gene_ids']))
        # any multi-copy?
        if any([x > 1 for x in per_genome_copy.values()]):
            to_remove.append(cluster_id)
    for cluster_id in to_remove:
        idx.pop(cluster_id, None)
    # status
    logging.info('  Number of multi-copy clusters removed: {}'.format(len(to_remove)))
    logging.info('  Number of single-copy clusters remaining: {}'.format(len(idx.keys())))
    if len(idx.keys()) < 1:
        logging.info('Exiting due to a lack of clusters')
        sys.exit(0)
    metadata_summary(idx)
    return idx    

def count_clades(idx, tax_idx):
    """
    input: {(taxonomy) : {genome_id : [gene_IDs]}}
    return: {(taxonomy) : n_genomes_with_gene / tota_genomes
    """
    cnt = {}
    for tax,v in idx.items():
        # total genomes in clade
        total_genomes = len(tax_idx[tax].keys())
        # genomes with the gene in the clade        
        n_genomes = len(idx[tax].keys())
        cnt[tax] = [n_genomes / float(total_genomes), n_genomes, total_genomes]
    return cnt

def dominant_clade(idx, prev_cutoff=1):
    """
    Which is the dominant clade?
    """
    for k in sorted(idx, key=lambda item: item[0], reverse=True):
        if idx[k][0] < prev_cutoff:
            return None
        else:
            return k

def get_external(tax_cnt, idx, clade):
    """
    Getting all external genes
    """
    external = []
    for tax,cnt in tax_cnt.items():
        if tax == clade:
            continue
        if cnt[0] > 0:
            # getting "external" genomes
            for genome_id in idx[tax].keys():
                external.append(genome_id)
    return external
        
def gene_tax_breadth(idx, tax_idx, ext_abs_cutoff = 10, ext_rel_cutoff=0.1):
    """ 
    Determining the taxonomic breadth of each gene cluster (unirefID).
    Determining which to use as species-specific markers.
    Return: {cluster_id : {clade : '', ext : [], taxon : ''}}
    """
    logging.info('Determining the taxonomic breadth of each gene cluster...')
    skipped = 0
    breadth_idx = {}
    for cluster_id in idx.keys():
        # fraction of genomes in each clade with the gene
        tax_cnt = count_clades(idx[cluster_id], tax_idx)
        # selecting the dominant
        clade = dominant_clade(tax_cnt)
        # skipping if the cluster is not present in all of the dominant clade
        if clade is None:
            skipped += 1
            continue
        # getting any genomes external to the target clade
        ext = get_external(tax_cnt, idx[cluster_id], clade)
        # cutoff for determining clade specificity
        clade_n_genomes = tax_cnt[clade][1]
        if len(ext) > ext_abs_cutoff or len(ext) / float(clade_n_genomes) > ext_rel_cutoff:
            skipped += 1
        # loading info
        breadth_idx[cluster_id] = {'clade' : clade[-1],
                                   'ext' : ext,
                                   'taxon' : '|'.join(list(clade)),
                                   'len' : int(cluster_id.split('|')[1])}
    msg = '  Number of clusters skipped due to lack of species-specificity: {}'
    logging.info(msg.format(skipped))
    return breadth_idx

def read_pkl(infile):
    """
    Reading in the metaphlan database pkl file
    """
    # loading
    if infile is not None:
        logging.info('No pkl provided. Creating a new object')
        pkl = pickle.load(bz2.open(infile, 'r'))
    else:
        logging.info('Reading in: {}'.format(infile))
        # creating new object
        pkl = {'taxonomy' : {}, 'markers' : {}}
    return pkl

def add_taxonomy(tax_idx, pkl):
    """
    Adding taxonomy to the metaphlan pkl.
    Example additions: 
      db['taxonomy']['taxonomy of genome1'] = ('NCBI taxonomy id of genome1', length of genome1)
      db['taxonomy']['taxonomy of genome2'] = ('NCBI taxonomy id of genome1', length of genome2)
    """
    for tax,v in tax_idx.items():
        for genome_id,genome_len in v.items():
            T = '|'.join(list(tax) + ['t__' + genome_id])
            pkl['taxonomy'][T] = ('', int(genome_len))
    return pkl

def add_markers(idx, pkl):
    """
    Adding markers to the metaphlan pkl.
    """
    for cluster_id,v in idx.items():
        marker_name = cluster_id.split('|')[0]
        pkl['markers'][marker_name] = copy.deepcopy(v)
    return pkl

def create_tax_id_index(infile):
    """
    Creating an index from an taxdump names file
    Return: {tax_name : tax_id}
    """
    logging.info('Loading: {}'.format(infile))
    regex = re.compile('\t\|\t')
    idx = {}
    with _open(infile) as inF:
        for line in inF:
            line = _decode(line, infile)
            line = regex.split(line)
            idx[line[1].replace(' ', '_')] = int(line[0])
    return idx
        
def main(args):
    # viewing pkl file
    if args.pkl_n > 0:
        view_mphn_pkl(args.metaphlan_pkl, N=args.pkl_n, viruses=args.pkl_viruses)

    # creating indices for all required data
    ## creating {marker_uuid : gene_length} index
    #marker_len_idx = create_marker_len_idx(args.gene_fasta)
    ## creating {tax_name : taxid} index
    tax_id_idx = create_tax_id_index(args.gtdb_taxdump_names)    
    ## creating {unirefID : [marker_uuid]} index
    cluster_idx = create_metadata_indices(args.gene_metadata, tax_id_idx)

    # filtering out multi-copy gene cluster (changes in place)
    filter_multicopy(cluster_idx)
    
    # determining the taxonomic breadth of each cluster
    cluster_tax_breadth(cluster_idx)
    
    
    #-- TO HERE --#
    # change to per cluster:
    ## for the cluster, is it species-specific (with a fraction of non-species genes)?
    
    ## getting the full taxonomy (w/ genomeIDs) for each cluster
    tax_idx = get_full_tax(cluster_idx)
    ## just clusters that are single copy for each genome
    cluster_idx = filter_multi_copy_clusters(cluster_idx)
        
    # taxonomic breadth of each cluster_id (unirefID)
    tax_breadth = gene_tax_breadth(cluster_idx, tax_idx,
                                   ext_abs_cutoff=args.ext_abs_cutoff,
                                   ext_rel_cutoff=args.ext_rel_cutoff)

    # loading the pkl
    pkl = read_pkl(args.pkl)    
    ## adding taxonomy to the pkl
    pkl = add_taxonomy(tax_idx, pkl)
    ## adding markers
    pkl = add_markers(tax_breadth, pkl)
    # formatting the fasta
    pprint(pkl)
    sys.exit(1)
    
if __name__ == '__main__':
    args = parser.parse_args()
    # potential input
    # gene metadata (unirefID, taxonomy + taxID, genomeID)
    main(args)
