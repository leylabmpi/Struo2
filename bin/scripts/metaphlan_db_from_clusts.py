#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import bz2
import pickle
import copy
import logging
import pickle
import argparse
import collections
import numpy as np
from itertools import chain
from pprint import pprint
# 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants

desc = 'Create metaphlan database from gene cluster data'
epi = """DESCRIPTION:
"markers" definition
* single-copy across all genomes
* 1 representative sequence per clade per cluster
# PROBLEMS
* the distribution of marker genes must be known across all genomes in each species,
  but struo is designed to just use one genome per species
* all genomes from the GTDB would have to be included
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('gene_fasta', metavar='gene_fasta', type=str,
                    help='all genes in fasta format (nucleotide)')
parser.add_argument('gene_metadata', metavar='gene_metadata', type=str,
                    help='tab-delim table of gene metadata')
parser.add_argument('cluster_membership', metavar='cluster_membership', type=str,
                    help='gene cluster membeship (cluster_rep<tab>cluster_member)')
parser.add_argument('gtdb_taxdump_names', metavar='gtdb_taxdump_names', type=str,
                    help='GTDB taxdump names.dmp file')
parser.add_argument('gtdb_taxdump_nodes', metavar='gtdb_taxdump_nodes', type=str,
                    help='GTDB taxdump nodes.dmp file')
parser.add_argument('-o', '--output-prefix', type=str, default='metaphlan_db',
                    help='Output file prefix (default: %(default)s)')
parser.add_argument('--ext-cutoff', type=float, default=0.9,
                    help='Min fraction of clade-specificity (default: %(default)s)')
parser.add_argument('--metaphlan-pkl', type=str, default=None,
                    help='Metaphlan database pkl file to view')
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
    # searching for markers by cluster
    queries =[ '140__A1QYH8__BHO_0112400', '141__A1QYH8__X966_00125', '142__A1QYH8__BT0025']
    for q in queries:
        try:
            print([q, db['markers'][q]])
        except KeyError:
            pass
    sys.exit()

    # taxonomy
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
    # marker
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

def create_tax_id_index(names_file, nodes_file):
    """
    Creating an index from an taxdump names file
    Return: {tax_name : tax_id}
    """
    logging.info('Creating taxonomy graph...')
    G = nx.DiGraph()
    regex = re.compile('\t\|\t')
    # adding edges
    logging.info('Loading: {}'.format(nodes_file))
    with _open(nodes_file) as inF:
        for line in inF:
            line = _decode(line, nodes_file)
            line = regex.split(line)
            G.add_edge(int(line[0]), int(line[1]))
    # adding node attributes (tax names)
    logging.info('Loading: {}'.format(names_file))
    tax_idx = {}
    with _open(names_file) as inF:
        for line in inF:
            line = _decode(line, names_file)
            line = regex.split(line)
            node_id = int(line[0])
            tax = line[1].replace(' ', '_')
            G.nodes[node_id]['taxonomy'] = tax
            tax_idx[tax] = node_id
    return tax_idx, G

def create_membership_index(infile):
    """
    Reading in mmseqs cluster membership file: cluster_id<tab>gene_id
    Return: {gene_id : cluster_id}
    """
    logging.info('Loading: {}'.format(infile))
    idx = {}
    with _open(infile) as inF:
        for line in inF:
            line = _decode(line, infile)
            line = line.rstrip().split('\t')
            idx[line[1]] = line[0]
    return idx

def get_tax_ids(tax, tax_idx):
    """
    Getting tax_ids based on taxon names
    Args:
      tax : vector of taxonomic classifications
      idx : {taxonomic_name : taxid} 
    Return
      vector of taxids 
    """
    tax_ids = []
    for x in tax:
        try:
            tax_ids.append(tax_idx[x])
        except KeyError:
            msg = 'WARNING: {} not found in taxdump names file'
            logging.info(msg.format(x))
            tax_ids.append('')
    if len(tax) != len(tax_ids):
        tax = ','.join(tax)
        raise ValueError('Problem finding taxids for {}'.format(tax))
    return tax_ids

def create_metadata_index(infile, mem_idx, tax_idx, min_genomes=5):
    """
    Creating indices from gene metadata
    Args:
      infile : struo2 gene metadata table
      mem_idx : cluster membership index
      tax_idx : taxonomy index from taxdump names file
    Return:
      genome_idx : genome metadata
      clust_idx : cluster metadata
    """
    logging.info('Loading: {}'.format(infile))
    header = {}
    tax_levs = ['domain', 'phylum', 'class', 'order',
                'family', 'genus', 'species']
    genome_idx = {}
    clust_idx = collections.defaultdict(dict)
    clust_genome_idx = collections.defaultdict(dict)
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            line = _decode(line, infile)
            line = line.rstrip().split('\t')
            # header
            if i == 0:
                header = {x:i for i,x in enumerate(line)}
                continue
            # body
            seq_uuid = line[header['seq_uuid']]
            seq_name = line[header['seq_orig_name']]
            clust_id = mem_idx[seq_uuid]
            tax = [line[header[x]] for x in tax_levs]
            tax_ids = get_tax_ids(tax, tax_idx)                
            genome_taxid = line[header['taxid']]
            genome_name = line[header['genome_name']]            
            genome_length_bp = line[header['genome_length_bp']]
            ## genome metadata
            genome_idx[genome_name] = {'taxid' : genome_taxid,
                                       'len' : genome_length_bp,
                                       'taxonomy' : tuple(tax),
                                       'taxids' : tax_ids}
            ## cluster info
            clust_idx[clust_id][tuple(tax_ids)] = {seq_uuid : {'seq_name' : seq_name,
                                                               'genome_name' : genome_name}}
            ## cluster copies per genome
            try:
                clust_genome_idx[clust_id][genome_name] += 1
            except KeyError:
                clust_genome_idx[clust_id][genome_name] = 1
    logging.info('No. of clusters: {}'.format(len(clust_idx.keys())))
        
    # filtering multi-copy clusters
    filtered = {'low-prev' : 0, 'multi-copy' : 0}
    for clust_id in clust_genome_idx.keys():
        # single copy?
        clust_cnt = [clust_genome_idx[clust_id][x] for x in clust_genome_idx[clust_id].keys()]
        if any([x > 1 for x in clust_cnt]):
            clust_idx.pop(clust_id, None)
            filtered['multi-copy'] += 1
    logging.info('  No. of multi-copy clusters filtered: {}'.format(filtered['multi-copy']))
    logging.info('  No. of clusters remaining: {}'.format(len(clust_idx.keys())))
    
    # return
    return genome_idx, clust_idx

def get_marker_prev(dom_clade, clade_genome, tax_graph):
    """
    Determining the markers' prevalence across all genomes for the dominant clade.
    Args:
      dom_clade : tax_id
    """
    # all genomes in the dominant clade    
    all_clade_genomes = [desc for desc in descendants(tax_graph, dom_clade[-1])]
    # all genomes for the cluster in the clade
    cluster_clade_genomes = clade_genome[dom_clade]
    prev = len(cluster_clade_genomes) / float(len(all_clade_genomes))
    return prev

def clade_specific(genome_idx, cluster_idx, tax_graph, ext_cutoff=0.9, prev_cutoff=0):
    """
    Determining which clusters are clade-specific (or nearly specific to 1 clade)
    Args:
      genome_idx : genome metadata 
      cluster_idx : cluster metadata
      ext_cutoff : min fraction of genes in cluster that must be cluster-specific
      prev_cutoff : marker must be found in >= X% of all genomes in the clade
    Return:
      clade_idx : clade metadata
      clade_seq_idx : {seq_uuid : clade}
    """
    logging.info('Finding clade-specific clusters...')
    logging.info('  Using clade specificity cutoff: >={}'.format(ext_cutoff))
    # creating clade index
    to_rm = []
    dom_clade_idx = collections.defaultdict(dict)
    clade_seq_idx = collections.defaultdict(dict)
    for clust_id,v1 in cluster_idx.items():
        # clade info
        clade_cnt = {}
        clade_genome = {}
        for clade,v2 in v1.items():
            # number of genes per clade (per cluster)
            clade_cnt[clade] = len(v2.keys())            
            # genomes per clade (per cluster)
            clade_genome[clade]  = [v['genome_name'] for k,v in v2.items()]            
        clade_cnt = collections.Counter(clade_cnt)
        # dominant clade for cluster?
        dom = clade_cnt.most_common(1)[0]
        total = sum(clade_cnt.values())
        # filter if not clade specific
        if dom[1] / float(total) < ext_cutoff:
            to_rm.append(clust_id)
            continue
        # is marker found in enough of the entire clade?
        if get_marker_prev(dom[0], clade_genome, tax_graph) < prev_cutoff:
            to_rm.append(clust_id)
            continue            
        # taxa (genomes) external to dominant clade
        if dom[1] != total:
            ext = [clade_genome[x] for x in clade_cnt.keys() if x != dom[0]]
            ext = list(chain.from_iterable(ext))
        else:
            ext = []
        # clade index
        dom_clade_idx[clust_id] = {'clade' : dom[0][-1],
                                   'ext' : ext,
                                   'taxon' : dom[0]}
        # using first seq listed for clust-clade marker seq
        for clade,v2 in v1.items():
            seq_uuid = list(v2.keys())[0] 
            clade_seq_idx[seq_uuid] = {'clust_id' : clust_id,
                                       'clade' : clade}
    # filtering clade non-specific clusters
    for clust_id in to_rm:
        cluster_idx.pop(clust_id, None)
        
    # status
    logging.info('  No. of clade-specific clusters: {}'.format(len(cluster_idx.keys())))
    return dom_clade_idx, clade_seq_idx


def write_fasta(fasta_file, cluster_idx, genome_idx,
                clade_seq_idx, out_prefix = 'file.fna'):
    """
    Filtering & renaming marker sequences.
    Rename format: genome-taxid__cluster-id__gene-name
    """
    logging.info('Creating marker sequence file...')
    outfile = out_prefix + '.fna'
    logging.info('  Loading: {}'.format(fasta_file))
    status = {'filtered' : 0, 'clusts' : {}, 'seqs' : 0}
    with _open(fasta_file) as inF, open(outfile, 'w') as outF:
        write_seq = False
        seq_uuid = None
        clust_id = None
        for line in inF:
            line = _decode(line, fasta_file)
            # renaming sequence 
            if line.startswith('>'):
                write_seq = False
                seq_uuid = line.lstrip('>').rstrip()
                # in clade_seq_idx?
                try:
                    clade = clade_seq_idx[seq_uuid]['clade']
                    clust_id = clade_seq_idx[seq_uuid]['clust_id']
                except KeyError:
                    status['filtered'] += 1
                    continue
                # is cluster in cluster_idx?
                try:
                    gene_info = cluster_idx[clust_id][clade][seq_uuid]
                except KeyError:
                    status['filtered'] += 1
                    continue
                # getting genome info
                try:
                    genome_info = genome_idx[gene_info['genome_name']]
                except KeyError:
                    msg = 'Cannot find {} in genome metadata'
                    raise KeyError(msg.format(gene_info['genome_name']))
                # new seq header
                genome_taxid = genome_info['taxid']
                outF.write('>' + seq_uuid + '\n')
                status['clusts'][clust_id] = 1                
                status['seqs'] += 1
                write_seq = True
            elif write_seq == True:
                outF.write(line)
                try:
                    cluster_idx[clust_id][clade][seq_uuid]['len'] += len(line.rstrip())
                except KeyError:
                    cluster_idx[clust_id][clade][seq_uuid]['len'] = len(line.rstrip())
    
    # status
    logging.info('  No. of seqs filtered: {}'.format(status['filtered']))
    logging.info('  No. of clusters included: {}'.format(len(status['clusts'].keys())))
    logging.info('  No. of sequences written: {}'.format(status['seqs']))
    logging.info('  File written: {}'.format(outfile))
                
def write_pkl(genome_idx, cluster_idx, clade_idx, out_prefix = 'table'):
    """
    writing metaphlan pkl
    """
    logging.info('Creating metaphlan pkl file...')
    db = collections.defaultdict(dict)
    
    # genome taxonomy
    for genome_name,v in genome_idx.items():
        # genome taxonomy
        taxonomy = ['t__' + genome_name]
        taxonomy += genome_idx[genome_name]['taxonomy']
        taxonomy = ';'.join(taxonomy)
        # genome taxonomy taxids
        taxids = '|'.join([str(x) for x in genome_idx[genome_name]['taxids']])
        # genome length
        genome_len  = genome_idx[genome_name]['len']
        # database
        db['taxonomy'][taxonomy] = (taxids, genome_len)
    logging.info('  No. of genomes in marker DB: {}'.format(len(db['taxonomy'].keys())))
    
    # marker info
    for clust_id,v1 in cluster_idx.items():
        # genomes that are not clade for the marker
        ext = clade_idx[clust_id]['ext']
        # full taxonomy
        taxon = clade_idx[clust_id]['taxon']
        # clade info
        for clade,v2 in v1.items():
            for seq_uuid,v3 in v2.items():
                gene_len = cluster_idx[clust_id][clade][seq_uuid]['len']
                db['markers'][seq_uuid] = {'clade' : clade[-1],
                                           'ext' : ext,
                                           'len' : gene_len,
                                           'taxon' : taxon}
    logging.info('  No. of markers in marker DB: {}'.format(len(db['markers'].keys())))

    # writing pickle file
    outfile = out_prefix + '.pkl'
    with open(outfile, 'wb') as outF:
        pickle.dump(db, outF)
    logging.info('File written: {}'.format(outfile))
    
def main(args):
    # viewing metaphlan
    if args.metaphlan_pkl is not None:
        view_mphn_pkl(args.metaphlan_pkl, N=10, viruses=False)
    
    # creating indices for all required data
    ## creating taxonomy index {tax_name : taxid}
    tax_idx,tax_graph = create_tax_id_index(args.gtdb_taxdump_names, args.gtdb_taxdump_nodes)    
    ## create membership index {gene_id : cluster_id}
    mem_idx = create_membership_index(args.cluster_membership)    
    ## creating {unirefID : [marker_uuid]} index
    genome_idx,cluster_idx = create_metadata_index(args.gene_metadata, mem_idx, tax_idx)
    mem_idx = None
    if len(cluster_idx) < 1:
        return None
    
    # getting clade-specific clusters
    clade_idx,clade_seq_idx = clade_specific(genome_idx, cluster_idx, tax_graph,
                                             ext_cutoff=args.ext_cutoff)
    if len(clade_idx) < 1:
        return None
    
    # writing output
    ## fasta
    write_fasta(args.gene_fasta, cluster_idx, genome_idx, clade_seq_idx,
                out_prefix = args.output_prefix)
    ## pkl
    write_pkl(genome_idx, cluster_idx, clade_idx,
              out_prefix = args.output_prefix)


# output format: `(NCBI_taxid)__(uniref90_ID)__(CDS_name) (uniref90_ID_full);(taxonomy)`
# >100053__V6HUX8__LEP1GSC062_2642 UniRef90_V6HUX8;k__Bacteria|p__Spirochaetes|c__Spirochaetia|o__Spirochaetia_unclassified|f__Leptospiraceae|g__Leptospira|s__Leptospira_alexanderi;t__GCA_000243815
    
if __name__ == '__main__':
    args = parser.parse_args()
    # potential input
    # gene metadata (unirefID, taxonomy + taxID, genomeID)
    main(args)
