#!/usr/bin/env python
from __future__ import print_function
# batteries
import os
import sys
import re
import argparse
import logging
import csv
import urllib.request
import codecs
from collections import Counter
# 3rd party
import networkx as nx
from networkx.algorithms.dag import descendants
from networkx.algorithms.lowest_common_ancestors import lowest_common_ancestor
from networkx.algorithms.shortest_paths.unweighted import bidirectional_shortest_path

desc = 'Converting NCBI taxonomy to GTDB taxonomy'
epi = """DESCRIPTION:
Using the GTDB metadata table, which contains both NCBI and GTDB taxIDs
to convert NCBI taxonomy to GTDB taxonomy. 

"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('gtdb_metadata', metavar='gtdb_metadata', type=str, nargs='+',
                    help='>=1 gtdb-metadata file (or url)')
parser.add_argument('tax_queries', metavar='tax_queries', type=str,
                    help='List of taxa to query (1 per line)')
parser.add_argument('-q', '--query-taxonomy', type=str, default='ncbi_taxonomy',
                    choices=['ncbi_taxonomy', 'gtdb_taxonomy'],
                    help='Taxonomy of the query list (Default: %(default)s)')
parser.add_argument('-f', '--fraction', type=float, default=0.9,
                    help='Homogeneity of LCA (fraction) in order to be used (Default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def add_taxonomy(line, line_num, header, G, tax='ncbi_taxonomy'):
    """
    Adding taxonomy nodes/edits to the graph
    """
    regex = re.compile(r'^[dpcofgs]__')
    hierarchy = ['domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']
    # checking taxonomy format
    acc = line[header['accession']]
    T = line[header[tax]].split(';')
    T = [regex.sub('', x) for x in T]
    if len(T) < 7:
        msg = 'WARNING: Line{}: "{}" length is <7'
        logging.info(msg.format(line_num, tax))
        for i in range(7):
            try:
                _ = T[i]
            except IndexError:
                T.append('unclassified')
    T.append(acc)
    # adding taxonomy to graph
    for i in range(8):
        # adding node
        G[tax].add_node(T[i], taxonomy=hierarchy[i])
        # adding edge
        if i == 0:
            G[tax].add_edge('root', T[i])
        else:
            G[tax].add_edge(T[i-1], T[i])
        
def load_gtdb_metadata(infile, G):
    """ loading gtdb taxonomy & adding to DAG """
    # input
    try:
        ftpstream = urllib.request.urlopen(infile)
        inF = csv.reader(codecs.iterdecode(ftpstream, 'utf-8'))
    except ValueError:
        inF = open(infile)

    header = {}
    for i,line in enumerate(inF):        
        # parsing
        try:
            line = line.rstrip()
        except AttributeError:
            line = line[0].rstrip()
        if line == '':
            continue
        line = line.split('\t')
        if len(line) < 2:
            msg = 'Line{} does not contain >=2 columns'
            raise ValueError(msg.format(i))
        # header
        if i == 0:
            header = {x:ii for ii,x in enumerate(line)}
            continue
        # taxonomies
        add_taxonomy(line, i, header, G, tax='gtdb_taxonomy')
        add_taxonomy(line, i, header, G, tax='ncbi_taxonomy')
    try:
        inF.close()
    except AttributeError:
        pass
    return G

def DiGraph_w_root():
    """
    directed graph with a root node
    """
    G = nx.DiGraph()
    G.add_node('root')
    return G


def lca_frac_pass(D, lca_frac):
    D = Counter(D)
    mc = D.most_common(1)
    frac = mc[0][1] / float(sum(D.values()))
    if frac >= lca_frac:
        return [mc[0][0], str(round(frac, 3))]
    else:
        return [None,None]

def lca_many_nodes(G, nodes, lca_frac=1.0):
    """
    algorithm: using closest distance to 'root'
    """
    hierarchy = ['root', 'domain', 'phylum', 'class', 'order',
                 'family', 'genus', 'species', 'strain']    
    T = [{} for x in hierarchy]
    # getting path from nodes to root
    for n in nodes:
        path = bidirectional_shortest_path(G, 'root', n)
        for i,node in enumerate(path):
            try:
                T[i][node] += 1
            except KeyError:
                T[i][node] = 1
    # from tip to root, which passess LCA cutoff?
    for i in range(len(T))[::-1]:
        lca = lca_frac_pass(T[i], lca_frac)
        if lca[0] is not None:
            return lca
    raise ValueError('Cannot find LCA for nodes: {}'.format(','.join(nodes)))

def query_tax(tax_queries, G, qtax, ttax, lca_frac=1.0):
    logging.info('Loading queries: {}'.format(tax_queries))
    idx = {}
    status = {'hit' : 0, 'no hit': 0}
    # iterating queries
    with open(tax_queries) as inF:
        for line in inF:
            line = line.rstrip()
            if line == '' or line == 'root':
                continue
            # is node in G[qtax]?
            tips = []
            try:
                # getting descendents of the node
                for desc in list(descendants(G[qtax], line)):
                    if G[qtax].nodes[desc]['taxonomy'] == 'strain':
                        tips.append(desc)
                status['hit'] += 1
            except nx.exception.NetworkXError:
                status['no hit'] += 1
            # if tips, getting LCA in target-taxonomy
            if len(tips) > 0:                
                LCA = lca_many_nodes(G[ttax], tips, lca_frac=lca_frac)
                idx[line] = LCA
            # status
            x = status['hit'] + status['no hit']
            if x % 10 == 0:
                #sys.stderr.write('Queries processed: {}\r'.format(x))
                logging.info('Queries processed: {}'.format(x))
                                
    # status
    msg = 'No. of hits to {}: {}'
    logging.info(msg.format(qtax, status['hit']))
    msg = 'No. lacking hits to {}: {}'
    logging.info(msg.format(qtax, status['no hit']))
    # return
    return idx

def write_table(idx, qtax, ttax):
    print('\t'.join([qtax, ttax, 'lca_frac']))
    for k,v in idx.items():
        print('\t'.join([k] + v))

def main(args):
    # loading the graphs
    G = {'ncbi_taxonomy' : DiGraph_w_root(),
         'gtdb_taxonomy' : DiGraph_w_root()}
    for F in args.gtdb_metadata:
        logging.info('Loading: {}'.format(F))
        load_gtdb_metadata(F, G)
    # querying
    ttax = 'ncbi_taxonomy' if args.query_taxonomy == 'gtdb_taxonomy' else 'gtdb_taxonomy'
    idx = query_tax(args.tax_queries, G, qtax=args.query_taxonomy, ttax=ttax,
                    lca_frac = args.fraction)
    # writing results
    write_table(idx, qtax=args.query_taxonomy, ttax=ttax)
             
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

# # networkx approach
# * for each taxon in the GTDB
#   * create DAG of NCBI taxonomy (taxon at tip, root => tips)
#   * create DAG of GTDB taxonomy (taxon at tip, root => tips)
# * For each query:
#   * find in (NCBI|GTDB) DAG
#     * G.nodes(taxonomy="foo") 
#   * find all extant children
#     * networkx.algorithms.dag.descendants()
#     * DiGraph.successors(n)  # until no more successor
#   * for extant children in other db, find LCA
#     * all_pairs_lowest_common_ancestor()
