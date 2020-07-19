#!/usr/bin/env python
from __future__ import print_function
import sys,os
import argparse
import logging
import tempfile
import csv
import urllib.request
import codecs
from distutils.spawn import find_executable
import subprocess

desc = 'Prune >=1 phylogeny'
epi = """DESCRIPTION:
Prune GTDB phylogeny (bacteria and/ archaea)
to just the list of genome accessions provided.

`nw_prune` from the newick_utils toolset is used
for pruning. 

Trees can be provided as files or urls to files.

If >1 tree is provided, then the trees are merged;
`--root-brlen` determines the brlens to the root.

Output is written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('accs_to_keep', metavar='accs_to_keep', type=str, 
                    help='File of genome accessions to keep on the tree. 1 acc per line')
parser.add_argument('tree_file', metavar='tree_file', type=str, nargs='+',
                    help='>=1 newick file (or url to the file)')
parser.add_argument('-r', '--root-brlen', type=float, default=0.0001,
                    help='Root node branch length (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def read_tree(file_or_url, root_brlen):
    """ Reading in tree file (downloading if url) """
    logging.info('Reading tree: {}'.format(file_or_url))
    try:
        line = urllib.request.urlopen(file_or_url).read().decode('utf-8')
        #inF =   #csv.reader(codecs.iterdecode(ftpstream, 'utf-8'))
    except ValueError:
        line = open(file_or_url).read()
              
    line = line.rstrip().rstrip(';')
    line += '100.0:{}'.format(root_brlen)
    
    return line

def read_trees(tree_files, root_brlen):
    """ reading in >= tree file (or url) """
    trees = []
    for F in tree_files:
        trees.append(read_tree(F, root_brlen))
    trees = '(' + ','.join(trees) + ')'
    trees += '100.0:{}'.format(root_brlen)
    return trees

def prune_tree(tree_file, taxa_to_keep):
    """ Pruning trees via nw_prune """
    cmd = ['nw_prune', '-v', '-f', tree_file, taxa_to_keep]
    cmd = ' '.join(cmd)
    logging.info('CMD: {}'.format(cmd))
    try:
        res = subprocess.run(cmd, check=True, shell=True,
                             stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
 
    res = res.stdout.decode().rstrip()
    print(res)
        
def main(args):
    # checking for newick_utils exe
    if find_executable('nw_prune') is None:
        msg = 'Cannot find "nw_prune" in PATH. Is newick_utils installed?'
        raise IOError(msg)
    
    # downloading/merging trees
    ## temp output file
    dirpath = tempfile.mkdtemp()
    tmpTree_name = os.path.join(dirpath, 'TMP.nwk')
    ## reading in trees
    with open(tmpTree_name, 'w') as tmpTree:        
        if len(args.tree_file) > 1:
            tree = read_trees(args.tree_file, args.root_brlen)
        else:
            tree = read_tree(args.tree_file[0], args.root_brlen)
        tmpTree.write(tree + ';')

    # pruning 
    prune_tree(tmpTree_name, args.accs_to_keep)
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
