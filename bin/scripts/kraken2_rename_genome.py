#!/usr/bin/env python
from __future__ import print_function
import sys,os
import gzip
import argparse
import logging

desc = 'Renaming sequence headers in a genome fasta to kraken2 db format'
epi = """DESCRIPTION:
Output format: `kraken:taxid|<taxID>|<seqID>`
"""

parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('genome_file', metavar='genome_file', type=str,
                    help='genome file (can be gzip\'ed)')
parser.add_argument('taxID', metavar='taxID', type=str,
                    help='taxonomy ID used for renaming sequences')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def main(args):
    if args.genome_file.endswith('.gz'):
        inF = gzip.open(args.genome_file)
    else:
        inF = open(args.genome_file)
    for line in inF:
        line = line.decode('utf8').rstrip()
        if line.startswith('>'):
            line = '>kraken:taxid|{}|{}'.format(args.taxID, line.lstrip('>'))
        print(line)

    inF.close()
    

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
