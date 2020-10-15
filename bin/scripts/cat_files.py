#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import gzip
import argparse
import logging

desc = 'Simple cat of files; files can be gzip\'d'
epi = """DESCRIPTION:
Simple concatentation of files that allows for a mixture
of gzip'd and uncompressed files.
Output written to STDOUT
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('input', metavar='input', type=str, nargs='+',
                    help='Input file(s)')
parser.add_argument('--header', action='store_true', default=False,
                    help='Input files have headers, so just keep first (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')


def main(args):    
    for i,infile in enumerate(args.input):
        if infile.endswith('.gz'):
            _open = lambda x: gzip.open(x, 'rb')
        else:
            _open = lambda x: open(x, 'r')
        with _open(infile) as inF:
            for ii,line in enumerate(inF):
                # skipping header (except for first table)
                if i > 0 and ii == 0 and args.header is True:
                    continue
                # writing line 
                if infile.endswith('.gz'):
                    line = line.decode('utf-8')
                print(line, end='')
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
