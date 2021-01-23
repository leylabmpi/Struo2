#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import tarfile
import argparse
import logging

desc = 'Uncompress tarball'
epi = """DESCRIPTION:
Simple script for smartly uncompressing a tarball.
All files extracted the same output directory, regardless
of the directory structure in the tarball. 
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('tarball', metavar='tarball', type=str,
                    help='tarball file to extract')
parser.add_argument('-o', '--output-directory', type=str, default='.',
                    help='Output directory location (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def main(args):
    # output location
    if not os.path.isdir(args.output_directory):
        os.makedirs(args.output_directory)
    # extracting
    ext = 'r:gz' if args.tarball.endswith('.gz') else 'r'
    with tarfile.open(args.tarball, ext) as inF:
        files = {k:v for k,v in zip(inF.getnames(), inF.getmembers())}
        for F,M in files.items():
            outfile = os.path.split(F)[1]
            logging.info('Extracting file: {}'.format(outfile))
            outfile = os.path.join(args.output_directory, outfile)
            with inF.extractfile(M) as inFa, open(outfile, 'w') as outF:
                for line in inFa:
                    outF.write(line.decode('utf-8'))
        
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
