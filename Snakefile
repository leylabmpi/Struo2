# import 
from __future__ import print_function
import os
import sys
import re
import glob
import socket
import getpass
import subprocess
from distutils.spawn import find_executable
import pandas as pd

# setup
def skipped(x):
    return str(x).strip().lower().startswith('skip')

def make_fasta_splits(n_jobs):
    if str(n_jobs).lstrip().startswith('Skip'):
        n_jobs = 1
    zero_pad = len(str(n_jobs))
    zero_pad = '{0:0' + str(zero_pad) + 'd}'
    return [str(zero_pad.format(x+1)) for x in range(n_jobs)]

## load
configfile: 'config.yaml'

## outdir
config['output_dir'] = config['output_dir'].rstrip('/') + '/'
config['db_name'] = config['db_name'].rstrip('/') + '/'

## Samples table
if not os.path.isfile(config['samples_file']):
    raise IOError('Cannot find file: {}'.format(config['samples_file']))
config['samples'] = pd.read_csv(config['samples_file'], sep='\t')
for f in [config['samples_col'], config['fasta_file_path_col'],
          config['taxID_col'], config['taxonomy_col']]:
    if f not in config['samples'].columns:
        raise ValueError('Cannot find column: {}'.format(f))
config['samples'][config['samples_col']] = config['samples'][config['samples_col']].str.replace('[^A-Za-z0-9]+', '_')
config['samples'] = config['samples'].set_index(config['samples'][config['samples_col']])

### check that files exist (skipping if not)
rowID = 0
to_rm = []
for index,row in config['samples'].iterrows():
    rowID += 1
    file_cols = [config['fasta_file_path_col']]
    for f in file_cols:
        if not os.path.isfile(str(row[f])):
           msg = 'Samples table (Row {}): Cannot find file: {}; Skipping\n'
           sys.stderr.write(msg.format(rowID, row[f]))
           to_rm.append(row[config['samples_col']])
sys.stderr.write('Total number of skipped rows: {}\n'.format(len(to_rm)))
config['samples'].drop(to_rm, inplace=True)
if config['samples'].shape[0] < 1:
    raise ValueError('No genomes remaining after filtering!')
config['samples_unique'] = config['samples'][config['samples_col']].unique().tolist()

## temp_folder
config['pipeline']['username'] = getpass.getuser()
config['pipeline']['email'] = config['pipeline']['username'] + '@tuebingen.mpg.de'
config['tmp_dir'] = os.path.join(config['tmp_dir'], config['pipeline']['username'])
config['tmp_dir'] = os.path.join(config['tmp_dir'], 'Struo2_' + str(os.stat('.').st_ino) + '/')
print('\33[33mUsing temporary directory: {} \x1b[0m'.format(config['tmp_dir']))

## batches
config['params']['humann3']['splits'] = \
  make_fasta_splits(config['params']['humann3']['batches'])

## including modular snakefiles
snake_dir = config['pipeline']['snakemake_folder']
include: snake_dir + 'bin/dirs'
include: snake_dir + 'bin/Snakefile'
if not skipped(config['databases']['genes']):
    include: snake_dir + 'bin/genes/Snakefile'
else:
    print('\33[33mSkipping creation of genes database; assuming the database already exist!\x1b[0m')
if not skipped(config['databases']['kraken2']):
    include: snake_dir + 'bin/kraken2/Snakefile'
    if not skipped(config['databases']['bracken']):
        include: snake_dir + 'bin/bracken/Snakefile'
if not skipped(config['databases']['humann3_bowtie2']) and \
   not skipped(config['databases']['humann3_diamond']):
    include: snake_dir + 'bin/humann3/query_dmnd/Snakefile'
    include: snake_dir + 'bin/humann3/query_mmseqs/Snakefile'
    include: snake_dir + 'bin/humann3/db_create/Snakefile'
if not skipped(config['databases']['metaphlan3']):
    include: snake_dir + 'bin/metaphlan3/Snakefile'


## pipeline main
wildcard_constraints:
    sample="[^/]+"

localrules: all

rule all:
    input:
        all_which_input

