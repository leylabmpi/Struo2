#!/usr/bin/env python
from __future__ import print_function
# batteries
import os
import sys
import re
import gzip
import bz2
import random
import argparse
import logging
from functools import partial
from multiprocessing import Pool

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

# argparse
desc = 'A simple simulation of genome misassemblies'
epi = """DESCRIPTION:
* Input is a table of genomes
  * tab-delimited table
  * columns must include:
    * "taxon" => name of taxon (change via --taxon-column)
    * "fasta" => path to assembly fasta file (change via --fasta-column)
* Simulated miassemblies
  * Contig breaks => slicing contigs into mulitiple
  * Relocations => intra genome fragment relocation
  * Chimeras => inter genome fragment relocation  
* Output genome contig naming
  * "_b" = break
  * "_r" = relocation
  * "_c" = chimera
  * "_R" = removal
* Table of [taxon,fasta] written to STDOUT
  * fasta = fasta file path for edited genomes
"""

parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('input_table', metavar='input_table', type=str,
                    help='Tab-delim table of genomes (see Description)')
parser.add_argument('-b', '--num-breaks', type=int, default=2,
                    help='No. of contig breaks to add per genome')
parser.add_argument('-r', '--num-relocates', type=int, default=2,
                    help='No. of fragment relocations within the same genome')
parser.add_argument('-c', '--num-chimeras', type=int, default=2,
                    help='No. of chimeras to add per genome')
parser.add_argument('-o', '--outdir', type=str, default='mis-asmbl',
                    help='Output directory')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Parallel processes')
parser.add_argument('-T', '--taxon-column', type=str, default='taxon',
                    help='Name of taxon column')
parser.add_argument('-F', '--fasta-column', type=str, default='fasta',
                    help='Name of fasta column')
parser.add_argument('-g', '--gzip-out', action='store_true', default=False,
                    help='gzip output?')
parser.add_argument('--version', action='version', version='0.0.1')


def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(x):
    """
    Decoding input, if needed
    """
    try:
        x = x.decode('utf-8')
    except AttributeError:
        pass
    return x

def read_fasta(filename):
    fasta = dict()
    header = None
    with _open(filename) as inF:
        for line in inF:            
            line = _decode(line).rstrip()
            if line == '':
                continue
            if line.startswith('>'):
                header = line.lstrip('>')
                continue
            try:
                fasta[header] += line
            except KeyError:
                fasta[header] = line
    return fasta

def load_input_table(filename, args):
    header = dict()
    D = dict()
    with open(filename) as inF:
        for i,line in enumerate(inF):
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                continue
            taxon = line[header[args.taxon_column]]
            fasta = line[header[args.fasta_column]]
            D[taxon] = fasta    
    return D

def insert_seq(str1, str2, at):
    """ (str1, str2, int) -> str

    Return the DNA sequence obtained by inserting the 
    second DNA sequence into the first DNA sequence 
    at the given index.
    """
    str1_split1 = str1[:at]
    str1_split2 = str1[at:]
    return str1_split1 + str2 + str1_split2

def remove_seq(seq, s_start, s_len):
    seq1 = seq[:s_start]
    seq2 = seq[s_start+s_len:]
    return seq1 + seq2

def select_donor_contig(genome, max_tries, min_len):
    for i in range(max_tries):
        contig_d = random.sample([x for x in genome.keys()], 1)[0]
        contig_d_seq = genome[contig_d]
        if len(contig_d_seq) < min_len:
            continue
        return contig_d,contig_d_seq
    return None,None

def add_break(focal_gen):
    """
    adding break point to contigs
    """
    # selecting contig
    contig = random.sample([x for x in focal_gen.keys()], 1)[0]
    contig_seq = focal_gen[contig]
    # selecting split location
    if len(contig_seq) < 100:
        focal_gen.pop(contig, None)
        return 0
    c_insert = random.randint(1, len(contig_seq)-1)
    logging.info('  BR: Insert location: {}'.format(c_insert))
    # splitting
    frag1 = contig_seq[c_insert:]
    frag2 = contig_seq[:c_insert]
    ## inserting
    focal_gen.pop(contig, None)
    focal_gen[contig + '_b1'] = frag1
    focal_gen[contig + '_b2'] = frag2

def add_relocate(focal_gen, min_len=1000, max_len=10000, max_tries=10):
    """
    adding break point to contigs
    """
    # selecting donor contig for insertion
    contig_d,contig_d_seq = select_donor_contig(focal_gen, max_tries, min_len)
    if contig_d is None:
        return 0
    # selecting location of fragment 
    if len(contig_d_seq) < min_len:
        c_start = 0
        c_len = len(contig_d_seq)
    else:
        max_end = len(contig_d_seq)-min_len
        if max_end < 1:
            c_start = 0
        else:
            c_start = random.randint(0, max_end) 
        c_len = random.randint(min_len, max_len)
    frag_d = contig_d_seq[c_start:c_start+c_len]
    logging.info('  RE: Insert fragment size: {}'.format(len(frag_d)))
    # removing from donor
    focal_gen.pop(contig_d, None)
    focal_gen[contig_d + '_R'] = remove_seq(contig_d_seq, c_start, c_len)
    if len(focal_gen[contig_d + '_R']) < 100:
        focal_gen.pop(contig_d + '_R', None)
    # adding fragment to focal genome
    ## selecting recipient contig
    contig_r = random.sample(focal_gen.keys(), 1)[0]
    contig_r_seq = focal_gen[contig_r]
    ## selecting insert location
    if len(contig_r_seq) < 100:
        focal_gen.pop(contig_r, None)
        return 0
    c_insert_loc = random.randint(0, len(contig_r_seq)-1)
    logging.info('  RE: Insert location: [{}, {}]'.format(contig_r, c_insert_loc+1))
    ## inserting
    focal_gen.pop(contig_r, None)
    focal_gen[contig_r + '_r'] = insert_seq(contig_r_seq, frag_d, c_insert_loc)

def add_chimera(focal_gen, taxon, genomes, min_len=1000, max_len=10000, max_tries=10):
    """
    adding chimeric fragments
    """
    # selecting chimeric genome
    x = set(genomes.keys()) - set(taxon)
    donor = random.sample(list(x), 1)[0]
    logging.info('  CH: Inserting fragment from taxon: {}'.format(donor))
    ## loading genome
    donor_gen = read_fasta(genomes[donor])
    # selecting contig
    # selecting donor contig for insertion
    contig_d,contig_d_seq = select_donor_contig(donor_gen, max_tries, min_len)
    if contig_d is None:
        return 0
    # selecting location
    if len(contig_d_seq) < min_len:
        c_start = 0
        c_len = len(contig_d_seq)
    else:
        max_end = len(contig_d_seq)-min_len
        if max_end < 1:
            c_start = 0
        else:
            c_start = random.randint(0, max_end)        
        c_len = random.randint(min_len, max_len)
    frag_d = contig_d_seq[c_start:c_start+c_len]
    logging.info('    CH: Insert fragment size: {}'.format(len(frag_d)))
    # removing fragment from donor
    donor_gen.pop(contig_d, None)
    donor_gen[contig_d + '_R'] = remove_seq(contig_d_seq, c_start, c_len)
    if len(donor_gen[contig_d + '_R']) < 100:
        donor_gen.pop(contig_d + '_R', None)
    # adding fragment to focal genome
    ## selecting contig
    contig_r = random.sample(focal_gen.keys(), 1)[0]
    contig_r_seq = focal_gen[contig_r]
    ## selecting insert location
    c_insert_loc = random.randint(0, len(contig_r_seq)-1)
    logging.info('    CH: Insert location: [{}, {}]'.format(contig_r, c_insert_loc+1))
    ## inserting
    focal_gen.pop(contig_r, None)
    focal_gen[contig_r + '_c'] = insert_seq(contig_r_seq, frag_d, c_insert_loc)
    
def per_genome(taxon, genomes, num_chimeras, num_breaks, num_relocates):
    """
    modifying genome
    """
    # load focal genome
    logging.info('Processing taxon: {}'.format(taxon))
    focal_gen = read_fasta(genomes[taxon])
    # break points
    for i in range(num_breaks):
        add_break(focal_gen)
    # relocations
    for i in range(num_relocates):
        add_relocate(focal_gen)
    # chimeras
    for i in range(num_chimeras):
        add_chimera(focal_gen, taxon, genomes)
    # result
    return [taxon,focal_gen]

def write_genomes(genomes, outdir, gzip_out=False):
    """
    writing out edited genomes
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfiles = dict()
    for taxon,seqs in genomes:
        outfile = os.path.join(outdir,
                               taxon.replace(' ', '_') + '.fna')
        if gzip_out is True:
            outfile += '.gz'
            with gzip.open(outfile, 'wb') as outF:
                for contig,seq in seqs.items():
                    if len(seq) < 100:
                        continue
                    x = '>' + contig + '\n' + seq + '\n'
                    outF.write(x.encode())
        else:
            with open(outfile, 'w') as outF:
                for contig,seq in seqs.items():
                    if len(seq) < 100:
                        continue
                    outF.write('>' + contig + '\n' + seq + '\n')
        logging.info('File written: {}'.format(outfile))
        outfiles[taxon] = outfile
    return outfiles

def write_table(outfiles):
    """
    writing table of output files
    """
    print('\t'.join(['taxon', 'fasta']))
    for taxon,filename in outfiles.items():
        print('\t'.join([taxon,filename]))

def main(args):
    # load table of genomes
    genomes = load_input_table(args.input_table, args)
    # per genome
    func = partial(per_genome,
                   genomes=genomes,
                   num_chimeras=args.num_chimeras,
                   num_breaks=args.num_breaks,
                   num_relocates=args.num_relocates)
    if args.threads > 1:
        with Pool(args.threads) as p:
            genomes = p.map(func, genomes.keys())
    else:
        genomes = map(func, genomes.keys())
    # writing out genomes
    outfiles = write_genomes(genomes, args.outdir, args.gzip_out)
    # writing table
    write_table(outfiles)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
