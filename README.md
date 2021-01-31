![Struo2](https://github.com/leylabmpi/Struo2/workflows/Struo2/badge.svg)

Struo2
======

**Struo2:** a pipeline for building custom databases for common metagenome profilers

> "Struo" --> from the Latin: “I build” or “I gather”

![](https://media.giphy.com/media/lPSMFWqKCJY5Let2v5/giphy.gif)

* Version: 2.1.2
* Authors:
  * Nick Youngblut <nyoungb2@gmail.com>
* Maintainers:
  * Nick Youngblut <nyoungb2@gmail.com>

# Citation

Struo version 2 has major changes from version 1, but verion 2 has not yet been published.
For now, please just cite Struo1:

## Struo version 1 

Cuesta-Zuluaga, Jacobo de la, Ruth E. Ley, and Nicholas D. Youngblut. 2019.
"Struo: A Pipeline for Building Custom Databases for Common Metagenome Profilers."
Bioinformatics , November.
[https://doi.org/10.1093/bioinformatics/btz899](https://doi.org/10.1093/bioinformatics/btz899)

# Overview

Efficiently create/update custom databases for the following metagenome profilers:

* [Kraken2](https://github.com/DerrickWood/kraken2)
* [Bracken](https://github.com/jenniferlu717/Bracken)
* [HUMAnN3](https://github.com/biobakery/humann)

You can also just use Struo2 for efficiently clustering genes via mmseqs and
generating gene & gene-cluster databases that can be efficiently updated via
mmseqs2.

# Changes from Version 1

* Support for HUMAnN3
  * HUMAnN3 uses a much more updated version of UniRef and a new version of DIAMOND
* All coding sequences from all genomes are now clustered prior to annotation,
  and then the annotations (UniRef IDs, by default) are propagated to each member
  of each cluster.
  * This is substantially faster than the per-genome annotation appraoch used for Struo1
  * It also allows for efficient database updates via `mmseqs clusterupdate`.
* By default `mmseqs search` is used for annotation instead of `diamond blastp`
  * `mmseqs search` can be a bit faster and more sensitive than DIAMOND
    * see [Steinegger and Soeding 2017](https://www.nature.com/articles/nbt.3988)
* Metadata is saved for each gene in the database (e.g., gene taxonomy)
* Each database can be updated incrementally
  * Intermediate files are saved for faster database reconstruction
  * For the database of all genes, `mmseqs clusterupdate` is used to add to the
    existing cluster database instead of re-clustering all sequences.
  * For HUMAnN3, only clusters lacking an annotation (to UniRef by default)
    are queried, which can save a great deal of time.
* Users can provide genes via either:
  * A set of genomes fasta files (genes called by prodigal)
  * A set of gene sequences (eg., produced by [PLASS](https://github.com/soedinglab/plass))
    * If only amino-acid genes provided, the sequences can be rev-translated or skipped
    * If only nucleotide genes provided, the sequences can be translated or skipped
* Experimental support for metaphlan3, but see the notes below

# Pre-built custom databases

Custom GTDB databases available at the [struo data ftp server](http://ftp.tue.mpg.de/ebio/projects/struo2/)

**GTDB releases available:**
* Release 95 (13.07.2020)
  * Number of genomes included: 30,989
  * GTDB taxdump
    * taxIDs assigned with [gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump)
  * Genome phylogeny
    * GTDB `ar122_r95.tree` & `bac120_r95.tree` grafted together
  
# Setup for installing Struo2

> Only Unix/Linux OS supported

## Download

To download the pipeline, clone the Git repository:

```
git clone --recurse-submodules git@github.com:leylabmpi/struo2.git 
```

Note the use of submodules. If needed, you can update the submodule(s) via:

```
git submodule update --remote --init --recursive
```

## conda env setup

> Versions listed are those that have been tested. Newer versions will likely work

* If just running the pipeline
  * python=3.6
  * snakemake=5.31.1
* If using the utility scripts
  * r-base=3.6
  * r-argparse=2.0.1
  * r-curl=4.2
  * r-data.table=1.12.4
  * r-dplyr=0.8.3
  * ncbi-genome-download=0.2.10
  * newick_utils=1.6

If you want email notifications upon pipeline success/failure, then you need
mutt installed on your OS.

## Setting a location for necessary files

```
# creating a directory to hold all of the necessary Struo2 data files
OUTDIR=./data/
mkdir -p $OUTDIR
```

## taxdump files

The taxdump files are used for creating/updating the Kraken2 database.

### GTDB taxdump

By default, the pipeline uses custom GTDB taxIDs generated with
[gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump). 
To download the custom taxdump files:

```
wget --directory-prefix $OUTDIR http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/taxdump/names.dmp
wget --directory-prefix $OUTDIR http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/taxdump/nodes.dmp
```

### NCBI taxdump

If you would rather use NCBI taxonomy instead of the GTDB taxonomy:

```
wget --directory-prefix $OUTDIR https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -pzxvf $OUTDIR/taxdump.tar.gz
```

## UniRef

UniRef databases are used for annotating genes.
UniRef IDs are required for HUMANnN3.
You do not need any UniRef databases if just creating Kraken2/Bracken databases.

### IF using `mmseqs search` for gene annotation

**mmseqs UniRef database(s)**

See the [mmseqs2 wiki on database downloading](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases)

```
# you must have mmseqs2 installed
# Example of downloading UniRef50
mmseqs databases UniRef50 $OUTDIR/UniRef50 DB_TMP
```

### IF using `diamond blastp` for gene annotation

**HUMAnN3 UniRef diamond database(s)**

See the "Download a translated search database" section of the
[humann3 docs](https://github.com/biobakery/biobakery/wiki/humann3#welcome-to-the-humann-30-tutorial).

```
# Example download of UniRef50 DIAMOND database
wget --directory-prefix $OUTDIR http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901.tar.gz
```

## UniRef50-90 index

> Optional, but recommended

This is needed to map annotations from UniRef90 clusters to UniRef50 clusters.
This allows for just annotating against UniRef90 and then mapping
those annotations to UniRef50 cluster IDs.
You then do not have to annotate against UniRef90 clusters and UniRef50 clusters,
which requires a lot more querying of genes against UniRef.

```
wget --directory-prefix $OUTDIR http://ftp.tue.mpg.de/ebio/projects/struo2/install/uniref_2019.01/uniref50-90.pkl
```

## Getting reference genomes for the custom databases

### Downloading genomes

* If using [GTDB](https://gtdb.ecogenomic.org/) genomes, run `GTDB_metadata_filter.R` to select genomes
* If downloading genomes from genbank/refseq, you can use `genome_download.R`
* You can also include your own genomes (e.g., newly created MAGs)

**An example:**

```
# Filtering GTDB metadata to certain genomes
./GTDB_metadata_filter.R -o gtdb-r95_bac-arc.tsv https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_metadata_r95.tar.gz

# Downloading all genomes (& creating tab-delim table of genome info)
./genome_download.R -o genomes -p 8 gtdb-r95_bac-arc.tsv > genomes.txt

# Note: the output of ./genome_download.R can be directly used for running the `Struo2` pipeline (see below)
# Note: genome fasta files can be compressed (gzip or bzip2) or uncompress for input to Struo2
```


# Input genome/gene data

The genomes and/or individual genes used for creating/updating databases

## Genome data

The table of input files/data can be created using the helper scripts described above
if you downloaded the genomes from the GTDB or NCBI.

The `samples.txt` file, by default.

* The pipeline requires a tab-delimited table that includes the following columns (column names specified in the `config.yaml` file):
  * `samples_col` (default = `ncbi_organism_name`)
    * This will usually just be the species/strain names
    * Do not include special characters
  * `accession_col` (default = `accession`)
    * Genome accession
    * You can fill with "blank" or "NA" if no accession is available
  * `fasta_file_path_col` (default = `fasta_file_path`)
    * Path to the genome fasta
    * The fasta can be uncompressed or compressed (gzip or bzip2)
  * `taxID_col` (default = `gtdb_taxid`)    
    * This should be the NCBI/GTDB taxonomy ID at the species/strain level
      * Needed for Kraken
    * For custom genomes (e.g., MAGs), you will need to taxonomically classify the genomes
      and get a taxid 
  * `taxonomy_col` (default = `gtdb_taxonomy`)
    * This should at least include `g__<genus>;s__<species>`
    * The taxonomy can include higher levels, as long as levels 6 & 7 are genus and species
    * Any taxonomy lacking genus and/or species levels will be labeled:
      * `g__unclassified`  (if no genus)
      * `s__unclassified`  (if no species)
    * This is needed for HUMAnN3

Other columns in the file will be ignored. The path to the samples file should be specified in the `config.yaml` file (see below)

## Gene data

This can only be used for updating existing databases
(e.g., existing custom GTDB Kraken2 or HUMAnN databases).

The gene data must include:

* Amino acid sequencing in fasta format
  * [optional] corresponding nucleotide sequences in fasta format
  * The fasta can be compressed via gzip
* A tab-delimited table of metadata with the following columns:
  * `seq_uuid`
    * A unique identifier for that sequence
  * `seq_orig_name`
    * The original name for that sequence
  * `genus`
    * The genus-level taxonomy of that sequences
    * Format: `g__<taxonomy>`
    * If unknown, use `g__unclassified` 
  * `species`
    * The species-level taxonomy of that sequences
    * Format: `s__<taxonomy>`
    * If unknown, use `s__unclassified`
  * `taxid`
    * The taxid corresponding to the genus-species taxonomy
    * [taxonkit](https://bioinf.shenwei.me/taxonkit/) can help with this


# Running the pipeline

> This applies to creating and updating the databases

## Edit the config file

* Select a config to edit
  * If creating a new database: edit `config.yaml`
  * If updating an existing database: edit `config-update.yaml`
* Specify the input/output files paths
  * See above for downloading input files
* Modify parameters as needed
  * Make sure to add the path to the UniRef diamond database for HUMAnN3
    * see above for instructions on retrieving this file
  * Can the taxdump files from GTDB to NCBI if needed
  * By defeault, `mmseqs search` is used instead of `diamond blastp`
* Modify `temp_folder:` if needed
  * This folder is used just for read/write of temporary files
  * Use SSDs, if available 


### Running locally

This is only recommended for test runs

It is always good to run snakemake in `dryrun` mode first:

```
snakemake --use-conda -j -Fqn
```

For an actual run with 2 cores:

```
snakemake --use-conda -j 2
```

### Running on a cluster

If SGE, then you can use the `snakemake_sge.sh` script. You can create a similar bash script
for other cluster architectures. See the following resources for help:

* [Ley Lab snakemake profiles](https://github.com/leylabmpi/snakemake_profiles)
* [Snakemake docs on cluster config](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)
* [Official snakemake profiles](https://github.com/Snakemake-Profiles)

### General info on using `snakemake`

Snakemake allows for easy re-running of the pipeline on just genomes that have not yet been processed.
You can just add more genomes to the input table and re-run the pipeline (test first with `--dryrun`).
Snakemake should just process the new genomes and then re-create the combined dataset files (this must be done each time).
Make sure to not mess with the files in the `nuc_filtered` and `prot_filtered` directories! Otherwise,
snakemake may try to run all genomes again through the computationally expensive gene annotation process.


## Using the resulting databases

Set the database paths in humann3, kraken2, etc. to the new, custom database files.

* genes
  * `genome_reps_filtered.fna.gz`
    * representative gene sequences (nucleotide)
  * `genome_reps_filtered.faa.gz`
    * representative gene sequences (amino acid)
  * `genome_reps_filtered.tsv.gz`
    * representative gene sequences (metadata)
  * `genes_db.tar.gz`
    * mmseqs genes database
  * clusters
    * `clusters_db.tar.gz`
      * mmseqs cluster database
    * `clusters_membership.tsv.gz`
      * mmseqs cluster membership table
    * `clusters_reps.faa.gz`
      * cluster representative sequences (amino acid)
* kraken2
  * `*.k2d`
* bracken
  * `database*mers.kmer_distrib`
* humann3
  * nucleotide
    * `genome_reps_filt_annot.fna.gz`
      * plus the `*.bt2` files 
  * amino acid
    * `uniref*_*.dmnd`
      * Example: `uniref90_201901.dmnd`
  
### Example of a HUMANnN3 run

Run HUMAnN3 with custom databases created by Struo2. Change that PATHs as necessary. 

```
STRUO2_OUT_DIR=./struo2_output/
NUC_DB=$STRUO2_OUT_DIR/humann3/uniref50/genome_reps_filt_annot.fna.gz
PROT_DB=$STRUO2_OUT_DIR/humann3/uniref50/protein_database/uniref50_201901.dmnd
MP_DB=$STRUO2_OUT_DIR/mpa_v30_CHOCOPhlAn_201901_marker_info.txt

READS=/path/to/example/read/files/CHANGE/THIS/reads.fq

# metaphlan database not actually used
touch $MP_DB

# humann3 run	
humann3 --bypass-nucleotide-index \
  --nucleotide-database $NUC_DB \
  --protein-database $PROT_DB \
  --metaphlan-options "--bowtie2db $MP_DB" \
  --output-basename humann \
  --input $READS \
  --output humann_output
```


# Utilities

## `GTDB_metadata_filter.R`

This tool is useful for selecting which GTDB genomes to include in a custom database.

Filter >=1 genome assembly metadata file (e.g., bac120_metadata_r89.tsv)
by assembly quality or other parameters.

## `genome_download.R`

This tool is useful for downloading genomes from NCBI. 

Download a set of genomes based on NCBI assembly accessions provided
in a table. The file paths of the downloaded genome fasta files will be
appended to the input table.

## `tree_prune.py`

This tool is useful for creating a GTDB archaea/bacteria phylogeny
of all genomes in your custom database. The phylogeny can be used
for phylogenetic analyses of metagenomes (e.g., Faith's PD or Unifrac).

Prune >=1 phylogeny to just certain taxa. If >1 phylogeny provided,
then the phylogenies are merged.

## `gtdb_to_taxdump`

This is a [separate repo](https://github.com/nick-youngblut/gtdb_to_taxdump).

This is useful for creating an NCBI taxdump (names.dmp and nodes.dmp)
from the GTDB taxonomy. Note that the taxIDs are arbitrary and don't
match anything in the NCBI! 

## `nbci-gtdb_map.py`

This is a [separate repo](https://github.com/nick-youngblut/gtdb_to_taxdump).

This is useful for mapping between GTDB and NCBI taxonomies.
The mapping is based on the GTDB archaeal and bacterial metadata tables,
which contain both GTDB and NCBI taxonomic lineages. 

