Struo2
======

**Struo2:** a pipeline for building custom databases for common metagenome profilers

> "Struo" --> from the Latin: “I build” or “I gather”


* Version: 2.0.1
* Authors:
  * Nick Youngblut <nyoungb2@gmail.com>
  * Jacobo de la Cuesta <jacobo.delacuesta@tuebingen.mpg.de>
* Maintainers:
  * Nick Youngblut <nyoungb2@gmail.com>
  * Jacobo de la Cuesta <jacobo.delacuesta@tuebingen.mpg.de>


# Citation

Struo version 2 has major changes from version 1, but verion 2 has not yet been published.

## Struo version 1 

Cuesta-Zuluaga, Jacobo de la, Ruth E. Ley, and Nicholas D. Youngblut. 2019.
"Struo: A Pipeline for Building Custom Databases for Common Metagenome Profilers."
Bioinformatics , November.
[https://doi.org/10.1093/bioinformatics/btz899](https://doi.org/10.1093/bioinformatics/btz899)

# Changes from Version 1

* All coding sequences from all genomes are now clustered prior to annotation
(DIAMOND hits to UniRef), and then the UniRef IDs are propagated to each member
of each cluster. This can greatly cut down on the amount of DIAMOND computation.

# Pre-built custom databases

## Version 1

Custom GTDB databases available at the [struo data ftp server](http://ftp.tue.mpg.de/ebio/projects/struo/)

**GTDB releases available:**
* Release 86 (14.03.2019)
  * Number of genomes included: 21,276
  * NCBI taxonomy/taxIDs used
* Release 89 (30.08.2019)
  * Number of genomes included: 23,361
  * GTDB taxonomy/taxIDs used
    * taxIDs assigned with [gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump)
  
# Tutorial

For a step-by-step example of how to prepare and execute Struo, see the notebook in the `./tutorial/` folder

# Description

## Struo’s workflow

![](./images/struo_workflow.png)
Struo's workflow encompasses the steps from genome download to database construction

## Setup

### Download

To download the pipeline, clone the Git repository:

```
git clone git@github.com:leylabmpi/Struo.git 
```

### conda env setup

> Versions listed are those that have been tested

* python=3.6
* snakemake=5.7.0
* r-base=3.6
* r-argparse=2.0.1
* r-curl=4.2
* r-data.table=1.12.4
* r-dplyr=0.8.3
* ncbi-genome-download=0.2.10
* newick_utils=1.6

### UniRef diamond database(s)

You will need a UniRef diamond database for the humann3 database construction (e.g., UniRef90).
See the "Download a translated search database" section of the
[humann3 docs](https://github.com/biobakery/biobakery/wiki/humann3#welcome-to-the-humann-30-tutorial)

## Getting reference genomes for the custom databases

### Downloading genomes

* If using [GTDB](https://gtdb.ecogenomic.org/) genomes, run `GTDB_metadata_filter.R` to select genomes
* If downloading genomes from genbank/refseq, you can use `genome_download.R`

Example:

```
# Filtering GTDB metadata to certain genomes
./GTDB_metadata_filter.R -o gtdb-r89_bac-arc.tsv https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/ar122_metadata_r89.tsv

# Downloading all genomes (& creating tab-delim table of genome info)
./genome_download.R -o genomes -p 8 gtdb-r89_bac-arc.tsv > genomes.txt

# Note: the output of ./genome_download.R can be directly used for running the `Struo` pipeline (see below)
```

### User-provided databases

Users can also provide genomes as compressed fasta files (`.fna.gz`).
This also requires adding the corresponding information to the `samples.txt` file (see below)

## Input data (`samples.txt` file)

The table of input files/data can be created using the helper scripts described above. 

* The pipeline requires a tab-delimited table that includes the following columns (column names specified in the `config.yaml` file):
  * Sample ID
    * This will usually just be the species/strain names
  * Path to the genome assembly fasta file
    * NOTE: these must be gzip'ed
  * taxonomy ID
    * This should be the NCBI taxonomy ID at the species/strain level
      * Needed for Kraken
  * taxonomy
    * This should at least include `g__<genus>;s__<species>`
    * The taxonomy can include higher levels, as long as levels 6 & 7 are genus and species
    * Any taxonomy lacking genus and/or species levels will be labeled:
      * `g__unclassified`  (if no genus)
      * `s__unclassified`  (if no species)
    * This is needed for humann3

Other columns in the file will be ignored. The path to the samples file should be specified in the `config.yaml` file (see below)

### Using the GTDB taxonomy instead of NCBI taxIDs

kraken2 & humann3 databases used NCBI taxIDs, and thus the NCBI taxonomy is used by default
for `Struo`. You can instead create custom taxIDs from the GTDB taxonomy with
[gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump). 

The resulting `names.dmp` and `nodes.dmp` files, along with a genome metadata file that includes the gtdb_taxids,
then you can modify the Struo pipeline to fully use the GTDB taxonomy & taxIDs.
You will need to modify the `config.yaml` file (see "If using GTDB taxIDs" below).


## Running the pipeline

### Edit the `config.yaml`

* Specify the input/output paths
* Modify parameters as needed
  * Make sure to add the path to the UniRef diamond database for HUMAnN3
    * see above for instructions on retrieving this file
* Modify `temp_folder:` if needed
  * This folder is used just for read/write of temporary files

#### If using GTDB taxIDs

If you have followed "Using the GTDB taxonomy instead of NCBI taxIDs" above, then
make the following modifications to the `config.yaml` file:

```
## column names in samples table
taxID_col: 'gtdb_taxid'
taxonomy_col: 'gtdb_taxonomy'

#-- if custom NCBI taxdump files --#
names_dmp: /YOUR/PATH/TO/names.dmp
nodes_dmp: /YOUR/PATH/TO/nodes.dmp
```


### Running locally

`snakemake --use-conda`

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

Set the database paths in humann2, kraken2, etc. to the new, custom database files.

* humann2
  * nucleotide
    * `all_genes_annot.fna.gz`
  * amino acid
    * `all_genes_annot.dmnd`
* kraken2
  * `database*mers.kraken`
  
### Example of a humann2 run

Run humann2 with custom databases created by Struo. Change that PATHs as necessary. 

```
STRUO_OUT_DIR=./struo_output/
NUC_DB=`dirname $STRUO_OUT_DIR"/all_genes_annot.fna.gz"`
PROT_DB=`dirname $STRUO_OUT_DIR"/all_genes_annot.dmnd"`
MTPHLN_BT2_DB=`dirname ./metaphlan2_db/mpa_v20_m200/mpa_v20_m200.1.bt2`
MTPHLN_PKL_DB=/ebio/abt3_projects/databases_no-backup/metaphlan2/mpa_v20_m200/mpa_v20_m200.pkl

humann2 --gap-fill on --bypass-nucleotide-index  \
  --nucleotide-database $NUC_DB  \
  --protein-database $PROT_DB \
  --metaphlan-options "Skip --mpa_pkl $MTPHLN_PKL_DB --bowtie2db $MTPHLN_BT2_DB" \
  --tmp-dir /dev/shm/humann2_temp/ \
  --threads 12 \
  --input-format fastq  \
  --output-basename SRS018656 \
  --input SRS018656_R1.fq
```
  

### Adding more samples (genomes) to an existing custom DB

If you set `keep_intermediate: True` for your initial run, then the
intermediate files from the computationally intensive steps are kept,
and so those genomes don't have to be reprocessed. Only new genomes will
be processed, and then the database(s) will be re-created with old + new
genomes.

To create a database with more genomes:

* Add new genomes to the input table.
* **If** you want to over-write your old databases:
  * DO NOT change the `db_name:` parameter in the config.yaml file
* **OR if** you want to create new database:
  * Change the `db_name:` parameter in the config.yaml file
* Re-run the snakemake pipeline.
  * Snakemake should skip the genomes that have already been processed.
  * Use `--dryrun` to see what snakemake is going to do before actually running the pipeline.
  * You may need to set `use_ancient: True` in order to have snakemake skip the diamond mapping for humann2
    * This is needed if the timestamps on the genome gene files have been (accidently) modified since the last run.


### Adding existing gene sequences to humann2 databases

If you have gene sequences already formatted for creating a humann2 custom DB,
and you'd like to include them with the gene sequences generated from the
input genomes, then just provide the file paths to the nuc/prot fasta files
(`humann2_nuc_seqs` and `humann2_prot_seqs` in the `config.yaml` file).

All genes (from genomes & user-provided) will be clustered altogether with `vsearch`.
See the `vsearch_all:` setting in the `config.yaml` for the default clustering parameters used.
You can use `vsearch_all: Skip` to skip the clustering and instead all of the sequences
will just be combined without removing redundancies.


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



# TODO

* Create a diamond DB using `diamond >=0.9` so that users can run humann2 with
the most up-to-date version of diamond
  * Note this will require creating an updated UniRef50 db
* Add support for humann3
* Create databases for GTDBr90
* Create [metaphlan3 marker database](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0#customizing-the-database)
  * [issue about creating the DB](https://github.com/biobakery/MetaPhlAn/issues/103)
  * Get taxonomy of the genomes (provided by GTDB in the metadata)
  * Call genes via prokka
  * Map to UniRef90 via DIAMOND (--evalue=1, --query-cover=80, --id=90)
  * For each species in the genome taxonomy:
    * determine which markers are shared by all members of the species (same uniref ID)
    * for markers that are just 'nearly' unique to a species:
      * `ext` => all genomes "external" to the clade
      * rule of thumb: <=10 "external" genomes
  * Rename sequences:
    * marker_name = `(NCBI_taxid)(UniRef90_cluster)(CDS_name)`
    * Example `>100053__V6HZP8__LEP1GSC062_3436`
  * Determine marker info
    * clade, ext, len, taxon
      * ext: list of genomes sharing that marker
      * clade = species
      * taxon = species?
    * ...also the taxonomy of each genome
    * ...then update the mpa pkl file
  * Create bowtie2 database
   
   