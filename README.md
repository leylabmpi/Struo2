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

## UniRef diamond database(s)

You will need a UniRef diamond database for the humann3 database construction (e.g., UniRef50).
See the "Download a translated search database" section of the
[humann3 docs](https://github.com/biobakery/biobakery/wiki/humann3#welcome-to-the-humann-30-tutorial)

## UniRef50-90 index

Needed to map annotations from UniRef90 clusters to UniRef50 clusters.
This allows for just annotating against UniRef90 and then UniRef50 cluster IDs
an be mapped based on the UniRef90 cluster IDs.
You then do not have to annotate against UniRef90 clusters and UniRef50 clusters,
which requires a lot more querying of genes against UniRef. 

```
wget http://ftp.tue.mpg.de/ebio/projects/struo2/install/uniref_2019.01/uniref50-90.pkl
```
## GTDB taxdump

The taxdump files are used for creating/updating the Kraken2 database.

By default, the pipeline uses custom GTDB taxIDs generated with
[gtdb_to_taxdump](https://github.com/nick-youngblut/gtdb_to_taxdump). 
To download the custom taxdump files:

```
wget http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/taxdump/names.dmp
wget http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/taxdump/nodes.dmp
```

## NCBI taxdump

If you would rather use NCBI taxonomy instead of the GTDB taxonomy:

```
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -pzxvf taxdump.tar.gz
```

## Getting reference genomes for the custom databases

### Downloading genomes

* If using [GTDB](https://gtdb.ecogenomic.org/) genomes, run `GTDB_metadata_filter.R` to select genomes
* If downloading genomes from genbank/refseq, you can use `genome_download.R`
* You can also include your own genomes (e.g., newly created MAGs)

An example:

```
# Filtering GTDB metadata to certain genomes
./GTDB_metadata_filter.R -o gtdb-r95_bac-arc.tsv https://data.gtdb.ecogenomic.org/releases/release95/95.0/ar122_metadata_r95.tar.gz https://data.gtdb.ecogenomic.org/releases/release95/95.0/bac120_metadata_r95.tar.gz

# Downloading all genomes (& creating tab-delim table of genome info)
./genome_download.R -o genomes -p 8 gtdb-r95_bac-arc.tsv > genomes.txt

# Note: the output of ./genome_download.R can be directly used for running the `Struo2` pipeline (see below)
# Note: genome fasta files can be compressed (gzip or bzip2) or uncompress for input to Struo2
```

## Input genome data (`samples.txt` file)

The table of input files/data can be created using the helper scripts described above
if you downloaded the genomes from the GTDB or NCBI.

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

## Running the pipeline

> This applies to creating and updating the databases

### Edit the config file

* Select a config to edit
  * If creating a new database: edit `config.yaml`
  * If updating an existing database: edit `config-update.yaml`
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
    * `all_genes_annot.fna.gz`
  * amino acid
    * `all_genes_annot.dmnd`
  
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
* metaphlan3 custom database using genes_db info
  * using gene cluster & gene taoxnomy to determine clade-level inclusion
  * method
    * build DAG of taxonomy (eg., taxdump + taxids + GTDB metadata)
      * attributes for tips (genomes)
        * ncbi_total_length
    * for each gene cluster:
      * place on taxonomy
        * cluster membership => genes => taxid
      * determine LCA
        * must be inclusive of N% of descendents
	* [optional] for each descendent internal node determine externals
    * format for metaphlan3
      * sequence headers
        * standard marker_name = `(NCBI_taxid)(UniRef90_cluster)(CDS_name)`
          * names are actually arbitrary but must match the keys in `db['markers'][new_marker_name]`	 
      * pkl
        * taxonomy
	  * `db['taxonomy'][FULL_TAXONOMY] = [genome_taxid, genome_length]`
	* markers
	  * clade: leaf of the taxonomy (genome name)
	  * ext: non-target genomes (full taxonomy)
	  * len: average length of the marker
	  * taxon: full taxonomy of target clade
* metaphlan3 custom database using humann_db info
  * using all genome-derep-annot genes:
  * for each species in the genome taxonomy (provided by user-samples)
    * determine which markers are shared by all members
      * sharing based on uniref IDs
    * for markers that are just 'nearly' unique to a species:
      * `ext` => all genomes "external" to the clade
      * rule of thumb: <=10 "external" genomes

