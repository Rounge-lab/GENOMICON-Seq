# GENOMICON-Seq: GENOmics Modeling, In-silico CONstruction and Sequencing

## Table of content

1. [Introduction](#Introduction)
2. [Quick set-up](#quick-set-up)
3. [Quick start](#quick-start)
4. [Computational power and processing time](#computational-power-and-processing-time)
4. [Simulation parameters](#simulation-parameters)
5. [Amplicon sequencing simulation parameters](#Amplicon-sequencing-simulation-parameters)
6. [WES sequencing simulation parameters](#WES-sequencing-simulation-parameters)
7. [Main inputs](#Main-inputs)
8. [Main outputs](#Main-outputs)
9. [Amplicon sequencing main outputs](#amplicon-sequencing-main-ouputs)
10. [WES sequencing main outputs](#WES-sequencing-main-ouputs)

## Introduction

Welcome to GENOMICON-Seq, a tool for genetic exploration through simulation. GENOMICON-Seq is designed to simulate both amplicon sequencing and whole exome sequencing (WES), providing a robust platform for users to experiment with virtual genetic samples.

At the core of GENOMICON-Seq is the ability to generate samples with varying mutation frequencies, which are then subjected to a simulated library preparation process. This unique feature allows users to introduce and adjust the level of noise within each simulation, effectively masking the detection of real mutations to various extents. Such simulations are critical in demonstrating how different noise levels can obscure mutation detection, offering valuable insights for refining laboratory techniques and enhancing mutation analysis.

Users can explore a range of different scenarios facilitating a deeper understanding of the dynamic interplay between genetic variations and the technical limitations of sequencing technologies.

Dive into GENOMICON-Seq and harness the power of simulation to advance your mutation detection strategies and guide improvements in laboratory practices. 


The GENOMICON-Seq allows users to:
A. control the generation of the initial sample with specified fasta sequences, their copy numbers, and mutation insertions. Users can define the mutation mode, mutation targets, and substitution types thereby fully controlling the mutational process; 
B. control the fraction of the initial sample that will be fragmented and the length of the produced fragments. As each fragmentation is a unique random process, users can specify the number of fragmentation replicates generated for each sample.
C. chose between amplicon or whole exome sequencing. 
- In amplicon sequencing users control the number of simulated PCR reactions for each fragmented replicate, the polymerase error rate in the PCR reactions, the number of cycles and the amplification efficiency.
- In WES users control the matching length between the probe and simulated fragments that influences the selection of fragments for sequencing. In addition, a simplified PCR mimics the short indexing PCR excluding the polymerase error rate.  
D. control the sequencing process, by specifying the number of reads that will be produced, and whether the fragment length will impact their probability of being sequenced (shorter fragment - higher probability of being sequenced). We modified the [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (1.6.0) tool enabling it to produce the reads from generated fragments but preserving its ability to simulate the sequencing error of different Illumina instruments.

The main outputs from each simulation are the generated FASTQ.GZ files, and overview of the mutations and their frequencies after the sample generation and after the simulated library preparation processes.  

Here we will show how to set up and run the simulations, briefly introduce the parameters, main inputs and main outputs. However, for a more detailed explanation about the design and theoretical background please see the full GENOMICON-Seq [manual](manual/GENOMICON_manual.pdf). 

## Quick set-up

The GENOMICON-Seq represents a collection of scripts written in R, C++, PowerShell Python, and well-known bioinformatic tools, while the whole process is orchestrated by Snakemake. Currently, it is only available as a Docker container. 

Prerequisites:
- Ensure that Docker is installed on your system before proceeding with the installation.
- Approximately 12 GB of free space (for testing purposes)
- A multi-core newer generation processor (at least 6 cores)
- At least 16GB RAM 

To quickly set up GENOMICON-Seq, execute the following commands in your terminal. This will download and run the script that sets up the Docker container, play dataset, Snakefiles, and their corresponding configuration files:

```
curl -L https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/genomicon_setup.sh -o genomicon_setup.sh
chmod +x genomicon_setup.sh
./genomicon_setup.sh
```

Use arguments `--ampliseq` or `--wes` when running the `genomicon_setup.sh` to download files/scripts necessary for either amplicon sequencing or WES simulation. Without the arguments, the script will download files/scripts for both simulations. 


A directory structure would look like this after running the `genomicon_setup.sh`

```
GENOMICON-Seq
├── config_ampliseq.yml
├── config_wes.yml
├── input_data_ampliseq
│   ├── all_primers.csv
│   ├── HPV16.fa
│   ├── HPV16.fa.fai
│   ├── multifasta_copy_number.csv
│   ├── primer_set_1.csv
│   ├── primer_set_2.csv
│   └── probability_variant_table.csv
├── input_data_wes
│   ├── chr1.fasta.gz
│   ├── dummy_probes.fasta
│   ├── fasta_lengths.csv
│   ├── multifasta_copy_number.csv
│   ├── probability_variant_table.csv
│   ├── SBS2_profile_short.csv
│   ├── xgen-exome-research-panel-v2-probes-hg38_filtered_chr1.bed
│   ├── xgen-exome-research-panel-v2-probes-hg38_short_all.bed
│   ├── xgen-exome-research-panel-v2-targets-hg38_filtered_chr1.bed
│   └── xgen-exome-research-panel-v2-targets-hg38_short_all.bed
├── Snakefile_ampliseq
├── Snakefile_wes
└── SQL_database
    ├── chr1.sqlite
    └── HPV16REF.sqlite
```

Play dataset is meant for testing. 
For amplicon sequencing simulation, human papillomavirus (HPV) type 16 will be used as an example. The necessary files will be in the `input_data_ampliseq`. The downloaded `config_ampliseq.yml` is set up for a test run. Check the chapter [Main inputs](#Main-inputs) for more information about each input file, and chapter [Amplicon sequencing simulation parameters](#Amplicon-sequencing-simulation-parameters) for detailed overview of the parameters and their usage.

For WES simulation, chromosome 1 sequence will be used as an example. The necessary files will be in the `input_data_wes`. The downloaded `config_wes.yml` is set up for a test run. Check the chapter [Main inputs](#Main-inputs) for more information about each input file, and [WES sequencing simulation parameters](#WES-sequencing-simulation-parameters) for detailed overview of the parameters and their usage.

The scripts necessary for running the simulations are integrated into the container, but can be accessed here: (scripts)[scripts/]. Read generation is governed by modified [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (v1.6.0), modified InSilicoSeq scripts can be found here: [modified_iss_scripts](modified_iss_scripts/).

In addition, during the simulation run, SQL-databases for each sequence in a fasta file will be generated. While for small genomes (several Mb) the generation time will not be long, for `WES simulation` on a `whole human genome`, generation time might take several hours. The premade SQL-database for each hg38 chromosome is available for downloading here [SQL_database_whole_human_genome](https://zenodo.org/records/12683302/files/SQL_database_human.tar.gz). The packed database size is 21.8 GB, when unpacked its size is 99.8 GB. Once made (or downloaded), the generation step of the SQL_database will be omitted as long as the headers of of the fasta sequences correspond to the names of the SQL-database.  

We also recommend downloading the human genome (hg38) fasta file where chromosome headers correspond to the SQL-database names here: [human_genome_hg38](https://zenodo.org/records/12683302/files/human_genome.fasta.gz)

## Quick start

As the simulation would require some substantial computational power (we will cover that in the next chapter [Computational power and processing time](#computational-power-and-processing-time)), a play dataset will be used for a test run. The genome for amplicon sequencing simulation is the double-stranded DNA, HPV16 genome, 7906bp. The parameters in the confing_ampliseq.yml are specified for testing purposes, each parameter will be explained in the next chapter. To run the simulation, run the following command from the main folder where the Snakemake_ampliseq and config_ampliseq.yml are:

```
docker run -it -v $(pwd):/usr/src/app/pipeline -w /usr/src/app/pipeline mimsto86/genomicon-seq:v1.0 snakemake -j 1 -p -s ./Snakefile_ampliseq --configfile ./config_ampliseq.yml
``` 

For simplicity and testing purposes, only human chromosome 1 will be used in the WES simulation. For the whole human genome, greater computational power will be needed. The command for the test run is similar to the amplicon sequencing simulation.

```
docker run -it -v $(pwd):/usr/src/app/pipeline -w /usr/src/app/pipeline mimsto86/genomicon-seq:v1.0 snakemake -j 1 -p -s ./Snakefile_wes --configfile ./config_wes.yml
```

A new folder `Output_data_ampliseq` or `Output_data_wes` will be created, depending on the chosen simulation. Main outputs for both simulations will be covered in the [Main outputs](#Main-outpus). 

## Computational power and processing time

GENOMICON-Seq is designed to run in High Performance Computing (HPC) environments due to its computational demands. For testing and lighter tasks, it can be also run on personal computers. Here are the details of the systems we have used for testing and development:

MacOS Sonoma 14.5

Processor: 2.6 GHz Six-core Intel Core i7 (12 logical cores)
RAM: 16 GB DDR4

Ubuntu 22.04.3 (jammy)

Processor: 2.6 GHz Six-core Intel Core i7 (12 logical cores)
RAM: 16 GB DDR4 

0n a HPC (we used [Educloud](https://www.uio.no/english/services/it/research/platforms/edu-research/)), for larger number of genome copies in amplicon sequencing simulation (we used 200 000 genome copies of HPV16) we used 60 CPUs and 4GB per CPU in our runs. In WES runs we used 60 CPUs and 8GB per run (when the copy number of complete human genome hg38 or only chromosome 1 was 2500). Depending on whether the complete genome was run or only chromosome 1, and the parameter values the simulation time was approximately 18-25h for the complete human genome, and approximately 1h for chromosome 1.  

GENOMICON-Seq includes a wide range of adjustable parameters that can significantly affect the tool's performance and resource demands. Due to this variability, it is challenging to predict the exact memory requirements, optimal number of processor cores, or the time needed to complete the processes. We recommend conducting initial tests for each new scenario to determine the appropriate configurations that best meet your computational and analysis needs.

## Simulation parameters

Most of the parameters are identical between simulations, however due to the nature of the simulated processes, some parameters are specific for amplicon sequencing simulation, while others can be found only in the WES simulation.

## Amplicon sequencing simulation parameters

For more details about parameters, see the user [manual](manual/GENOMICON_manual.pdf).
Parameters can be `OPTIONAL`, `MANDATORY`, `MUTUALLY EXCLUSIVE`, or `ADVANCED`
If parameter requires the double quotation mark enclosure, it will be noted. 
When the parameter should be omitted, set it up as as en empty string `""`.
The values in the parameters presented here are only examples.
Pre-configured config_ampliseq.yml can be used for testing purposes.

### I. General parameters

`NUM_SAMPLES: 1`

Specify the number the samples to be generated (`integer` ,do not enclose the value with double quotation marks) - `MANDATORY` parameter

`FRAG_REPLICATES: 1` 

Specify the number of fragmentation technical replicates to be generated from each sample (`integer`, do not enclose the value with double quotation marks) -`MANDATORY`

`NCORE: 8`  

Specify the number of CPUs that the simulation will use (`integer`, do not enclose the value with double quotation marks) - `MANDATORY`

`SEED: "23"` 

Specify the seed for reproducibility of the simulation - `OPTIONAL` 
If not specified, random seed will be assigned and recorded (`integer`, enclose the value with quotation marks)

`OUTPUT_DIRECTORY: "/usr/src/app/pipeline/Output_data_ampliseq" `

Specify the output directory - requires absolute path, edit only the name of the directory not the path `/usr/src/app/pipeline/` as it points the absolute path in the container. 
The parameter is `MANDATORY` (enclose the path with double quotation marks).
Change the name of the output for each new run! Or remove the previously generated directory before the new run to avoid rewriting, and possible mix-up of the outputs!

`INPUT_DIRECTORY: "input_data_ampliseq"` 

Specify the input directory - requires only the name of the directory
When using the `genomiconseq_setup.sh`, the directory name by default is `input_data_ampliseq` 
The parameter is `MANDATORY` (enclose the name with double quotation marks).

`FASTA_FILE: "HPV16.fa"` 

Specify the name of the fasta file - file needs to be placed in the input directory
The parameter is `MANDATORY` (enclose the name with double quotation marks).
Multifasta files are allowed.


### II. Sample parameters

#### 1. Genome topology

`CIRCULAR: "Y"` 

Specify the topology of the given genome (`"Y"` or `"N"`) - `OPTIONAL` (enclose the option with double quotation marks)
If not specified, the genome is considered to be circular (`"Y"`)

#### 2. Genome copies in the sample

`NR_COPIES: "10"` 

Specify the copy number of genomes to be generated - `OPTIONAL` - default is `1000` (`integer`, enclose the value with double quotation marks).
Specified number of copies applies on all sequences in the fasta file.

`MULTIFASTA_COPIES_TABLE: "multifasta_copy_numbers.csv"`

Specify the `csv` table containing the copy number of each sequence in the fasta file - `OPTIONAL` (enclose the table name with double quotation marks)
`NR_COPIES` and `MULTIFASTA_COPIES_TABLE` are `MUTUALLY EXCLUSIVE`
Read the full [manual](manual/GENOMICON_manual.pdf) for detailed description of the table format.

#### 3. Mutation targets

`MUT_CONTEXT: "ATC,CTA"` 

Specify mutational targets as nucleotide contexts (single nucleotide, tri or pentanucleotide context) - `OPTIONAL`
If not specified, any nucleotide in the sequence can mutate with the equal probability.
Accepted format: multiple context can be specified as a comma-separated string - All specified contexts must be of the same type (single, tri or penta)

#### 4. Substitution type

`TS_TV_RATIO: "3:1"`

Specify the transition vs transversion ratio (`colon-separated integers`, enclosed with double double quotation marks) - `OPTIONAL`
If none are specified (empty string ""), every nucleotide can be replaced with another with equal probability.

`SUBSTITUTION_PROBABILITY_TABLE: "probability_variant_table.csv"` 

Specify the name of the table with pre-determined substitution probabilities - `OPTIONAL` (enclose the table name with double quotation marks)
Example of the table can be found in the folder `input_data_ampliseq` - `probability_variant_table.csv`
`TS_TV_RATIO` and `SUBSTITUTION_PROBABILITY_TABLE` parameters are `MUTUALLY EXCLUSIVE`
If none are specified (empty string `""`), every nucleotide can be replaced with another with equal probability.

#### 5. Mutation mode

Two available mutation modes: `Deterministic` & `Specific mutation rate`.
Modes are `MUTUALLY EXCLUSIVE` - leave the empty strings (`""`) for the unused mode's parameters. 
When no mutations should be inserted - leave the `deterministic mode`'s parameters empty, and set `MUTATION_RATE` parameter to `"0"`
For detailed description on how these modes operates, see the full [manual](manual/GENOMICON_manual.pdf).

#### 5.a Deterministic Mode

`MUT_GENOME_FRACTION: "0.2"`

Specify the fraction of all genome copies that will be mutated (`0 <= float <= 1`, enclose the value with double quotation marks)

`FRACTION_POSITION: "0.1"` 

Specify the fraction of all targeted positions that will mutate (`0 <= float <= 1`, enclose the value with quotation marks).
If `MUT_CONTEXT` is specified, only the targeted contexts will be taken into account.

#### 5.b Specific mutation rate mode

`MUTATION_RATE: "1e-6"` 

Specify the mutation rate (`scientific number`, e.g `"1e-6"`, enclose the value with double quotation marks).

### III. Fragmentation parameters

`FRAGMENT_FRACTION: "1"`

Specify the fraction of all genome copies in a sample that will be fragmented (`0 <= float <= 1`, enclose the value with double quotation marks) - `OPTIONAL`
If not specified, the default is `1`.

`FRAGMENT_LENGTH: "250:550"` or `350`

Specify the length of the fragments to be produced - `OPTIONAL` (enclose the value with double quotation marks).
It can be a single non-negative `integer` value where all produced fragments will have the same length, e.g. `"450"`.
Or it can be the range of values written as a colon-separated `integers` values (e.g `"250:1000"`)
If not specified (empty string "") - by default range of fragment length is `"250:450"`.

### IV. PCR parameters

`PRIMERS: "all_primers.csv"` 

Specify the `csv` table containing the primer information - `MANDATORY` (enclose the table name with double quotation marks).
See the full [manual](manual/GENOMICON_manual.pdf) for the detailed description of the file, chapter Input files.
Or check the `all_primers.csv` file in the input directory that contains the HPV16-primers that can be used in test runs.
`NOTE:` In addition to this file PCR requires additional csv files not specified in parameters, one for each specified PCR reaction (see the [Main inputs](#Main-inputs) or full [manual](manual/GENOMICON_manual.pdf) chapter Input files).

`NUM_PCR: 2` 

Specify the number of PCR reactions to be produced (`integer`, do not enclose the value with double quotation marks) - `OPTIONAL`
By default the number of PCR reactions is 2

`AMPLICON_LENGTH: 250` 

Specify the minimum amplicon length to be produced - `MANDATORY` (`integer`, do not enclose the value with double quotation marks)

`NUM_CYCLES: "25"`

Specify the number of PCR cycles to be performed - `OPTIONAL` (`integer`, enclose the value with double quotation marks).
The default is 25.

`POL_ERROR_RATE: "1e-5"`

Specify the polymerase error rate - `OPTIONAL` (`scientific number`, e.g. `"1e-5"`, enclose the value with double quotation marks).
If not specified, the default is 1e-6 

`FRACTION_PCR_SEQUENCING: "0.005"` 

Specify the fraction of PCR reaction that will be sequenced - `MANDATORY` (`0 <= float <= 1`, enclose the value with double quotation marks).

#### Advanced PCR parameters 
See to [manual](manual/GENOMICON_manual.pdf) for detailed description of these parameters!

`K_PARAMETER_pcr: "1"` 

The k-parameter shapes the steepness of the simulated PCR efficiency drop.
`0 < float`, enclose the value with the double quotation marks.
`K_PARAMETER` and `MIDPOINT_CYCLE` (next parameter) are responsible for the amplification efficiency drop during the PCR cycling - `OPTIONAL`.
The default is 1.

`MIDPOINT_CYCLE: "20"`

Specify the PCR cycle at which the amplification efficiency is 50% - `OPTIONAL` 
The default is a cycle corresponding to the 60% of all specified cycles.

### V. Sequencing parameters

`NUM_READS: "100000"` 

Specify the number of reads to be produced (`integer`, enclose the value with double quotation marks) - `MANDATORY`

`MODE: "kde"`

Specify the sequencing mode - `OPTIONAL`
[InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (v1.6.0) supports `2 modes`, `"basic"` / `"perfect"` or `"kde"` (enclose the given option with double quotation marks).
`basic/perfect` mode initiate the production of error-free reeds.
`kde` is the default option, requires the specification of the ERROR_MODEL.

`ERROR_MODEL: "novaseq"`

Specify the sequencing error model - `OPTIONAL`
[InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (v1.6.0) supports 4 build-in error models `HiSeq`, `NextSeq`, `NovaSeq`, and `MiSeq` (enclose the given option with quotation marks)
The default in `None`.

`FASTA_GZ_OUTPUT: "--compress"` 

Specify whether the produced fastq files should be gz-compressed - `OPTIONAL`
Leave the empty string (such as "") to omit the compression

`GC_BIAS: "--gc_bias"`

Specify whether the GC-bias should be applied
If applied, reads with abnormal GC-counts might not be produced and the number of generated reads might be lower than specified.
The parameter is `OPTIONAL`, leave the empty string "" to omit the bias.

For more details about `MODE`, `ERROR_MODEL`, `FASTA_GZ_OUTPUT`, and `GC_BIAS` parameters please refer to the [InSilicoSeq documentation](https://insilicoseq.readthedocs.io/en/latest/index.html).

`OPT_FRAG_LENGTH: "NO"` 

Specify the optimal fragment length mode - `MANDATORY` 
`OPT_FRAG_LENGTH` alters the probability of a fragment to be sequenced based on its length
For each specified mode, value represents the length at which fragment has 50% probability to be sequenced.
Fragments with greater length than the one specified will have less chance for sequencing and vice versa.
Each chosen mode must be enclosed by double quotation marks
5 available modes:
   1. `"NO"` - mode will not be applied
   2. `"fixed {integer}"` (e.g. `"fixed 450"`) - sets the length to 450 bp, fragments with this length will have 50% to be sequenced
   3. `"default"` - the length of 350 bp will have 50% chance to sequenced
   4. `"medain"`  - median length based on the length of all fragments selected to be sequenced is determined. Fragments with this length will have 50% chance to be sequenced.
   5. `"quartil {1|2|3|4}"` (e.g. `"quartil 3"`) - value of the selected quartile will be determined based on the length of all fragments. Fragments with this length will have 50% chance to be sequenced.

#### Advanced Sequencing parameter 
See to [manual](manual/GENOMICON_manual.pdf) for detailed description of the parameter.

`K_PARAMETER_seq: "0.005"`

The parameter is `MANDATORY` only when `OPT_FRAG_LENGTH` is other than `“NO”`.
We recommend the value 0.005
`0 < float`, enclose the value with the double quotation marks.

## WES sequencing simulation parameters

For more details about parameters, see the user [manual](manual/GENOMICON_manual.pdf).
Parameters can be `OPTIONAL`, `MANDATORY`, `MUTUALLY EXCLUSIVE`, or `ADVANCED`
If parameter requires the double quotation mark enclosure, it will be noted. 
When the parameter should be omitted, set it up as as en empty string `""`.
The values in the parameters presented here are only examples.
Pre-configured config_wes.yml can be used for testing purposes.

### I. General parameters

`NUM_SAMPLES: 1`

Specify the number the samples to be generated (`integer` ,do not enclose the value with double quotation marks) - `MANDATORY` 

`FRAG_REPLICATES: 1` 

Specify the number of fragmentation technical replicates to be generated from each sample (`integer`, do not enclose the value with double quotation marks) -`MANDATORY`

`NCORE: 8`  

Specify the number of CPUs that the simulation will use (`integer`, do not enclose the value with double quotation marks) - `MANDATORY` 

`SEED: "23"` 

Specify the seed for reproducibility of the simulation - `OPTIONAL` 
If not specified random seed will be assigned and recorded (`integer`, enclose the value with quotation marks)

`OUTPUT_DIRECTORY: "/usr/src/app/pipeline/Output_data_ampliseq" `

Specify the output directory - requires absolute path, edit only the name of the directory not the `path /usr/src/app/pipeline/`
The parameter is `MANDATORY` (enclose the path with double quotation marks).
Change the name of the output for each new run! Or remove the previously generated directory before the new run to avoid rewriting, and possible mix-up of the outputs!

`INPUT_DIRECTORY: "input_data_ampliseq"` 

Specify the input directory - requires only the name of the directory
When using the `genomiconseq_setup.sh`, the directory name by default is `input_data_ampliseq` 
The parameter is `MANDATORY` (enclose the name with double quotation marks).

`FASTA_FILE: "HPV16.fa"` 

Specify the name of the fasta file - file needs to be placed in the input directory
The parameter is `MANDATORY` (enclose the name with double quotation marks).
Multifasta files are allowed.

### II. Sample parameters

#### 1. Genome copies in the sample

`NR_COPIES: "50"` 

Specify the copy number of genomes to be generated - `OPTIONAL` - default is `100` (`integer`, enclose the value with double quotation marks).
Specified number of copies applies on all sequences in the fasta file.

`MULTIFASTA_COPIES_TABLE: "multifasta_copy_numbers.csv"`

Specify the `csv` table containing the copy number of each sequence in the fasta file - `OPTIONAL` (enclose the table name with double quotation marks)
`NR_COPIES` and `MULTIFASTA_COPIES_TABLE` are `MUTUALLY EXCLUSIVE`
Read the full [manual](manual/GENOMICON_manual.pdf) for detailed description of the table format.

#### 2. Mutation targets

`MUT_CONTEXT: "ATC,CTA"` 

Specify mutational targets as nucleotide contexts (single, tri or pentanucleotide context) - `OPTIONAL`
If not specified, any nucleotide in the sequence can mutate with the equal probability.
Accepted format: multiple context can be specified as a comma-separated string - All specified contexts must be of the same type (single, tri or penta)

#### 3. Substitution type

`TS_TV_RATIO: "3:1"`

Specify the transition vs transversion ratio (`colon-separated integers`, enclosed with double double quotation marks) - `OPTIONAL`
If none are specified (empty string ""), every nucleotide can be replaced with another with equal probability.

`SUBSTITUTION_PROBABILITY_TABLE: "probability_variant_table.csv"` 

Specify the name of the table with pre-determined substitution probabilities - `OPTIONAL` (enclose the table name with double quotation marks)
Example of the table can be found in the folder `input_data_ampliseq` - `probability_variant_table.csv`
`TS_TV_RATIO` and `SUBSTITUTION_PROBABILITY_TABLE` parameters are `MUTUALLY EXCLUSIVE`
If none are specified (empty string ""), every nucleotide can be replaced with another with equal probability.

#### 4. Mutation mode

Three mutation modes are available: `Deterministic`, `Specific mutation rate`, and `SBS-mimicry` 
`Deterministic` and `Specific mutation rate` modes are `MUTUALLY EXCLUSIVE` - leave the empty strings ("") for unused mode's parameters. 
`SBS-mimicry` mode refers to Single Nucleotide Substitution, and is the `extension of the Deterministic mode`, in addition to specifying the `SBS_TABLE`, `Deterministic mode`'s parameters needs to be specified!
When `SBS-mimicry` mode is used, also set `MUT_CONTEXT`, `TS_TV_RATIO`, `SUBSTITUTION_PROBABILITY_TABLE` to empty strings "".
When no mutations should be inserted - leave the `Deterministic` and `SBS-mimicry` mode parameters empty, and set `MUTATION_RATE` parameter to `"0"`
For detailed description on how these modes operates, see the [manual](manual/GENOMICON_manual.pdf).

#### 4.a Deterministic Mode

`MUT_GENOME_FRACTION: "0.2"`

Specify the fraction of all genome copies that will be mutated (`0 <= float <= 1`, enclose the value with double quotation marks)

`FRACTION_POSITION: "0.1"` 
Specify the fraction of all targeted positions that will mutate (`0 <= float <= 1`, enclose the value with quotation marks).
If `MUT_CONTEXT` is specifed, only the targeted contexts will be taken into account.

#### 4.b Specific mutation rate mode

`MUTATION_RATE: "1e-6"` 

Specify the mutation rate (`scientific number`, e.g `"1e-6"`, enclose the value with double quotation marks).

#### 4.c SBS-mimicry

`SBS_SIGNATURES: "sbs2_short.csv"`

Specify the name of the SBS-signature table. See the user [manual](manual/GENOMICON_manual.pdf) for the detailed description of the table format, and please use the [sbs_table_trimming.py](additional_scripts/sbs_table_trimming.py) that converts the original SBS-signature table to the simulation-required format. 

Download the SBS-signature table from [COSMIC](https://cancer.sanger.ac.uk/signatures/sbs/), save it as a csv file, for example SBS2_table.csv. Each table has values for each 96 mutational contexts stored in 5 columns `SBS2_GRCh37`, `SBS2_GRCh38`, `SBS2_mm9`, `SBS2_mm10`, `SBS2_rn6`, pick the column name (for hg37 or hg38), we used SBS2_GRCh38, and run the script

```
python sbs_table_trimming.py --input_sbs path/the/table/SBS2.csv --column_name SBS2_GRCh38 --output_sbs path/to/output/SBS2_selected.csv
``` 
Place the generated table into the `input_data_wes` directory and set up a table name in the `config_wes.yml`.

### III. Fragmentation parameters

`FRAGMENT_FRACTION: "1"`

Specify the fraction of all genome copies in a sample that will be fragmented (`0 <= float <= 1`, enclose the value with double quotation marks) - `OPTIONAL`
If not specified, the default is `1`.

`MIN_FRAGMENT_LENGTH: "250"` 

Specific the minimum fragment length - `MANDATORY`

`MAX_FRAGMENT_LENGTH: "1000"` 

Specific the maximum fragment length - `MANDATORY`


### IV. Probe capture enrichment parameters

`PROBES_FASTA: ""`

Specify the fasta file containing probe sequences. The file must be placed in the specified input folder (input_data_wes). The input folder does contain the dummy-probe sequences for testing purposes, that binds only to chromosome 1.

`PROBES_BED: "xgen-exome-research-panel-v2-probes-hg38_short_all.bed"`

Specify the BED file containing the binding coordinates of probes to the selected chromosomes.
`PROBES_FASTA` and `PROBES_BED` are `MUTUALLY EXCLUSIVE`, one or the other has to be specified. For testing purposes we are using the probe BED file from [xGen Exome Research Panel v2 hg38](https://eu.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-hybridization-capture/pre-designed-hyb-cap-panels/exome-hyb-panel-v2) placed in the input_data_wes folder when the `genomiconseq_intall.sh` is run.

PROBES_TARGES: "xgen-exome-research-panel-v2-targets-hg38_short_all.bed" 
Specify BED file containing coordinates of either probe targeted exons or all exons.
For testing purposes we are using BED file of [xGen Exome Research Panel v2 probe-targeting exons](https://eu.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-hybridization-capture/pre-designed-hyb-cap-panels/exome-hyb-panel-v2). BED file is placed in the `input_data_wes` folder when the `genomiconseq_intall.sh` is run

`MATCHING_LENGTH: "100"`

Specify the minimum matching length between the probe and the fragment - `MANDATORY`
Matching length selects the fragments for sequencing. 

### V. PCR Parameters 

`NUM_CYCLES: "25"`

Specify the number of PCR cycles to be performed - `OPTIONAL` (`integer`, enclose the value with double quotation marks).
The default is 25.

#### Advanced PCR parameters 
See to [manual](manual/GENOMICON_manual.pdf) for detailed description of these parameters!

`K_PARAMETER_pcr: "1"` 

The k-parameter shapes the steepness of the simulated PCR efficiency drop.
`0 < float`, enclose the value with the double quotation marks.
`K_PARAMETER` and `MIDPOINT_CYCLE` (next parameter) are responsible for the amplification efficiency drop during the PCR cycling - `OPTIONAL`.
The default is 1.

`MIDPOINT_CYCLE: "20"`

Specify the PCR cycle at which the amplification efficiency is 50% - `OPTIONAL` 
The default is a cycle corresponding to the 60% of all specified cycles.

### V. Sequencing parameters

`NUM_READS: "100000"` 

Specify the number of reads to be produced (`integer`, enclose the value with double quotation marks) - `MANDATORY`

`MODE: "kde"`

Specify the sequencing mode - `OPTIONAL`
[InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (v1.6.0) supports `2 modes`, `"basic"` / `"perfect"` or `"kde"` (enclose the given option with double quotation marks).
`basic/perfect` mode initiate the production of error-free reeds.
`kde` is the default option, requires the specification of the ERROR_MODEL.

`ERROR_MODEL: "novaseq"`

Specify the sequencing error model - `OPTIONAL`
[InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/index.html) (v1.6.0) supports 4 build-in error models `HiSeq`, `NextSeq`, `NovaSeq`, and `MiSeq` (enclose the given option with quotation marks)
The default in `None`.

`FASTA_GZ_OUTPUT: "--compress"` 

Specify whether the produced fastq files should be gz-compressed - `OPTIONAL`
Leave the empty string (such as "") to omit the compression

`GC_BIAS: "--gc_bias"`

Specify whether the GC-bias should be applied
If applied, reads with abnormal GC-counts might not be produced and the number of generated reads might be lower than specified.
The parameter is `OPTIONAL`, leave the empty string "" to omit the bias.

For more details about `MODE`, `ERROR_MODEL`, `FASTA_GZ_OUTPUT`, and `GC_BIAS` parameters please refer to the [InSilicoSeq documentation](https://insilicoseq.readthedocs.io/en/latest/index.html).

`OPT_FRAG_LENGTH: "NO"` 

Specify the optimal fragment length mode - `MANDATORY` 
`OPT_FRAG_LENGTH` alters the probability of a fragment to be sequenced based on its length
For each specified mode, value represents the length at which fragment has 50% probability to be sequenced.
Fragments with greater length than the one specified will have less chance for sequencing and vice versa.
Each chosen mode must be enclosed by double quotation marks
5 available modes:
   1. `"NO"` - mode will not be applied
   2. `"fixed {integer}"` (e.g. `"fixed 450"`) - sets the length to 450, fragments with this length will have 50% to be sequenced
   3. `"default"` - the length of 350 will have 50% chance to sequenced
   4. `"medain"`  - median length based on the length of all fragments selected to be sequenced is determined. Fragments with this length will have 50% chance to be sequenced.
   5. `"quartil {1|2|3|4}"` (e.g. `"quartil 3"`) - value of the selected quartile will be determined based on the length of all fragments. Fragments with this length will have 50% chance to be sequenced.

#### Advanced Sequencing parameter 
See to [manual](manual/GENOMICON_manual.pdf) for detailed description of the parameter.

`K_PARAMETER_seq: "0.005"`

The parameter is `MANDATORY` only when `OPT_FRAG_LENGTH` is other than `“NO”`.
We recommend the value 0.005
`0 < float`, enclose the value with the double quotation marks.


## Main inputs

For the more detailed info about each put file and their formats please see the user [manual](manual/GENOMICON_manual.pdf).

### a. Necessary input

`FASTA` or `FASTA.gz` file with the genome sequence(s) - multiple genomes/chromosomes can be stored in one FASTA file. use as few as possible characters for the genome name in the header of the FASAT file.

`CSV file with primer names and their sequences` required for `amplicon sequencing`, play dataset uses `all_primers.csv` file. This is the whole set of primers used for the whole genome sequencing of HPV16 in the [TaME-Seq method](https://github.com/jean-marc-costanzi/TaME-seq). The format of the table must be the same. The `“name”` column contains primer names, the `“F_R”` column contains information about the primer orientation, and `“seq”` is the primer sequence.

`CSV file with primer names`, one file per PCR reaction required for `amplicon sequencing` - Depending on the specified number of PCR reactions, for each reaction, a separate primer set has to be pre-made, and placed in `input_data_ampliseq`. In our `input_data_ampliseq` from the play dataset, we have 2 primer_set CSV files, `primer_set_1.csv` and `primer_set_2.csv`, the index number at the end specifies in which PCR reaction, primers should be used. This is a single-column table (column “name”) with ONLY primer names. In case of only one PCR reaction, all names stored in the all_primers.csv should also be in the primer_set_1.csv. The number of `primer_set CSV tables` must match the number of specified PCR reactions in `NUM_PCR`.

`FASTA` or `BED probe file` - required for `WES` - Either one of the files needs to be present in the `input_data_wes`. A BED file must have three columns, chromosome name, and start and end coordinates of the probe binding sites. 

`BED file with exon regions targeted by probes` - required for `WES` - The file is required for mutation generation. Generated mutations will be found only in the regions specified by this file.

`SQL database` - Besides the input files in the main input directory both `amplicon sequencing` and `WES`, the simulation requires an SQL database. In the chapter [Quick set-up](#quick-set-up), we specified the link from which the complete hg38 SQL_database can be downloaded, saving time time to generate it, as well as the h38 fasta file with all the chromosomes which names match the SQL database names. 

`SBS table` containing the percentage of single base substitutions (SBS) required only in `WES` when the `SBS-mimicry` mutation mode is used - For each specific SBS signature, a table in numeric form can be downloaded from [COSMIC](https://cancer.sanger.ac.uk/signatures/sbs/) and trimmed to contain only the required information. Non-trimmed tables contain the SBS percentage for hg37, hg38, mm9, mm10, and rn6. In chapter [WES sequencing simulation parameters](#WES-sequencing-simulation-parameters), we demonstrated how to convert the original table format to the format accepted by the simulation.

### b. Optional input

`CSV table with the copy number of different genomes` (`multifasta_copy_number.csv` in the `input_data` folder) - In the case where FASTA file contains many different genomes/chromosomes, their individual copy numbers can be specified in this table- otherwise, parameter `NR_COPIES` is applied on all genomes/chromosomes in the FASTA file. The table has two columns, `“fasta_name”` containing headers of each fasta sequence in the FASTA file, and `“copy_number”` containing the number of copies for each of the sequences.

`CSV table with the probability of each nucleotide being substituted with other nucleotides` - this table replaces the parameter `TS_TV_RATIO` (transitions vs transversions ratio). It enables the specification of the mutation probability of each nucleotide to mutate to any other nucleotide but itself. The table has 3 columns, `“Nucleotide”`, `“Variant”` and `“Probability”`, where for each nucleotide the probability to be mutated to any of the three variants is specified. 

## Main outputs

The complete list of the output files made in both amplicon sequencing and WES are presented in the tool [manual](manual/GENOMICON_manual.pdf). Here we outline the most important ones that will enable the mutation tracking.

### Amplicon sequencing main outputs

For amplicon sequencing, the output folder should look like this

```
Output_data_ampliseq
├── cleanup_2_complete.txt
├── cleanup_3_complete.txt
├── cleanup_4_complete.txt
├── cleanup_5_complete.txt
├── sample_1
│   ├── HPV16REF_inserted_mutations_overview.csv
│   ├── HPV16REF_mutation_counts.csv
│   ├── cleanup_1_complete.txt
│   ├── fasta_lengths.csv
│   ├── generated_reads
│   │   ├── tech_replicate_1_PCR_1_R1.fastq.gz
│   │   ├── tech_replicate_1_PCR_1_R2.fastq.gz
│   │   ├── tech_replicate_1_PCR_2_R1.fastq.gz
│   │   └── tech_replicate_1_PCR_2_R2.fastq.gz
│   ├── log_folder
│   │   ├── Defining_seq_and_mut.log
│   │   └── Making_SQL_database.log
│   └── tech_replicate_1
│       ├── HPV16REF_fragments.bed.gz
│       ├── PCR
│       │   ├── HPV16REF_splitting_complete.log
│       │   ├── PCR_1
│       │   │   ├── HPV16REF_seq_mutations_overview.csv
│       │   │   ├── PCR_1-polymerase_error_mutations.csv
│       │   │   └── PCR_1_sequenced_frags.fasta.gz
│       │   ├── PCR_1.log
│       │   ├── PCR_2
│       │   │   ├── HPV16REF_seq_mutations_overview.csv
│       │   │   ├── PCR_2-polymerase_error_mutations.csv
│       │   │   └── PCR_2_sequenced_frags.fasta.gz
│       │   ├── PCR_2.log
│       │   ├── fragments_pcr1.bed.gz
│       │   ├── fragments_pcr2.bed.gz
│       │   └── merging_complete.log
│       └── log_folder
│           └── Fragmentation.log
└── seed_log.txt
```

Here, the amplicon sequencing simulation encompassed 2 PCR reactions to achieve the whole genome sequencing of the HPV16. For each generated sample (specified by `NUM_SAMPLES`), generated FASTQ.GZ (or none-compressed FASTQ) can be found in a folder `generated_reads`. 
For each sequence in `FASTA_FILE` a `*_inserted_mutations_overview.csv` will be generated in `sample_*` folder containing all mutations that were inserted during the sample generation. The file will include the info about position, nucleotide, complementary nucleotide, trinucleotide, reverse complementary trinucleotide, pentanucleotide, reverse complementary pentanucleotide contexts, and the frequency of each mutation (number of initial genomes the mutation has been introduced to). File will not be produced if no mutations were introduced during the sample generation process. 
In addition, again for each sequence in `FASTA_FILE` a `*_mutation_counts.csv` provides the information about how many genome copies have been mutated and how many mutations each mutated copy contained. 


In `PCR_*` folder produced for each PCR reaction specified, `*_seq_mutations_overview.csv` will be generated for each sequence in the `FASTA_FILE`. The csv has the same format as the `*_inserted_mutations_overview.csv` file but the frequency of each mutation now represent the number of unique fragments that are sequenced. File enable the comparison of mutation frequencies before and after the simulated library preparation process and provides the insight into the mutations that can be expected to be sequenced. If no mutations were introduced during the sample generation process, or if all mutations has been lost during the library preparation process, only the simple log file will be produced.

If polymerase error rate has been specified (`POL_ERROR_RATE`), the information about the error mutations will be stored in `PCR_*-polymerase_error_mutations.csv` for each PCR reaction generated. The file holds the info about the mutation position, and to which nucleotide the current nucleotide at this position has been mutated to. 

### WES sequencing main outputs

For WES, the output folder should look like this

```
Output_data_wes
├── cleanup_1_complete.txt
├── cleanup_3_complete.txt
├── sample_1
│   ├── chr1_inserted_mutations_overview.csv
│   ├── chr1_mutation_counts.csv
│   ├── cleanup_2_complete.txt
│   ├── generated_reads
│   │   ├── tech_replicate_1_PCR_R1.fastq.gz
│   │   └── tech_replicate_1_PCR_R2.fastq.gz
│   ├── log_folder
│   │   ├── 0_Making_SQL_database.log
│   │   └── 1_Defining_seq_and_mut.log
│   └── tech_replicate_1
│       ├── PCR_filtered.bed.gz
│       ├── PCR_reaction
│       │   ├── chr1_seq_mutations_overview.csv
│       │   └── sequenced_frags.fasta.gz
│       └── log_folder
│           ├── 2_Fragmentation_prep.log
│           ├── PCR_reaction
│           │   └── PCR_reaction.log
│           ├── chr1_all_converted_to_csv.log
│           ├── chr1_all_filtered.log
│           ├── chr1_conversion_to_bedgz_complete.log
│           └── chr1_fragmentation_complete.log
└── seed_log.txt
```

Here, the WES simulation included chromosome 1. For each generated sample (specified by `NUM_SAMPLES`), generated FASTQ.GZ (or none-compressed FASTQ) can be found in a folder `generated_reads`. 
For each chromosome sequence in `FASTA_FILE` a `*_inserted_mutations_overview.csv` will be generated in `sample_*` folder, containing all mutations that were inserted during the sample generation. The file will include the info about position, nucleotide, complementary nucleotide, trinucleotide, reverse complementary trinucleotide, pentanucleotide, reverse complementary pentanucleotide contexts, and the frequency of each mutations (number of initial genomes the mutation has been introduced to). File will not be produced if no mutations were introduced during the sample generation process. 
In addition, again for each sequence in `FASTA_FILE` a `*_mutation_counts.csv` provides the information about how many genome copes has been mutated and how many mutations each mutated copy contained.

In `PCR_reaction` folder `*_seq_mutations_overview.csv` will be generated for each chromosome sequence in the `FASTA_FILE`. The csv has the same format as the `*_inserted_mutations_overview.csv` file but the frequency of each mutation now represent the number of unique fragments holding the mutation which were sequenced. File enable the comparison of mutation frequencies before and after the simulated library preparation process and provides the insight into the mutations that can be expected to be sequenced. If no mutations were introduced during the sample generation process, or if all mutations has been lost during the library preparation process, only the simple log file will be produced.
