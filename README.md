# LinkSeq-Demux -- basecall/demultiplex and trim linked-reads

[![Docker build](https://img.shields.io/badge/Docker%20build-Available-informational)](https://hub.docker.com/repository/docker/fargen/linkseq-demux)

## Table of Contents
* [Overview](#overview)
* [Trimming](#trimming)
* [Setup](#setup)
* [Running on tiny-bcl](#running-on-tiny-bcl)
    * [Running demux pipeline](#run-demux-pipeline)
    * [Output](#output)
* [Barcode whitelist](#barcode-whitelist)

## Overview

This Nextflow pipeline basecalls and demultiplexes [linked-reads from 10x Genomics](https://www.10xgenomics.com/linked-reads/) and trims sequences, including barcode contamination. To run this pipeline, the 8-base sample indexes are needed, corresponding to the 10x Genomics indexes (e.g. `SI-GA-A1`).

This pipeline makes some assumptions about the input data. For example, it makes the assumption that it is paired-end sequencing, and therefore uses `--use-bases-mask=Y*,I*,Y*` in `bcl2fastq`, and assumes that the read lengths (and index length) is found in `RunInfo.xml`.

## Trimming

There are four trimming steps in the nextflow pipeline, each of which is listed below.

* `trim_adapters` trims adapter sequences from 3' end using [BBtools/BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).
* `bctrim` trims linked-read barcodes from read 2. When the insert size is small and read 1 and 2 overlap, read 2 may be contaminated by the barcode attached to read 1. The `bin/trimR2bc.py` was taken from release `v1.0` of [trimR2bc](https://github.com/ElisabetThomsen/trimR2bc/releases/tag/v1.0) by Elisabet Thomsen. See the GitHub repository for more details: https://github.com/ElisabetThomsen/trimR2bc
* `polyG_trim` trims poly-G tails from reads using [fastp](https://github.com/OpenGene/fastp).
* `quality_trim_read1`/`quality_trim_read2` trims low quality bases from read 1 and 2 respectively.

## Setup

Install dependencies with `conda` using the `environment.yml` file:

```
conda env create -f environment.yml
```

Activate the environment (check the name of the environment, it should be `linkseq-demux`):

```
conda activate linkseq-demux
```

Pull this project with `nextflow`:

```
nextflow pull https://github.com/fargenfo/linkseq-demux
```

## Running on tiny-bcl

Here's how to run this pipeline on the "tiny-bcl" example dataset from 10x Genomics. First of all, download the tiny-bcl tar file and the Illumina Experiment Manager sample sheet: tiny-bcl-samplesheet-2.1.0.csv.

> Download the tiny-bcl data:
> https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/mkfastq#example_data

Next, edit the `[Data]` part of the samplesheet from the following:

```
Lane,Sample_ID,index,Sample_Project
5,Sample1,SI-GA-C5,tiny_bcl
```

To the following:

```
Lane,Sample_ID,index
5,Sample1,CGACTTGA
5,Sample1,TACAGACT
5,Sample1,ATTGCGTG
5,Sample1,GCGTACAC
```

Using that the index `SI-GA-C5` corresponds to the four octamers `CGACTTGA,TACAGACT,ATTGCGTG,GCGTACAC` (10X Genomics have a tool for this on [their website](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/bcl2fastq-direct)).

Also add adapter sequences to the `[Settings]` part of the samplesheet:
```
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

### Run demux pipeline

Use the `tiny-bcl.config` file to get an idea how to define input parameters. For the `tiny-bcl` dataset, you should set `rundir` to `tiny-bcl-2.0.0` and `samplesheet` to `tiny-bcl-samplesheet-2.1.0.csv`. Additionally, you need the barcode whitelist, see the [barcode whitelist](https://github.com/fargenfo/linkseq#barcode-whitelist) section on how to obtain this list.

Change the memory and CPU specifications in the configuration (under `process` and `executor`) to suit your needs before continuing.

When you've made the config file and activated the `linkseq-demux` environment, run the pipeline like this:

```
nextflow run fargenfo/linkseq-demux -c [your config]
```

### Output

The output from the pipeline, run on the `tiny-bcl` data, is shown below. The compressed FASTQ data is in `outs/Sample1/fastqs`. There are various logs, from the basecalling itself via `Bcl2Fastq`, from the various trimming steps, from read synchronization, and from `FastQC`. There is also a HTML report `outs/multiqc/multiqc_report.html` that combines the statistics from `Bcl2Fastq` with the `FastQC` report.

```
$ tree outs/
outs/
├── Sample1
│   ├── fastqc
│   │   ├── Sample1_L005_R1_fastqc.html
│   │   ├── Sample1_L005_R2_fastqc.html
│   │   ├── fastqc.log
│   │   └── zips
│   │       ├── Sample1_L005_R1_fastqc.zip
│   │       └── Sample1_L005_R2_fastqc.zip
│   ├── fastqs
│   │   ├── Sample1_L005_R1.fastq.gz
│   │   └── Sample1_L005_R2.fastq.gz
│   └── logs
│       ├── adapter_trim
│       │   └── L005.log
│       ├── bctrim
│       │   └── L005.log
│       ├── polyG_trim
│       │   └── L005.log
│       ├── quality_trim
│       │   ├── L005_R1.log
│       │   └── L005_R2.log
│       └── sync_reads
│           └── L005.log
├── bcl2fastq
│   ├── Stats.json
│   └── log.log
└── multiqc
    ├── multiqc_data
    │   ├── multiqc.log
    │   ├── multiqc_bcl2fastq_bylane.txt
    │   ├── multiqc_bcl2fastq_bysample.txt
    │   ├── multiqc_data.json
    │   ├── multiqc_fastqc.txt
    │   ├── multiqc_general_stats.txt
    │   ├── multiqc_qualimap_bamqc_genome_results.txt
    │   └── multiqc_sources.txt
    └── multiqc_report.html
```

## Barcode whitelist

The whitelist contains barcodes in the 10X Genomics GemCode technology, and can be obtained in the LongRanger software bundle here:

longranger-2.y.z/longranger-cs/2.y.z/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt

LongRanger can be obtained on the 10X Genomics website:

https://support.10xgenomics.com/genome-exome/software/overview/welcome

It's easier to download the file from this link though:
```
http://cb.csail.mit.edu/cb/ema/data/4M-with-alts-february-2016.txt
```
