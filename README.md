---
title: "README"
author: "Ram Krishna Shrestha"
date: "18/01/2017"
output: html_document
---


## Introduction

This is a small python project that brings some functions together while analysing next generation sequencing reads. In my experience, many times while analysing NGS data, I have to pick up a tool for one function. But I wanted one tool with all functions and can be controlled by options.

The script process_ngs.py has functions for general statistics of the reads, excess reads by interval/sequence ids, get sequence and length (two column data), get subsequence, clip sequence reads from 5' or 3' or both, reverse the sequence reads or reverse complement the sequence reads.

Users can combine any functions to get the output results. Please see the usage/examples below in Usage section.

## Requirement

1. python v2.6.7+
2. Biopython

## Usage

To get help about the options

> process_ngs.py -h

```
usage: process_fasta_fastq.py [-h] [-i INPUT] [--stats] [--interval INTERVAL]
                              [--seqid SEQID] [-l] [--filterbylength]
                              [--subseq] [-x MIN] [-y MAX]
                              [--leftclip LEFTCLIP] [--rightclip RIGHTCLIP]
                              [-r] [--reversecomp] [-o OUTPUT]

Program to process fasta or fastq file in wide range of aspects Check the help
for different types of available options.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Fasta or fastq input file
  --stats               Prints the basic statistics of the reads
  --interval INTERVAL   gets the sequences in interval specified. Provide
                        start point and then interval. e.g. 3,2
  --seqid SEQID         comma separated list of sequence id from input file
  -l, --getlength       Outputs the sequence ID and sequence length (in tab-
                        delimited)
  --filterbylength      filter reads by length
  --subseq              get subsequence from sequence reads. Default: gets
                        first 100 bps in every sequence
  -x MIN, --min MIN     provide minimum position for subsequence. Must provide
                        --subseq
  -y MAX, --max MAX     provide maximum position for subsequence. Must provide
                        --subseq
  --leftclip LEFTCLIP   removes [int] bases from left end or 5' end
  --rightclip RIGHTCLIP
                        removes [int] bases from right end or 3' end
  -r, --reverse         reverses the sequence, not reverse complement
  --reversecomp, --reverse_complement
                        reverse complements sequence reads
  --translate           Translate nucleotide seqeunce to amino acid sequence
  
  -o OUTPUT, --output OUTPUT
                        output filename
```

Some examples of usage:
```
To get the general statistics of the data

> process_ngs.py --input your_file --stats

To extract the sequence reads and get general stats of the extract reads only

> process_ngs.py --input your_file --seqid seqid1,seqid2 --stats

> process_ngs.py --input your_file --seqid file_with_list_of_seqids --stats

To extract reads with seqid and output reads as subreads from 10th base to 100th bp

> process_ngs.py --input your_file --seqid seqid1,seqid2 --subseq --min 10 --max 100

To extract reads with seqid and output reads as subreads from 10th base to 100th bp, reverse complement and then translate

> process_ngs.py --input your_file --seqid seqid1,seqid2 --subseq --min 10 --max 100 --reverse_complement --translate

To extract reads with seqid and output reads as subreads from 10th base

> process_ngs.py --input your_file --seqid seqid1,seqid2 --subseq --min 10

To extract reads with seqid and output reads as subreads from 10th base and get stats of the final output reads

> process_ngs.py --input your_file --seqid seqid1,seqid2 --subseq --min 10 --stats

To filter the reads by length

> process_ngs.py --input your_file --filterbylength --min 10 --max 100

Clip reads from 5' end

> process_ngs.py --input your_file --leftclip 5

Clip reads from 5' and 3' end

> process_ngs.py --input your_file --leftclip 5 --rightclip 3

Clip reads from 5' and 3' end and get statistics

> process_ngs.py --input your_file --leftclip 5 --rightclip 3 --stats

Clip reads from 5' and 3' and reverse the reads

> process_ngs.py --input your_file --leftclip 4 --rightclip 5 --reverse

Clip reads from 5' and 3' and reverse complement the reads

> process_ngs.py --input your_file --leftclip 4 --rightclip 5 --reversecomp

Extract reads occurring at every 5th position starting from 2nd read. If starting point is not provided, it starts from 1st read.

> process_ngs.py --input you_file --interval 2,5
```

As shown above, combination of some functions have been used. All functions that make sense can be combined.


## Further development

I would like to further add some more functions to the script (like quality trimming) in the future. Let me know if anyone wants to add functions to this script.