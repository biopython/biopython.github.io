---
title: Reading from unix pipes
permalink: wiki/Reading_from_unix_pipes
layout: wiki
tags:
 - Cookbook
---

Problem
-------

There are many circumstances when reading data from a Unix pipe is
preferable to reading data from a file. One example is reading sequences
from a compressed file, which is often preferable to uncompressing the
file and then reading from it.

You can also use pipes at the Windows command line, but this isn't quite
as flexible.

Solution
--------

This example script reads a Solexa/Illumina FASTQ from stdin, converts
the data to Sanger FASTQ (using PHRED scores) and writes it to stdout.

``` python
import sys
from Bio import SeqIO

recs = SeqIO.parse(sys.stdin, "fastq-solexa")
SeqIO.write(recs, sys.stdout, "fastq")
```

The following bash command can be used to extract the compressed
sequence and pipe it to the script (solexa2sanger\_fq.py).

<bash> gunzip -c some\_solexa.fastq.gz | python solexa2sanger\_fq.py
</bash>

This will write the sequence in Sanger FASTQ format to stdout - in this
case the screen.

You might can also use it like this, where the python script is fed its
data from an input file, and the output which would have been printed to
screen is instead sent to an output file:

<bash> python solexa2sanger\_fq.py &lt; some\_solexa.fastq &gt;
some\_phred.fastq </bash>
