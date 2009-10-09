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

This example script uses [Bio.SeqIO](SeqIO "wikilink") to read a
Solexa/Illumina FASTQ from stdin, converts the data to Sanger FASTQ
(using PHRED scores) and writes it to stdout. See this more general page
on [converting sequence files](Converting_sequence_files "wikilink") for
some background.

``` python
import sys
from Bio import SeqIO
SeqIO.convert(sys.stdin, "fastq-solexa", sys.stdout, "fastq")
```

Pipes are a feature of the command line that enable the stdout output of
a program or command to be directed to the stdin input of another
command or program. For example, the following shell command can be used
to extract the compressed sequence and pipe it to the script
(solexa2sanger\_fq.py).

<bash> gunzip -c some\_solexa.fastq.gz | python solexa2sanger\_fq.py
</bash>

This will write the sequence in Sanger FASTQ format to stdout - in this
case the screen.

Redirection is similar to using pipes, but instead of directing the
output of one program to the input of another, redirection redirects the
contents of a file to a program's stdin, and/or the output of a
program's stdout to a file. In this example, the python script is fed
its data redirected from an input file, and the output which would have
been printed to screen is instead redirected to an output file:

<bash> python solexa2sanger\_fq.py &lt; some\_solexa.fastq &gt;
some\_phred.fastq </bash>

Redirection can also be used to redirect a program or command's stderr
to a file. Further examples of using redirection can be found
[here](http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO-3.html).
