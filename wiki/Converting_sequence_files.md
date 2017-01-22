---
title: Converting sequence files with the Bio.SeqIO module.
permalink: wiki/Converting_sequence_files
layout: wiki
tags:
 - Cookbook
---

Problem
-------

Many bioinformatics tools take different input file formats, so there is
a common need to interconvert between sequence file formats. One useful
option is the commandline tool [seqret from
EMBOSS](http://emboss.sourceforge.net/apps/cvs/emboss/apps/seqret.html),
but here we'll show how to tackle this problem with
[Bio.SeqIO](SeqIO "wikilink").

Solution
--------

Suppose you have a GenBank file which you want to turn into a Fasta
file. For example, let's consider the file
[`cor6_6.gb`](https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/cor6_6.gb)
(which is included in the Biopython unit tests under the GenBank
directory):

``` python
from Bio import SeqIO

input_handle = open("cor6_6.gb", "rU")
output_handle = open("cor6_6.fasta", "w")

sequences = SeqIO.parse(input_handle, "genbank")
count = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()
print("Converted %i records" % count)
```

In this example the GenBank file contained six records and started like
this:

```
LOCUS       ATCOR66M      513 bp    mRNA            PLN       02-MAR-1992
DEFINITION  A.thaliana cor6.6 mRNA.
ACCESSION   X55053
VERSION     X55053.1  GI:16229
...
```

The resulting Fasta file also contains all six records and looks like
this:

```
>X55053.1 A.thaliana cor6.6 mRNA.
AACAAAACACACATCAAAAACGATTTTACAAGAAAAAAATA...
...
```

Note that all the Fasta file can store is the identifier, description
and sequence.

By changing the format strings, that code could be used to convert
between any supported file formats.

You don't have to work with file handles - see this example using [stdin
and stdout pipes with Bio.SeqIO](Reading_from_unix_pipes "wikilink").

How it works
------------

See the [Bio.SeqIO](SeqIO "wikilink") page.
