---
title: AlignIO
permalink: wiki/AlignIO
layout: wiki
---

This page describes Bio.AlignIO, a new multiple sequence Alignment
Input/Output interface for BioPython which is currently only in our
source code repository. It should be included in the next release of
Biopython.

Aims
----

You may already be familiar with the [Bio.SeqIO](SeqIO "wikilink")
module which deals with files containing one or more sequences
represented as [SeqRecord](SeqRecord "wikilink") objects. The purpose of
the SeqIO module is to provide a simple uniform interface to assorted
file formats.

Similarly, Bio.AlignIO deals with files containing one or more sequence
alignments represented as Alignment objects. Bio.AlignIO uses the same
set of functions for input and output as in Bio.SeqIO, and the same
names for the file formats supported.

Note that the inclusion of Bio.AlignIO does lead to some duplication or
choice in how to deal with some file formats. For example, Bio.SeqIO and
Bio.Clustalw will both read sequences from Clustal files - but
Bio.Clustalw also includes a command line wrapper to call the program.

My vision is that for manipulating sequence alignments you should try
Bio.AlignIO as your first choice. In some cases you may only care about
the sequences themselves, in which case try using
[Bio.SeqIO](Bio.SeqIO "wikilink") on the alignment file directly. Unless
you have some very specific requirements, I hope this should suffice.

Peter

File Formats
------------

This table lists the file formats that Bio.SeqIO can read and write. The
format name is a simple lowercase string. Where possible we use the same
name as [BioPerl's
SeqIO](http://www.bioperl.org/wiki/HOWTO:SeqIO#Formats) and
[EMBOSS](http://emboss.sourceforge.net/docs/themes/SequenceFormats.html).

| Format name | Reads | Writes | Notes                                                                                                                                                                                              |
|-------------|-------|--------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| fasta       | Yes   | Yes    | This refers to the input file format introduced for Bill Pearson's FASTA tool, where each record starts with a "&gt;" line. Note that storing more than one alignment in this format is ambiguous. |
| clustal     | Yes   | Yes    | See also Bio.Clustalw for calling the command line tool.                                                                                                                                           |
| nexus       | Yes   | No     | Also known as PAUP format. Uses Bio.Nexus                                                                                                                                                          |
| phylip      | Yes   | Yes    | Truncates names at 10 characters.                                                                                                                                                                  |
| stockholm   | Yes   | Yes    | Also known as PFAM format.                                                                                                                                                                         |
||

In addition, you can store the (gapped) sequences from an alignment in
any of the [file formats supported by
Bio.SeqIO](SeqIO#File_Formats "wikilink"). The most common example of
this is storing multiple alignments in the simple fasta format. However,
storing more than one alignment in a single such file is ambiguous - see
the section below on alignment input.

Alignment Input
---------------

As in [Bio.SeqIO](SeqIO "wikilink"), there are two functions for
alignment input. These are **Bio.AlignIO.read()** for when the file
contains one and only one alignment, and the more general
**Bio.AlignIO.parse()** when the file may contain multiple separate
alignments.

Alignment Output
----------------

As in [Bio.SeqIO](SeqIO "wikilink"), there is a single output function
**Bio.AlignIO.write()**.

File Format Conversion
----------------------

Suppose you have a file containing PHYLIP alignment(s) that you want to
convert into the PFAM/Stockholm format:

``` python
from Bio import AlignIO

input_handle = open("example.phy", "rU")
output_handle = open("example.sth", "w")

alignments = AlignIO.parse(input_handle, "phylip")
AlignIO.write(alignments, output_handle, "stockholm")

output_handle.close()
input_handle.close()
```

By changing the format strings, that code could be used to convert
between any supported file formats.
