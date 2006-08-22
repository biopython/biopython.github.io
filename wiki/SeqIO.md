---
title: SeqIO
permalink: wiki/SeqIO
layout: wiki
---

This page is for the documentation for a proposed new Sequence
Input/Output interface for BioPython.

The code is available on [Bug
2059](http://bugzilla.open-bio.org/show_bug.cgi?id=2059) and is being
discussed on the [Development mailing
list](http://biopython.org/wiki/Mailing_lists).

We would like to recreate the simplicity of [BioPerl's
SeqIO](http://www.bioperl.org/wiki/HOWTO:SeqIO), and in the long term
its [impressive list of supported file
formats](http://www.bioperl.org/wiki/Sequence_formats).

As currently implemented, the BioPython code covers multiple alignment
file formats as well. Alignment specific handling may be required in the
future should the BioPython alignment object become capable of holding
more than just sequence level annotation. See also BioPerl's list of
[multiple alignment
formats](http://www.bioperl.org/wiki/Multiple_alignment_formats).

Peter

Helper Functions
----------------

There are four helper functions which all take a filename, and optional
format. Each sequence is returned as a SeqRecord object.

-   File2SequenceIterator, returns a iterator (low memory, forward
    access only)
-   File2SequenceList, returns a list of sequences (high memory,
    random access)
-   File2SequenceDict, returns a dictionary of sequences (high memory,
    random access by ID)
-   File2Alignment, returns an alignment object (for use with multiple
    sequence alignment file formats)

For sequential file formats (like Fasta, GenBank, EMBL etc) the file can
be read record by record, and thus the iterator interface can save a
significant amount of memory (RAM) which allows you to deal with very
large files.

For interlaced file formats (like Clustal/Clustalw or annotated
Stockholm files) the entire file must be read in one go. You may not
save much memory by using the SequenceIterator in this case, but it is
provided so that you shouldn't need to re-write your code if you change
the input file format.

Examples using the Helper Functions
-----------------------------------

To do...

SequenceIterator
----------------

...
