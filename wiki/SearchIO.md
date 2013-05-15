---
title: SearchIO
permalink: wiki/SearchIO
layout: wiki
---

Matching the names in BioPerl, Biopython has a [SeqIO](SeqIO "wikilink")
module for sequence file input/output, and [AlignIO](AlignIO "wikilink")
for multiple sequence alignment input/output. The third member of the
BioPerl trio is SearchIO, and a Biopython equivalent is being worked on
in summer 2012 by [Google Summer of
Code](Google_Summer_of_Code "wikilink") student Wibowo Arindrarto
([blog](http://bow.web.id/blog/tag/gsoc/)).

The final module name in Biopython isn't settled, but Bio.SearchIO is
being used initially. This will cover pairwise sequence search file
input/output, for example from BLAST, HMMER or Bill Pearson's FASTA
suite. See the [BioPerl SearchIO
HOWTO](http://www.bioperl.org/wiki/HOWTO:SearchIO) for background.

This was included in Biopython 1.61 onwards as an *experimental* module
if you want to test it.
