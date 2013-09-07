---
title: SearchIO
permalink: wiki/SearchIO
layout: wiki
tags:
 - Wiki Documentation
---

Matching the names in BioPerl, Biopython has a [SeqIO](SeqIO "wikilink")
module for sequence file input/output, and [AlignIO](AlignIO "wikilink")
for multiple sequence alignment input/output. The third member of the
BioPerl trio is SearchIO, and a Biopython equivalent was written during
summer 2012 by [Google Summer of Code](Google_Summer_of_Code "wikilink")
student Wibowo Arindrarto ([blog](http://bow.web.id/blog/tag/gsoc/)).

This covers pairwise sequence search file input/output, for example from
BLAST, HMMER, BLAT, or Bill Pearson's FASTA suite. See the [BioPerl
SearchIO HOWTO](http://www.bioperl.org/wiki/HOWTO:SearchIO) for
background.

It is included in Biopython 1.61 onwards as an *experimental* module if
you want to test it. An chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) on
Bio.SearchIO is also published alongside the 1.61 release.

This wiki describes the important bits with some small examples. For a
full reference, consult the [API
documentation](http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html).

Supported File Formats
----------------------

The table below lists all formats supported by Bio.SearchIO. Note that
for writing support, the writer assumes that all the necessary
attributes of the objects being written are present. It is not possible,
for example, to write BLAST XML data to a HMMER 3.0 plain text output
straight away.

<table style="width:100%;">
<caption>Table 1: Bio.SearchIO supported file formats</caption>
<colgroup>
<col width="12%" />
<col width="22%" />
<col width="22%" />
<col width="22%" />
<col width="22%" />
</colgroup>
<thead>
<tr class="header">
<th><p>Format name</p></th>
<th><p>Read</p></th>
<th><p>Write</p></th>
<th><p>Index</p></th>
<th><p>Notes</p></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><p>blast-tab</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>BLAST+ tabular output (both <code>-m</code> <code>6</code> and <code>-m</code> <code>7</code> flags are supported).</p></td>
</tr>
<tr class="even">
<td><p>blast-text</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>n/a</p></td>
<td><p>BLAST+ plain text output (up to version 2.2.26+). Newer versions may not always work.</p></td>
</tr>
<tr class="odd">
<td><p>blast-xml</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>BLAST+ XML output.</p></td>
</tr>
<tr class="even">
<td><p>blat-psl</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>BLAT default output (PSL format). Variants with or without header are both supported. PSLX (PSL + sequences) is also supported.</p></td>
</tr>
<tr class="odd">
<td><p>exonerate-text</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>Exonerate plain text output. Due to the way Biopython stores its sequences, at the moment support is limited to text outputs without split codons (for protein queries). If you are parsing a text output of protein queries containing split codon alignments (for example, from the <code>protein2genome</code> alignment mode), the parser will fail.</p></td>
</tr>
<tr class="even">
<td><p>exonerate-cigar</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>Exonerate cigar string.</p></td>
</tr>
<tr class="odd">
<td><p>exonerate-vulgar</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>Exonerate vulgar string.</p></td>
</tr>
<tr class="even">
<td><p>fasta-m10</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>Bill Pearson's FASTA <code>-m</code> <code>10</code> output.</p></td>
</tr>
<tr class="odd">
<td><p>hmmer3-domtab</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>HMMER3.0 domain table output format. The name <code>hmmer3-domtab</code> per se is in fact not used, since the program name has to be specified. For example, when parsing hmmscan output, the format name would be <code>hmmscan-domtab</code>.</p></td>
</tr>
<tr class="even">
<td><p>hmmer3-tab</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>1.61</p></td>
<td><p>HMMER 3.0 table output format.</p></td>
</tr>
<tr class="odd">
<td><p>hmmer3-text</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>HMMER 3.0 plain text output format.</p></td>
</tr>
<tr class="even">
<td><p>hmmer2-text</p></td>
<td><p>1.61</p></td>
<td><p>n/a</p></td>
<td><p>1.61</p></td>
<td><p>HMMER 2.x plain text output format.</p></td>
</tr>
<tr class="odd">
</tr>
</tbody>
</table>

Format-specific Arguments
-------------------------

Although mostly similar to Biopython's SeqIO and AlignIO modules, there
is a small difference in the main Bio.SearchIO functions. Depending on
the file format being used, you may pass additional keyword arguments
that determines how the parser / indexer / writer behaves. Shown below
are some formats which accepts extra keyword arguments.

| Format name | Argument name                            | Default value              | Applicable for                                                                              | Explanation                                                             |
|-------------|------------------------------------------|----------------------------|---------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| blast-tab   | comments                                 | False                      | reading, writing, indexing                                                                  | Boolean, whether the input/output file is the commented variant or not. |
| fields      | Default BLAST tabular output field names | reading, writing, indexing | Space-separated string, list of fields / columns in the input/output file.                  |
| blast-xml   | encoding                                 | "utf-8"                    | writing                                                                                     | XML encoding name.                                                      |
| indent      | " " (empty space)                        | writing                    | Character(s) to use for indenting the XML.                                                  |
| increment   | 2                                        | writing                    | How many times the character defined in `indent` are printed when printing a child element. |
| blat-psl    | pslx                                     | False                      | reading, writing, indexing                                                                  | Boolean, whether the input/output file contains sequences or not.       |
| header      | False                                    | writing                    | Boolean, whether to write PSL header or not.                                                |
||

Conventions
-----------

The main goal of creating Bio.SearchIO is to have a common, easy to use
interface across different search output files. As such, we have also
created some conventions / standards for Bio.SearchIO that extend beyond
the common object model. These conventions apply to all files parsed by
Bio.SearchIO, regardless of their individual formats.

### Python-style sequence coordinates

When storing sequence coordinates (start and end values), Bio.SearchIO
uses the Python-style slice convention: zero-based and half-open
intervals. For example, if in a BLAST XML output file the start and end
coordinates of an HSP are 10 and 28, they would become 9 and 28 in
Bio.SearchIO. The start coordinate becomes 9 because Python indices
start from zero, while the end coordinate remains 28 as Python slices
omit the last item in an interval.

Beside giving you the benefits of standardization, this convention also
makes the coordinates usable for slicing sequences. For example, given a
full query sequence and the start and end coordinates of an HSP, one can
use the coordinates to extract part of the query sequence that results
in the database hit.

When these objects are written to an output file using
Bio.SearchIO.write, the coordinate values are restored to their
respective format's convention. Using the example above, if the HSP
would be written to an XML file, the start and end coordinates would
become 10 and 28 again.

### Sequence coordinate order

Some search output format reverses the start and end coordinate
sequences according to the sequence's strand. For example, in BLAST
plain text format if the matching strand lies in the minus orientation,
then the start coordinate will always be bigger than the end coordinate.

In Bio.SearchIO, start coordinates are always smaller than the end
coordinates, regardless of their originating strand. This ensures
consistency when using the coordinates to slice full sequences.

Note that this coordinate order convention is only enforced in the
HSPFragment level. If an HSP object has several HSPFragment objects,
each individual fragment will conform to this convention. But the order
of the fragments within the HSP object follows what the search output
file uses.

Similar to the coordinate style convention, the start and end
coordinates' order are restored to their respective formats when the
objects are written using Bio.SearchIO.write.

### Frames and strand values

Bio.SearchIO only allows -1, 0, 1 and None as strand values. For frames,
the only allowed values are integers from -3 to 3 (inclusive) and None.
Both of these are standard Biopython conventions.

FAQ
---

1.  *How does Bio.SearchIO differ from Bio.Blast.NCBIXML*?
    Both modules are based on completely different object models and are
    not compatible with each other. Not only that, the underlying
    parsers and writers are also different (indexing is not possible
    with Bio.Blast.NCBIXML). Finally, Bio.SearchIO is planned to be the
    replacement of Bio.Blast.NCBIXML.

<!-- -->

1.  *How does Bio.SearchIO differ from Bio.Blast.NCBIStandalone*?
    Again, they provide different object models. However, Bio.SearchIO
    currently uses the parser from Bio.Blast.NCBIStandalone internally,
    but that old module will be deprecated.


