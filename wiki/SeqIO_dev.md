---
title: SeqIO dev
permalink: wiki/SeqIO_dev
layout: wiki
---

This page is aimed at any developers or coders interesting in
understanding or extending the new Sequence Input/Output interface for
BioPython, [SeqIO](SeqIO "wikilink").

The code has now been checked into
[CVS](http://cvs.biopython.org/cgi-bin/viewcvs/viewcvs.cgi/biopython/Bio/SeqIO/?cvsroot=biopython#dirlist).
Related [Bug 2059](http://bugzilla.open-bio.org/show_bug.cgi?id=2059)
has been resolved.

The code is already available in BioPython 1.43.

Reading new file formats
------------------------

**Note:** The details are still subject to change

To add support for reading a new file format, you must implement an
iterator that expects a just file handle and returns SeqRecord objects.
You may do this using:

-   An iterator class subclassing something from Bio.SeqIO.Interfaces
-   A generator function (using the yield keyword; suitable for
    simple formats)
-   An ordinary function which returns an iterator. For example, you
    could build a list of SeqRecords and then turn it into an iterator
    using the iter() function.

You may accept additional *optional* arguments (an alphabet for
example). However there *must* be one and only one required argument
(the input file handle).

What you use as the SeqRecord's id, name and description will depend on
the file format. Ideally you would use the accesion number for the id.
This id should also be unique for each record (unless the records in the
file are in themselves ambiguous).

When storing any annotations in the record's annotations dictionary
follow the defacto standard laid down by the GenBank parser... I should
try and document this more.

If the supplied file seems to be invalid, raise a ValueError exception.

Finally, the new format must be added to the relevant dictionary mapping
in Bio/SeqIO/\_\_init\_\_.py so that the **Bio.SeqIO.parse()** and
**Bio.SeqIO.read()** functions are aware of it.

Writing new file formats
------------------------

**Note:** The details are still subject to change

To add support for writing a new file format you should write a sub
class of one of the writer objects in Bio.SeqIO.Interfaces

Then, the new format must be added to the relevant dictionary mappings
in Bio/SeqIO/\_\_init\_\_.py so that the **Bio.SeqIO.write()** function
is aware of your code.

If the supplied records cannot be written to this file format, raise a
ValueError exception. Where appropriate, please use the following
wording:

``` python
raise ValueError("Must have at least one sequence")
raise ValueError("Sequences must all be the same length")
raise ValueError("Duplicate record identifier: %s" % ...)
...
```

ToDo - Defined standard exceptions in Bio.SeqIO itself?

Possible additional formats
---------------------------

There are existing parsers in BioPython for the following file formats,
which could be integrated into Bio.SeqIO if appropriate.

### NBRF / PIR format

Bio.NBRF has a Martel parser for [PIR sequence
format](bp:PIR_sequence_format "wikilink"), which is similar to the
[FASTA format](bp:FASTA_sequence_format "wikilink"). It would need
addition work to return SeqRecords. It might be easier to extend to
reuse the Bio.SeqIO fasta code instead.

There is also [PSC
documentation](http://www.psc.edu/general/software/packages/seq-intro/nbrffile.html).

### KEGG format

Can Bio.KEGG parse files in [KEGG
format](bp:KEGG_sequence_format "wikilink")?

### MASE alignment format

Bio.IntelliGenetics seems to use Martel parse MASE format files into its
own record object. It could be extended to return SeqRecord objects.
See:

<http://pbil.univ-lyon1.fr/help/formats.html>

### MEME format

Bio.MEME has a parser for this file format, which at first glance looks
like it could be treated like an alignment format.

<http://meme.sdsc.edu>

### BLAST results

Pairwise alignments from the BLAST suite could be turned into two
SeqRecord objects with gapped sequences. Is this useful?

### COMPASS pairwise alignment format

Bio.Compass can parse the pairwise alignments from COMPASS. The output
is similar to BLAST in many ways. Again, is getting the results as
SeqRecord objects useful?
