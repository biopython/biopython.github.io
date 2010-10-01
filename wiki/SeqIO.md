---
title: SeqIO
permalink: wiki/SeqIO
layout: wiki
---

This page describes Bio.SeqIO, the standard Sequence Input/Output
interface for BioPython 1.43 and later. For implementation details, see
the [SeqIO development page](SeqIO_dev "wikilink").

There is a whole chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) on
Bio.SeqIO, and although there is some overlap it is well worth reading
in addition to this WIKI page. There is also the [API
documentation](http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html)
(which you can read online, or from within Python with the help
command).

Aims
----

Bio.SeqIO provides a simple uniform interface to input and output
assorted sequence file formats (including multiple sequence alignments),
but will *only* deal with sequences as [SeqRecord](SeqRecord "wikilink")
objects. There is a sister interface [Bio.AlignIO](AlignIO "wikilink")
for working directly with sequence alignment files as Alignment objects.

The design was partly inspired by the simplicity of [BioPerl's
SeqIO](bp:HOWTO:SeqIO "wikilink"). In the long term we hope to match
BioPerl's impressive list of supported [sequence file
formats](bp:Sequence_formats "wikilink") and [multiple alignment
formats](bp:Multiple_alignment_formats "wikilink").

Note that the inclusion of Bio.SeqIO (and
[Bio.AlignIO](AlignIO "wikilink")) in Biopython does lead to some
duplication or choice in how to deal with some file formats. For
example, Bio.Nexus will also read sequences from Nexus files - but
Bio.Nexus can also do much more, for example reading any phylogenetic
trees in a Nexus file.

My vision is that for manipulating sequence data you should try
Bio.SeqIO as your first choice. Unless you have some very specific
requirements, I hope this should suffice.

[Peter](User%3APeter "wikilink")

File Formats
------------

This table lists the file formats that Bio.SeqIO can read, write and
index, with the Biopython version where this was first supported (or
[git](git "wikilink") to indicate this is supported in our latest in
development code). The format name is a simple lowercase string. Where
possible we use the same name as [BioPerl's
SeqIO](bp:HOWTO:SeqIO#Formats "wikilink") and
[EMBOSS](http://emboss.sourceforge.net/docs/themes/SequenceFormats.html).

| Format name           | Read | Write       | Index | Notes                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|-----------------------|------|-------------|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ace                   | 1.47 | No          | 1.52  | Reads the contig sequences from an ACE assembly file. Uses Bio.Sequencing.Ace internally                                                                                                                                                                                                                                                                                                                                                                  |
| clustal               | 1.43 | 1.43        | No    | The [alignment format of Clustal X and Clustal W](bp:ClustalW_multiple_alignment_format "wikilink").                                                                                                                                                                                                                                                                                                                                                      |
| embl                  | 1.43 | 1.54        | 1.52  | The [EMBL flat file format](bp:EMBL_sequence_format "wikilink"). Uses Bio.GenBank internally.                                                                                                                                                                                                                                                                                                                                                             |
| fasta                 | 1.43 | 1.43        | 1.52  | This refers to the input [FASTA file format](bp:FASTA_sequence_format "wikilink") introduced for Bill Pearson's FASTA tool, where each record starts with a "&gt;" line. Resulting sequences have a generic alphabet by default.                                                                                                                                                                                                                          |
| fastq-sanger or fastq | 1.50 | 1.50        | 1.52  | [FASTQ files](http://en.wikipedia.org/wiki/FASTQ_format) are a bit like FASTA files but also include sequencing qualities. In Biopython, "fastq" (or the alias "fastq-sanger") refers to Sanger style FASTQ files which encode PHRED qualities using an ASCII offset of 33. See also the incompatible "fastq-solexa" and "fastq-illumina" variants.                                                                                                       |
| fastq-solexa          | 1.50 | 1.50        | 1.52  | [FASTQ files](http://en.wikipedia.org/wiki/FASTQ_format) are a bit like FASTA files but also include sequencing qualities. In Biopython, "fastq-solexa" refers to (early) Solexa/Illumina style FASTQ files which encode Solexa qualities using an ASCII offset of 64. See also what we call the "fastq-illumina" format.                                                                                                                                 |
| fastq-illumina        | 1.51 | 1.51        | 1.52  | [FASTQ files](http://en.wikipedia.org/wiki/FASTQ_format) are a bit like FASTA files but also include sequencing qualities. In Biopython, "fastq-illumina" refers to recent Solexa/Illumina style FASTQ files (from pipeline version 1.3+) which encode PHRED qualities using an ASCII offset of 64. For *good* quality reads, PHRED and Solexa scores are approximately equal, so the "fastq-solexa" and "fastq-illumina" variants are almost equivalent. |
| genbank or gb         | 1.43 | 1.48 / 1.51 | 1.52  | The [GenBank or GenPept flat file format](bp:GenBank_sequence_format "wikilink"). Uses Bio.GenBank internally for parsing. Biopython 1.48 to 1.50 wrote basic GenBank files with only minimal annotation, while 1.51 onwards will also write the features table (see [Bug 2294](http://bugzilla.open-bio.org/show_bug.cgi?id=2294)).                                                                                                                      |
| ig                    | 1.47 | No          | 1.52  | This refers to the IntelliGenetics file format, apparently the same as the [MASE alignment format](bp:Mase_multiple_alignment_format "wikilink").                                                                                                                                                                                                                                                                                                         |
| imgt                  | git  | git         | git   | This refers to the IMGT variant of the EMBL plain text file format.                                                                                                                                                                                                                                                                                                                                                                                       |
| nexus                 | 1.43 | 1.48        | No    | The [NEXUS multiple alignment](bp:NEXUS_multiple_alignment_format "wikilink") format, also known as PAUP format. Uses Bio.Nexus internally.                                                                                                                                                                                                                                                                                                               |
| phd                   | 1.46 | 1.52        | 1.52  | [PHD files](bp:PHD_sequence_format "wikilink") are output from PHRED, used by PHRAP and CONSED for input. Uses Bio.Sequencing.Phd internally.                                                                                                                                                                                                                                                                                                             |
| phylip                | 1.43 | 1.43        | No    | An alignment format. Truncates names at 10 characters.                                                                                                                                                                                                                                                                                                                                                                                                    |
| pir                   | 1.48 | No          | 1.52  | A "FASTA like" format introduced by the National Biomedical Research Foundation (NBRF) for the [Protein Information Resource (PIR)](bp:PIR_sequence_format "wikilink") database, now part of UniProt.                                                                                                                                                                                                                                                     |
| stockholm             | 1.43 | 1.43        | No    | The [Stockholm alignment format](bp:Stockholm_multiple_alignment_format "wikilink") is also known as PFAM format.                                                                                                                                                                                                                                                                                                                                         |
| swiss                 | 1.43 | No          | 1.52  | [Swiss-Prot](bp:Swissprot_sequence_format "wikilink") aka UniProt format. Uses Bio.SwissProt internally.                                                                                                                                                                                                                                                                                                                                                  |
| tab                   | 1.48 | 1.48        | 1.52  | [Simple two column tab separated sequence files](bp:Tab_sequence_format "wikilink"), where each line holds a record's identifier and sequence. For example, this is used by Aligent's eArray software when saving microarray probes in a minimal tab delimited text file.                                                                                                                                                                                 |
| qual                  | 1.50 | 1.50        | 1.52  | [Qual files](bp:Qual_sequence_format "wikilink") are a bit like FASTA files but instead of the sequence, record space separated integer sequencing values as PHRED quality scores. A matched pair of FASTA and QUAL files are often used as an alternative to a single FASTQ file.                                                                                                                                                                        |
||

With Bio.SeqIO you can treat sequence alignment file formats just like
any other sequence file, but the new [Bio.AlignIO](AlignIO "wikilink")
module is designed to work with such alignment files directly. You can
also convert a set of [SeqRecord](SeqRecord "wikilink") objects from any
file format into an alignment - provided they are all the same length.
Note that when using Bio.SeqIO to write sequences to an alignment file
format, all the (gapped) sequences should be the same length.

Sequence Input
--------------

The main function is **Bio.SeqIO.parse()** which takes a file handle and
format name, and returns a [SeqRecord](SeqRecord "wikilink") iterator.
This lets you do things like:

``` python
from Bio import SeqIO
handle = open("example.fasta", "rU")
for record in SeqIO.parse(handle, "fasta") :
    print record.id
handle.close()
```

In the above example, we opened the file using the built-in python
function **open**. The argument 'rU' means open for **r**eading using
**u**niversal readline mode - this means you don't have to worry if the
file uses Unix, Mac or DOS/Windows style newline characters.

Note that you *must* specify the file format explicitly, unlike
[BioPerl's SeqIO](bp:HOWTO:SeqIO "wikilink") which can try to guess
using the file name extension and/or the file contents. See *Explicit is
better than implicit* ([The Zen of
Python](http://www.python.org/dev/peps/pep-0020/)).

If you had a different type of file, for example a Clustalw alignment
file such as
'[opuntia.aln](http://biopython.open-bio.org/SRC/biopython/Tests/Clustalw/opuntia.aln)'
which contains seven sequences, the only difference is you specify
"clustal" instead of "fasta":

``` python
from Bio import SeqIO
handle = open("opuntia.aln", "rU")
for record in SeqIO.parse(handle, "clustal") :
    print record.id
handle.close()
```

Iterators are great for when you only need the records one by one, in
the order found in the file. For some tasks you may need to have random
access to the records in any order. In this situation, use the built in
python **list** function to turn the iterator into a list:

``` python
from Bio import SeqIO
handle = open("example.fasta", "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
print records[0].id  #first record
print records[-1].id #last record
```

Another common task is to index your records by some identifier. For
small files we have a function **Bio.SeqIO.to\_dict()** to turn a
[SeqRecord](SeqRecord "wikilink") iterator (or list) into a dictionary
(in memory):

``` python
from Bio import SeqIO
handle = open("example.fasta", "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()
print record_dict["gi:12345678"] #use any record ID
```

The function **Bio.SeqIO.to\_dict()** will use the record ID as the
dictionary key by default, but you can specify any mapping you like with
its optional argument, **key\_function**.

For larger files, it isn't possible to hold everything in memory, so
**Bio.SeqIO.to\_dict()** is not suitable. Biopython 1.52 inwards
includes the **Bio.SeqIO.index()** function for this situation, but you
might also consider [BioSQL](BioSQL "wikilink").

``` python
from Bio import SeqIO
record_dict = SeqIO.index("example.fasta", "fasta")
print record_dict["gi:12345678"] #use any record ID
```

Biopython 1.45 introduced another function, **Bio.SeqIO.read()**, which
like **Bio.SeqIO.parse()** will expect a handle and format. It is for
use when the handle contains one and only one record, which is returned
as a single [SeqRecord](SeqRecord "wikilink") object. If there are no
records, or more than one, then an exception is raised:

``` python
from Bio import SeqIO
record = SeqIO.read(open("single.fasta"), "fasta")
```

For the related situation where you just want the first record (and are
happy to ignore any subsequent records), you can use the iterator's
**.next()** method:

``` python
from Bio import SeqIO
first_record = SeqIO.parse(open("example.fasta", "rU"), "fasta").next()
```

Sequence Output
---------------

For writing records to a file use the function **Bio.SeqIO.write()**,
which takes a SeqRecord iterator (or list), output handle and format
string:

``` python
from Bio import SeqIO
sequences = ... # add code here
output_handle = open("example.fasta", "w")
SeqIO.write(sequences, output_handle, "fasta")
output_handle.close()
```

There are more examples in the following section on converting between
file formats.

Note that if you are writing to an alignment file format, all your
sequences must be the same length.

If you supply the sequences as a **SeqRecord** iterator, then for
sequential file formats like Fasta or GenBank, the records can be
written one by one. Because only one record is created at a time, very
little memory is required. See the example below filtering a set of
records.

On the other hand, for interlaced or non-sequential file formats like
Clustal, the **Bio.SeqIO.write()** function will be forced to
automatically convert an iterator into a list. This will destroy any
potential memory saving from using an generator/iterator approach.

File Format Conversion
----------------------

Suppose you have a GenBank file which you want to turn into a Fasta
file. For example, lets consider the file
'[cor6\_6.gb](http://biopython.open-bio.org/SRC/biopython/Tests/GenBank/cor6_6.gb)'
which is included in the Biopython unit tests under the GenBank
directory.

You could read the file like this, using the **Bio.SeqIO.parse()**
function:

``` python
from Bio import SeqIO
input_handle = open("cor6_6.gb", "rU")
for record in SeqIO.parse(input_handle, "genbank") :
    print record
input_handle.close()
```

Notice that this file contains six records. Now instead of printing the
records, let's pass the SeqRecord iterator to the **Bio.SeqIO.write()**
function, to turn this GenBank file into a Fasta file:

``` python
from Bio import SeqIO

input_handle = open("cor6_6.gb", "rU")
output_handle = open("cor6_6.fasta", "w")

sequences = SeqIO.parse(input_handle, "genbank")
count = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()

print "Converted %i records" % count
```

Or more concisely using the **Bio.SeqIO.convert()** function (in
Biopython 1.52 or later), just:

``` python
from Bio import SeqIO
count = SeqIO.convert("cor6_6.gb", "genbank", "cor6_6.fasta", "fasta")
print "Converted %i records" % count
```

In this example the GenBank file started like this:

`LOCUS       ATCOR66M      513 bp    mRNA            PLN       02-MAR-1992`  
`DEFINITION  A.thaliana cor6.6 mRNA.`  
`ACCESSION   X55053`  
`VERSION     X55053.1  GI:16229`  
`...`

The resulting Fasta file looks like this:

`>X55053.1 A.thaliana cor6.6 mRNA.`  
`AACAAAACACACATCAAAAACGATTTTACAAGAAAAAAATA...`  
`...`

Note that all the Fasta file can store is the identifier, description
and sequence.

By changing the format strings, that code could be used to convert
between any supported file formats.

Examples
--------

### Input/Output Example - Filtering by sequence length

While you may simply want to convert a file (as shown above), a more
realistic example is to manipulate or filter the data in some way.

For example, let's save all the "short" sequences of less than 300
nucleotides to a Fasta file:

``` python
from Bio import SeqIO

short_sequences = [] # Setup an empty list
for record in SeqIO.parse(open("cor6_6.gb", "rU"), "genbank") :
    if len(record.seq) < 300 :
        # Add this record to our list
        short_sequences.append(record)

print "Found %i short sequences" % len(short_sequences)

output_handle = open("short_seqs.fasta", "w")
SeqIO.write(short_sequences, output_handle, "fasta")
output_handle.close()
```

If you know about **list comprehensions** then you could have written
the above example like this instead:

``` python
from Bio import SeqIO

input_seq_iterator = SeqIO.parse(open("cor6_6.gb", "rU"), "genbank")

#Build a list of short sequences:
short_sequences = [record for record in input_seq_iterator \
                   if len(record.seq) < 300]

print "Found %i short sequences" % len(short_sequences)

output_handle = open("short_seqs.fasta", "w")
SeqIO.write(short_sequences, output_handle, "fasta")
output_handle.close()
```

I'm not convinced this is actually any easier to understand, but it is
shorter.

However, if you are using Python 2.4 or later, and you are dealing with
very large files with thousands of records, you could benefit from using
a **generator expression** instead. This avoids creating the entire list
of desired records in memory:

``` python
from Bio import SeqIO

input_seq_iterator = SeqIO.parse(open("cor6_6.gb", "rU"), "genbank")
short_seq_iterator = (record for record in input_seq_iterator \
                      if len(record.seq) < 300)

output_handle = open("short_seqs.fasta", "w")
SeqIO.write(short_seq_iterator, output_handle, "fasta")
output_handle.close()
```

Remember that for sequential file formats like Fasta or GenBank,
**Bio.SeqIO.write()** will accept a **SeqRecord** iterator. The
advantage of the code above is that only one record will be in memory at
any one time.

However, as explained in the output section, for non-sequential file
formats like Clustal **write** is forced to automatically turn the
iterator into a list, so this advantage is lost.

If this is all confusing, *don't panic* and just ignore the fancy stuff.
For moderately sized datasets having too many records in memory at once
(e.g. in lists) is probably not going to be a problem.

### Using the SEGUID checksum

In this example, we'll use **Bio.SeqIO** with the
**Bio.SeqUtils.CheckSum** module (in Biopython 1.44 or later). First of
all, we'll just print out the checksum for each sequence in the GenBank
file
[ls\_orchid.gbk](http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk):

``` python
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
for record in SeqIO.parse(open("ls_orchid.gbk"), "genbank") :
    print record.id, seguid(record.seq)
```

You should get this output:

`Z78533.1 JUEoWn6DPhgZ9nAyowsgtoD9TTo`  
`Z78532.1 MN/s0q9zDoCVEEc+k/IFwCNF2pY`  
`...`  
`Z78439.1 H+JfaShya/4yyAj7IbMqgNkxdxQ`

Now lets use the checksum function and **Bio.SeqIO.to\_dict()** to build
a [SeqRecord](SeqRecord "wikilink") dictionary using the SEGUID as the
keys. The trick here is to use the python lambda syntax to create a
temporary function to get the SEGUID for each SeqRecord - we can't use
the **seguid** function directly as it only works on
[Seq](Seq "wikilink") objects or strings.

``` python
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
seguid_dict = SeqIO.to_dict(SeqIO.parse(open("ls_orchid.gbk"), "genbank"),
                            lambda rec : seguid(rec.seq))
record = seguid_dict["MN/s0q9zDoCVEEc+k/IFwCNF2pY"]
print record.id
print record.description
```

Giving this output:

`Z78439.1 `  
`P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA.`

### Random subsequences

This script will read a Genbank file with a whole mitochondrial genome
(e.g. the tobacco mitochondrion, *Nicotiana tabacum mitochondrion*
[NC\_006581](http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=NC_006581)),
create 500 records containing random fragments of this genome, and save
them as a fasta file. These subsequences are created using a random
starting points and a fixed length of 200.

``` python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

#There should be one and only one record, the entire genome:
mito_record = SeqIO.read(open("NC_006581.gbk"), "genbank")

mito_frags=[]
limit=len(mito_record.seq)
for i in range(0, 500) :
    start=randint(0,limit-200)
    end=start+200
    mito_frag=mito_record.seq[start:end]
    record=SeqRecord(mito_frag,'fragment_%i' % (i+1),'','')
    mito_frags.append(record)

output_handle = open("mitofrags.fasta", "w")
SeqIO.write(mito_frags, output_handle, "fasta")
output_handle.close()
```

That should give something like this as the output file,

`>fragment_1 `  
`TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAAT`  
`CCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTC`  
`TGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTT`  
`CTTTGGCACCTACGGTAGAG`  
`...`  
`>fragment_500 `  
`ACCCAGTGCCGCTACCCACTTCTACTAAGGCTGAGCTTAATAGGAGCAAGAGACTTGGAG`  
`GCAACAACCAGAATGAAATATTATTTAATCGTGGAAATGCCATGTCAGGCGCACCTATCA`  
`GAATCGGAACAGACCAATTACCAGATCCACCTATCATCGCCGGCATAACCATAAAAAAGA`  
`TCATTAAAAAAGCGTGAGCC`

### Writing to a string

Sometimes you won't want to write your [SeqRecord](SeqRecord "wikilink")
object(s) to a file, but to a string. For example, you might be
preparing output for display as part of a webpage. If you want to write
multiple records to a single string, use **StringIO** to create a
string-based handle. The
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) has an
example of this in the **SeqIO** chapter.

For the special case where you want a single record as a string in a
given file format, Biopython 1.48 added a new format method:

``` python
from Bio import SeqIO
mito_record = SeqIO.read(open("NC_006581.gbk"), "genbank")
print mito_record.format("fasta")
```

The format method will take any output format supported by **Bio.SeqIO**
where the file format can be used for a single record (e.g. "fasta",
"tab" or "genbank"). Note that we don't recommend you use this for file
output - using **Bio.SeqIO.write()** is faster and more general.

Help!
-----

If you are having problems with **Bio.SeqIO**, please join the
discussion mailing list (see [mailing lists](mailing_lists "wikilink")).

If you think you've found a bug, please report it on
[bugzilla](http://bugzilla.open-bio.org/).
