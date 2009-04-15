---
title: SeqRecord
permalink: wiki/SeqRecord
layout: wiki
---

This page describes the **SeqRecord** object used in BioPython to hold a
sequence (as a [Seq](Seq "wikilink") object) with identifiers (ID and
name), description and optionally annotation and sub-features.

Most of the sequence file format parsers in BioPython can return
**SeqRecord** objects (and may offer a format specific record object
too, see for example Bio.SwissProt). The [SeqIO](SeqIO "wikilink")
system will *only* return SeqRecord objects.

In addition to the **SeqRecord** object's [API
documentation](http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html),
there is more information in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)), and the
[SeqIO](SeqIO "wikilink") page is also very relevant.

Creating a SeqRecord object
---------------------------

Most of the time you'll create **SeqRecord** objects by parsing a
sequence file with [Bio.SeqIO](SeqIO "wikilink"). However, it is useful
to know how to create a **SeqRecord** directly. For example,

``` python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
                   IUPAC.protein),
                   id="YP_025292.1", name="HokC",
                   description="toxic membrane protein, small")
print record
```

This would give the following output:

`   ID: YP_025292.1`  
`   Name: HokC`  
`   Description: toxic membrane protein, small`  
`   Number of features: 0`  
`   Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())`

Extracting information from a SeqRecord
---------------------------------------

Lets look in closer detail at the well annotated **SeqRecord** objects
Biopython creates from a GenBank file, such as
[ls\_orchid.gbk](http://biopython.org/DIST/docs/tutorial/examples/ls_orchid.gbk),
which we'll load using the [SeqIO](SeqIO "wikilink") module. This file
contains 94 records:

``` python
from Bio import SeqIO
for index, record in enumerate(SeqIO.parse(open("ls_orchid.gbk"), "genbank")) :
    print "index %i, ID = %s, length %i, with %i features" \
          % (index, record.id, len(record.seq), len(record.features))
```

And this is some of the output. Remember python likes to count from
zero, so the 94 records in this file have been labelled 0 to 93:

`index 0, ID = Z78533.1, length 740, with 5 features`  
`index 1, ID = Z78532.1, length 753, with 5 features`  
`index 2, ID = Z78531.1, length 748, with 5 features`  
`...`  
`index 92, ID = Z78440.1, length 744, with 5 features`  
`index 93, ID = Z78439.1, length 592, with 5 features`

Lets look in a little more detail at the final record:

``` python
print record
```

That should give you a hint of the sort of information held in this
object:

`ID: Z78439.1`  
`Name: Z78439`  
`Desription: P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA.`  
`Number of features: 5`  
`/source=Paphiopedilum barbatum`  
`/taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', ..., 'Paphiopedilum']`  
`/keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']`  
`/references=[`<Bio.SeqFeature.Reference ...>`, `<Bio.SeqFeature.Reference ...>`]`  
`/data_file_division=PLN`  
`/date=30-NOV-2006`  
`/organism=Paphiopedilum barbatum`  
`/gi=2765564`  
`Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACTTTGGTC ...', IUPACAmbiguousDNA())`

Lets look a little more closely... and use python's **dir()** function
to find out more about the SeqRecord object and what it does:

``` python
>>> dir(record)
[..., 'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'name', 'seq']
```

If you didn't already know, the **dir()** function returns a list of all
the methods and properties of an object (as strings). Those starting
with underscores in their name are "special" and we'll be ignoring them
in this discussion. We'll start with the **seq** property:

``` python
>>> print record.seq
Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACTTTGGTC ...', IUPACAmbiguousDNA())
>>> print record.seq.__class__
Bio.Seq.Seq
```

This is a [Seq](Seq "wikilink") object, another important object type in
Biopython, and worth of its own page on the wiki documentation.

The following three properties are all simple strings:

``` python
>>> print record.id
Z78439.1
>>> print record.name
Z78439
>>> print record.description
P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA.
```

Have a look at the raw GenBank file to see where these came from.

Next, we'll check the **dxrefs** property, which holds any database
cross references:

``` python
>>> print record.dbxrefs
[]
>>> print record.dbxrefs.__class__
<type 'list'>
```

An empty list? Disappointing... if we'd used a more recent GenBank file
the genome sequencing project reference would show up here.

How about the **annotations** property? This is a python dictionary...

``` python
>>> print record.annotations
{'source': 'Paphiopedilum barbatum', 'taxonomy': ...}
>>> print record.annotations.__class__
<type 'dict'>
>>> print record.annotations["source"]
Paphiopedilum barbatum
```

In this case, most of the values in the dictionary are simple strings,
but this isn't always the case - have a look at the references entry for
this example - its a list of **Reference** objects:

``` python
>>> print record.annotations["references"].__class__
<type 'list'>
>>> print len(record.annotations["references"])
2
>>> for ref in record.annotations["references"] : print ref.authors
Cox,A.V., Pridgeon,A.M., Albert,V.A. and Chase,M.W.
Cox,A.V.
```

Next is **features** which is another list property, and it contains
SeqFeature objects:

``` python
>>> print record.features.__class__
<type 'list'>
>>> print len(record.features)
5
```

SeqFeature objects are complicated enough to warrant their own page...

If you are using Biopython 1.48 or later, there will be a **format**
method. This lets you convert the **SeqRecord** into a string using one
of the output formats supported by [Bio.SeqIO](SeqIO "wikilink"), for
example:

``` python
>>> print record.format("fasta")
```

This should give:

`>Z78439.1 P.barbatum 5.8S rRNA gene and ITS1 and ITS2 DNA.`  
`CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACTTTGGTC`  
`ACCCATGGGCATTTGCTGTTGAAGTGACCTAGATTTGCCATCGAGCCTCCTTGGGAGCTT`  
`TCTTGTTGGCGAGATCTAAACCCCTGCCCGGCGGAGTTGGGCGCCAAGTCATATGACACA`  
`TAATTGGTGAAGGGGGTGGTAATCCTGCCCTGACCCTCCCCAAATTATTTTTTTAACAAC`  
`TCTCAGCAACGGATATCTCGGCTCTTGCATCGATGAAGAACGCAGCGAAATGCGATAATG`  
`GTGTGAATTGCAGAATCCCGTGAACATCGAGTCTTTGAACGCAAGTTGCGCCCGAGGCCA`  
`TCAGGCCAAGGGCACGCCTGCCTGGGCATTGCGAGTCATATCTCTCCCTTAATGAGGCTG`  
`TCCATACATACTGTTCAGCCGGTGCGGATGTGAGTTTGGCCCCTTGTTCTTTGGTACGGG`  
`GGGTCTAAGAGCTGCATGGGCTTTGGATGGTCCTAAATACGGAAAGAGGTGGACGAACTA`  
`TGCTACAACAAAATTGTTGTGCAAATGCCCCGGTTGGCCGTTTAGTTGGGCC`

If you are using Biopython 1.50 or later, there will also be a
**letter\_annotations** property. Again this is a dictionary but for
per-letter-annotation such as sequence quality scores or secondary
structure predictions. This kind of information isn't found in GenBank
files, so in this case the dictionary is empty:

``` python
>>> print record.letter_annotations
{}
```

Have a look at FASTQ or QUAL files to see how quality scores are
represented. Stockholm (PFAM) alignment files also often include
per-letter-annotation.
