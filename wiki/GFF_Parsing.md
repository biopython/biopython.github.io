---
title: Parsing GFF Files
permalink: wiki/GFF_Parsing
layout: wiki
tags:
 - Wiki Documentation
---

*Note:* GFF parsing is not yet integrated into Biopython. This
documentation is work towards making it ready for inclusion. You can
retrieve the current version of the GFF parser from:
<http://github.com/chapmanb/bcbb/tree/master/gff>, which in turn
led to <https://github.com/daler/gffutils>. Comments are very
welcome.

Generic Feature Format (GFF) is a biological sequence file format for
representing features and annotations on sequences. It is a tab
delimited format, making it accessible to biologists and editable in
text editors and spreadsheet programs. It is also well defined and can
be parsed via automated programs. GFF files are available from many of
the large sequencing and annotation centers. The
[specification](http://www.sequenceontology.org/gff3.shtml) provides
full details on the format and its uses.

Biopython provides a full featured GFF parser which will handle several
versions of GFF: GFF3, GFF2, and GTF. It supports writing GFF3, the
latest version.

GFF parsing differs from parsing other file formats like GenBank or PDB
in that it is not record oriented. In a GenBank file, sequences are
broken into discrete parts which can be parsed as a whole. In contrast,
GFF is a line oriented format with support for nesting features. GFF is
also commonly used to store only biological features, and not the
primary sequence.

These differences have some consequences in how you will deal with GFF:

-   Files are first examined to determine available annotations and
    define items of interest.
-   Sequences will often be parsed separately, from an associated
    FASTA file.
-   Parsing needs to consider available memory, which can be quickly
    used up on files with many annotations. This problem can be solved
    via two methods:
    -   Limiting parsing to features of interest.
    -   Iterating over portions of the file.

The documentation below provides a practical guide to examining, parsing
and writing GFF files in Python.

## Examining your GFF file

Since GFF is a very general format, it is extremely useful to start by
getting a sense of the type of data in the file and how it is
structured. `GFFExaminer` provides an interface to examine and query the
file. To examine relationships between features, examine a dictionary
mapping parent to child features:

``` python
import pprint
from BCBio.GFF import GFFExaminer

in_file = "your_file.gff"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
in_handle.close()
```

This file contains a flexible three level description of coding
sequences: genes have mRNA trasncripts; those mRNA transcripts each
contain common features of coding sequence, the CDS itself, exon, intron
and 5' and 3' untranslated regions. This is a common GFF structure
allowing representation of multiple transcripts:

``` bash
{('Coding_transcript', 'gene'): [('Coding_transcript', 'mRNA')],
 ('Coding_transcript', 'mRNA'): [('Coding_transcript', 'CDS'),
                                 ('Coding_transcript', 'exon'),
                                 ('Coding_transcript', 'five_prime_UTR'),
                                 ('Coding_transcript', 'intron'),
                                 ('Coding_transcript', 'three_prime_UTR')]}
```

Another item of interest for designing your parse strategy is
understanding the various tags used to label the features. These consist
of:

-   `gff_id` -- The record identifier being described. This will often
     refer to a chromosome or other scaffold sequence.
-   `gff_source` -- The source description from the second column of the
     GFF file, which specifies how a feature was generated.
-   `gff_type` -- The type of the feature, pulled from the 3rd column of
     the GFF file.
-   `gff_source_type` -- All combinations of sources and types in
     the file.

The `available_limits` function in the examiner gives you a high level
summary of these feature attributes, along with counts for the number of
times they appear in the file:

``` python
import pprint
from BCBio.GFF import GFFExaminer

in_file = "your_file.gff"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.available_limits(in_handle))
in_handle.close()
```

``` bash
{'gff_id': {('I',): 159,
            ('II',): 3,
            ('III',): 2,
            ('IV',): 5,
            ('V',): 2,
            ('X',): 6},
 'gff_source': {('Allele',): 1,
                ('Coding_transcript',): 102,
                ('Expr_profile',): 1,
                ('GenePair_STS',): 8,
                ('Oligo_set',): 1,
                ('Orfeome',): 8,
                ('Promoterome',): 5,
                ('SAGE_tag',): 1,
                ('SAGE_tag_most_three_prime',): 1,
                ('SAGE_tag_unambiguously_mapped',): 12,
                ('history',): 30,
                ('mass_spec_genome',): 7},
 'gff_source_type': {('Allele', 'SNP'): 1,
                     ('Coding_transcript', 'CDS'): 27,
                     ('Coding_transcript', 'exon'): 33,
                     ('Coding_transcript', 'five_prime_UTR'): 4,
                     ('Coding_transcript', 'gene'): 2,
                     ('Coding_transcript', 'intron'): 29,
                     ('Coding_transcript', 'mRNA'): 4,
                     ('Coding_transcript', 'three_prime_UTR'): 3,
                     ('Expr_profile', 'experimental_result_region'): 1,
                     ('GenePair_STS', 'PCR_product'): 8,
                     ('Oligo_set', 'reagent'): 1,
                     ('Orfeome', 'PCR_product'): 8,
                     ('Promoterome', 'PCR_product'): 5,
                     ('SAGE_tag', 'SAGE_tag'): 1,
                     ('SAGE_tag_most_three_prime', 'SAGE_tag'): 1,
                     ('SAGE_tag_unambiguously_mapped', 'SAGE_tag'): 12,
                     ('history', 'CDS'): 30,
                     ('mass_spec_genome', 'translated_nucleotide_match'): 7},
 'gff_type': {('CDS',): 57,
              ('PCR_product',): 21,
              ('SAGE_tag',): 14,
              ('SNP',): 1,
              ('exon',): 33,
              ('experimental_result_region',): 1,
              ('five_prime_UTR',): 4,
              ('gene',): 2,
              ('intron',): 29,
              ('mRNA',): 4,
              ('reagent',): 1,
              ('three_prime_UTR',): 3,
              ('translated_nucleotide_match',): 7}}
```

## GFF Parsing

### Basic GFF parsing

Generally, the GFF parser works similar to other parsers in Biopython.
Calling `parse` with a handle to a GFF file returns a set of `SeqRecord`
objects corresponding to the various IDs referenced in the file:

``` python
from BCBio import GFF

in_file = "your_file.gff"

in_handle = open(in_file)
for rec in GFF.parse(in_handle):
    print rec
in_handle.close()
```

The rec object is a Biopython [`SeqRecord`](SeqRecord "wikilink")
containing the features described in the GFF file. The features are
ordered into parent-child relationships based on the line by line
information in the original GFF file. See the detailed documentation on
[`SeqRecord`](SeqRecord "wikilink") and
[`SeqFeature`](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:seq_features)
objects for more details on accessing the information in these objects.

Since a GFF file is not broken down into an explicit record structure,
this requires reading the entire file, parsing all of the features, and
then returning those as records. This will be fine for small files, but
for most real life cases you will want to restrict parsing to a set of
features of interest or a section of lines at once to conserve memory.

### Limiting to features of interest

A GFF file will commonly contain many types of features, and you will be
interested in retrieving a subset of these. The `limit_info` argument to
`GFF.parse` allows exact specification of which features to parse, turn
into objects and retrieve. An example is retrieving all coding sequence
on chromosome 1:

``` python
from BCBio import GFF

in_file = "your_file.gff"

limit_info = dict(gff_id=["chr1"], gff_source=["Coding_transcript"])

in_handle = open(in_file)
for rec in GFF.parse(in_handle, limit_info=limit_info):
    print rec.features[0]
in_handle.close()
```

You will get back a single record for chromosome 1 which contains all of
the coding features in memory for further manipulation. Depending on
your memory requirements and workflow, it may make sense to do analysis
over each chromosome or set of features you are interested in.

### Iterating over portions of a file

Another way to break up a large GFF file parse into sections is to limit
the number of lines that are read at once. This is a useful workflow for
GFF files in which you don't need all of the features at once and can do
something useful with a few at a time. To do this, pass the
`target_lines` argument to `GFF.parse`:

``` python
from BCBio import GFF

in_file = "your_file.gff"

in_handle = open(in_file)
for rec in GFF.parse(in_handle, target_lines=1000):
    print rec
in_handle.close()
```

The parser will attempt to smartly break up the file at your requested
number of lines. For instance, if 1000 lines happens to come in the
middle of a nested coding feature (gene -&gt; transcript -&gt;
CDS/exon/intron), the parser would continue until the entire feature
region is read. This helps ensure that you have fully formed features
for analysis.

If your file has no nesting of features, or you just want a single line
at once, you can set `target_lines=1` and the parser will happily give
you back a `SeqRecord` object with a single `SeqFeature` for every line.

### Providing initial sequence records

GFF records normally contain annotation data, while sequence information
is available in a separate FASTA formatted file. The GFF parser can add
annotations to existing records. First parse the sequence file with
`SeqIO`, then feed the resulting sequence dictionary to the GFF parser:

``` python
from BCBio import GFF
from Bio import SeqIO

in_seq_file = "seqs.fa"
in_seq_handle = open(in_seq_file)
seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()

in_file = "your_file.gff"
in_handle = open(in_file)
for rec in GFF.parse(in_handle, base_dict=seq_dict):
    print rec
in_handle.close()
```

Note that this just adds directly to the existing dictionary. If you
apply filters to the GFF parser these are only applied to annotations;
records will not be removed from the initial sequence dictionary.

## Writing GFF3

The `GFF3Writer` takes an iterator of `SeqRecord objects`, and writes each
`SeqFeature` as a GFF3 line:

-   `seqid` -- SeqRecord ID
-   `source` -- Feature qualifier with key "source"
-   `type` -- Feature type attribute
-   `start`, `end` -- The Feature Location
-   `score` -- Feature qualifier with key "score"
-   `strand` -- Feature strand attribute
-   `phase` -- Feature qualifier with key "phase"

The remaining qualifiers are the final key/value pairs of the attribute.

A feature hierarchy is represented as `sub_features` of the parent
feature. This handles any arbitrarily deep nesting of parent and child
features.

### Converting other formats to GFF3

``` python
from BCBio import GFF
from Bio import SeqIO

in_file = "your_file.gb"
out_file = "your_file.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")

GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

in_handle.close()
out_handle.close()
```

### Writing GFF3 from scratch

You can create Biopython `SeqRecord` and `SeqFeature` objects from scratch
and use these to generate GFF output. [Chapter 4 of the
Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html) goes
into detail about the objects, and this example demonstrates the major
features:

``` python
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

out_file = "your_file.gff"
seq = Seq("GATCGATCGATCGATCGATC")
rec = SeqRecord(seq, "ID1")
qualifiers = {
    "source": "prediction",
    "score": 10.0,
    "other": ["Some", "annotations"],
    "ID": "gene1",
}
sub_qualifiers = {"source": "prediction"}
top_feature = SeqFeature(
    FeatureLocation(0, 20), type="gene", strand=1, qualifiers=qualifiers
)
top_feature.sub_features = [
    SeqFeature(FeatureLocation(0, 5), type="exon", strand=1, qualifiers=sub_qualifiers),
    SeqFeature(
        FeatureLocation(15, 20), type="exon", strand=1, qualifiers=sub_qualifiers
    ),
]
rec.features = [top_feature]

with open(out_file, "w") as out_handle:
    GFF.write([rec], out_handle)
```

This generates the following GFF:

``` bash
##gff-version 3
##sequence-region ID1 1 20
ID1     prediction      gene    1       20      10.0    +       .       other=Some,annotations;ID=gene1
ID1     prediction      exon    1       5       .       +       .       Parent=gene1
ID1     prediction      exon    16      20      .       +       .       Parent=gene1
```
