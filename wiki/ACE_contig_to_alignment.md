---
title: Represent an alignment from contig archived in ACE files.
permalink: wiki/ACE_contig_to_alignment
layout: wiki
tags:
 - Cookbook
---

Problem
-------

Sometimes it is useful to be able represent a contig produced as part of
genome or EST assembly as an alignment (eg to search for potential SNPs
in runs from mixed samples or to be able to write a contig out in a way
it can be viewed more easily). For assemblies that use the ACE file
format we can use Biopython's ACE handling to add reads that make a
contig to a generic alignment.'

Solution
--------

Let's represent contig in the ACE file that is used in Biopython's
testing framework:
[Tests/Ace/contig1.ace](https://github.com/biopython/biopython/blob/master/Tests/Ace/contig1.ace)
as an example

``` python
from Bio.Sequencing import Ace
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC, Gapped

def cut_ends(read, start, end):
  '''Replace residues on either end of a sequence with gaps.

  In this case we want to cut out the sections of each read which the assembler has
  decided are not good enough to include in the contig and replace them with gaps
  so the other information about positions in the read is maintained
  '''
  return (start-1) * '-' + read[start-1:end] + (len(read)-end) * '-'

def pad_read(read, start, conlength):
  ''' Pad out either end of a read so it fits into an alignment.

  The start argument is the position of the first base of the reads sequence in
  the contig it is part of. If the start value is negative (or 0 since ACE
  files count from 1, not 0) we need to take some sequence off the start
  otherwise each end is padded to the length of the consensus with gaps.
  '''
  if start < 1:
    seq = read[-1*start+1:]
  else:
    seq = (start-1) * '-' + read
  seq = seq + (conlength-len(seq)) * '-'
  return seq


# We will use the Ace parser to read individual contigs from file. Be aware
# that using this iterator can miss WA, CT, RT and WR tags (which can be
# anywhere in the file, e.g. the end). Read the file specification here:
# http://bozeman.mbt.washington.edu/consed/distributions/README.14.0.txt
# If you need these tags you'll need to use  Ace.read() (and lots of RAM).

ace_gen = Ace.parse(open("contig1.ace", 'r'))
contig = ace_gen.next()
align = Alignment(Gapped(IUPAC.ambiguous_dna, "-"))

# Now we have started our alignment we can add sequences to it,
# we will loop through contig's reads and get quality clipping from
# .reads[readnumber].qa and the position of each read in the contig
# .af[readnumber].padded_start and use the functions above to cut and
# pad the sequences before they are added

for readn in range(len(contig.reads)):
    clipst = contig.reads[readn].qa.qual_clipping_start
    clipe = contig.reads[readn].qa.qual_clipping_end
    start = contig.af[readn].padded_start
    seq = cut_ends(contig.reads[readn].rd.sequence, clipst, clipe)
    seq = pad_read(seq, start, len(contig.sequence))
    align.add_sequence("read%i" % (readn + 1), seq)
```

and when you print the alignment, or the sequences within it

``` python
>>>print align
Gapped(IUPACAmbiguousDNA(), '-') alignment with 2 rows and 856 columns
--------------------------------------------...--- read1
------GGATTGCCCTagtaacGGCGAGTGAAGCGGCAACAGCT...--- read2

>>> for read in align:
...    print read.seq[80:159]
tt*gtagagggaTGCTTCTGGGTAGCGACCGGTCTAAGTTCCTCGGAACAGGACGTCATAGAGGGTGAGAATCCCGTAT
TTTGTAGAGG*ATGCTTCTGGGTAGCGACCGGTCTAAGTTCCTCGGAACAGGACGTCATAGAGGGTGAGAATCCCGTAT
```

You can also now write the alignment to any format you want to using
[AlignIO](AlignIO "wikilink").

Discussion
----------

The details are given in the comments above, in broad strokes the ACE
contig is read in with Ace.parse(), a generic alignment is started then
the reads from the contig are added to the new alignment.
