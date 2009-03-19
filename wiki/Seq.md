---
title: Seq
permalink: wiki/Seq
layout: wiki
---

The Bio.Seq module in Biopython provides a few sequence related classes,
the Seq object and the MutableSeq object, plus some general purpose
sequence functions. In addition to this wiki page, there is a whole
chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) on this
topic.

In Biopython, sequences are usually held as Seq objects, which hold the
sequence string and an associated alphabet.

If you need to store additional information like a sequence identifer or
name, or even more details like a description or annotation, then we use
a [SeqRecord](SeqRecord "wikilink") object instead. These are the
sequence records used by the [SeqIO](SeqIO "wikilink") module for
reading and writing sequence files.

The Seq Object
==============

The Seq object essentially combines a Python string with an (optional)
biological alphabet. For example:

``` python
>>> from Bio.Seq import Seq
>>> my_seq = Seq("AGTACACTGGT")
>>> my_seq
Seq('AGTACACTGGT', Alphabet())
>>> my_seq.alphabet
Alphabet()
```

In the above example, we haven't specified an alphabet so we end up with
a default generic alphabet. Biopython doesn't know if this is a
nucleotide sequence or a protein rich in alanines, glycines, cysteines
and threonines. If you know, you should supply this information:

``` python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import generic_dna, generic_protein
>>> my_seq = Seq("AGTACACTGGT")
>>> my_seq
Seq('AGTACACTGGT', Alphabet())
>>> my_dna = Seq("AGTACACTGGT", generic_dna)
>>> my_dna
Seq('AGTACACTGGT', DNAAlphabet())
>>> my_protein = Seq("AGTACACTGGT", generic_protein)
>>> my_protein
Seq('AGTACACTGGT', ProteinAlphabet())
```

Why is this important? Well it can catch some errors for you - you
wouldn't want to accidentally try and combine a DNA sequence with a
protein would you:

``` python
>>> my_protein + my_dna
Traceback (most recent call last):
...
TypeError: Incompatable alphabets ProteinAlphabet() and DNAAlphabet()
```

Biopython will also catch things like trying to use nucleotide only
methods like translation (see below) on a protein sequence.
