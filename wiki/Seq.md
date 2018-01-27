---
title: Handling sequences with the Seq class.
permalink: wiki/Seq
layout: wiki
tags:
 - Wiki Documentation
---

In Biopython, sequences are usually held as ` Seq` objects, which hold
the sequence string and an associated alphabet.

This page describes the Biopython `Seq` object, defined in the `Bio.Seq`
module (together with related objects like the `MutableSeq`, plus some
general purpose sequence functions). In addition to this wiki page,
there is a whole chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) on the
`Seq` object - plus its [API
documentation](http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html)
(which you can read online, or from within Python with the help
command).

If you need to store additional information like a sequence identifier
or name, or even more details like a description or annotation, then we
use a [`SeqRecord`](SeqRecord "wikilink") object instead. These are the
sequence records used by the [`SeqIO`](SeqIO "wikilink") module for
reading and writing sequence files.

The Seq Object
==============

The `Seq` object essentially combines a Python string with an (optional)
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
and threonines. If *you* know, you should supply this information:

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
protein, would you?

``` python
>>> my_protein + my_dna
Traceback (most recent call last):
...
TypeError: Incompatable alphabets ProteinAlphabet() and DNAAlphabet()
```

Biopython will also catch things like trying to use nucleotide only
methods like translation (see below) on a protein sequence.

General methods
---------------

The `Seq` object has a number of methods which act just like those of a
Python string, for example the find method:

``` python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import generic_dna
>>> my_dna = Seq("AGTACACTGGT", generic_dna)
>>> my_dna
Seq('AGTACACTGGT', DNAAlphabet())
>>> my_dna.find("ACT")
5
>>> my_dna.find("TAG")
-1
```

There is a count method too:

``` python
>>> my_dna.count("A")
3
>>> my_dna.count("ACT")
1
```

However, watch out because just like the Python string's count, this is
a *non-overlapping* count!

``` python
>>> "AAAA".count("AA")
2
>>> Seq("AAAA", generic_dna).count("AA")
2
```

In some biological situations, you might prefer an overlapping count
which would give three for this example.

Nucleotide methods
------------------

If you have a nucleotide sequence (or a sequence with a generic
alphabet) you may want to do things like take the reverse complement, or
do a translation. Note some of these methods described here are only
available in Biopython 1.49 onwards.

### Complement and reverse complement

These are very simple - the methods return a new `Seq` object with the
appropriate sequence and the same alphabet:

``` python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import generic_dna
>>> my_dna = Seq("AGTACACTGGT", generic_dna)
>>> my_dna
Seq('AGTACACTGGT', DNAAlphabet())
>>> my_dna.complement()
Seq('TCATGTGACCA', DNAAlphabet())
>>> my_dna.reverse_complement()
Seq('ACCAGTGTACT', DNAAlphabet())
```

### Transcription and back transcription

If you have a DNA sequence, you may want to turn it into RNA. In
bioinformatics we normally assume the DNA is the coding strand (not the
template strand) so this is a simple matter of replacing all the
thymines with uracil:

``` python
>>> my_dna
Seq('AGTACACTGGT', DNAAlphabet())
>>> my_dna.transcribe()
Seq('AGUACACUGGU', RNAAlphabet())
```

Naturally, given some RNA, you might want the associated DNA - and again
Biopython does a simple U/T substitution:

``` python
>>> my_rna = my_dna.transcribe()
>>> my_rna
Seq('AGUACACUGGU', RNAAlphabet())
>>> my_rna.back_transcribe()
Seq('AGTACACTGGT', DNAAlphabet())
```

If you actually do want the template strand, you'd have to do a reverse
complement on top:

``` python
>>> my_rna
Seq('AGUACACUGGU', RNAAlphabet())
>>> my_rna.back_transcribe().reverse_complement()
Seq('ACCAGTGTACT', DNAAlphabet())
```

The chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) goes into
more detail on this strand issue.

### Translation

You can translate RNA:

``` python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import generic_rna
>>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", generic_rna)
>>> messenger_rna.translate()
Seq('MAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
```

Or DNA - which is assumed to be the coding strand:

``` python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import generic_dna
>>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)
>>> coding_dna.translate()
Seq('MAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
```

In either case there are several useful options - by default as you will
notice the in example above translation continues through any stop
codons, but this is optional:

``` python
>>> coding_dna.translate(to_stop=True)
Seq('MAIVMGR', ExtendedIUPACProtein())
```

Then there is the translation table, for which you can give an [NCBI
genetic code number or
name](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi):

``` python
>>> coding_dna.translate(table=2)
Seq('MAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
>>> coding_dna.translate(table="Vertebrate Mitochondrial")
Seq('MAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
```

You can of course combine these options:

``` python
>>> coding_dna.translate(table=2, to_stop=True)
Seq('MAIVMGRWKGAR', ExtendedIUPACProtein())
```

Consult the tutorial for more examples and arguments (e.g. specifying a
different symbol for a stop codon), or see the built in help:

``` python
>>> help(coding_dna.translate)
...
```

### Using nucleotide methods on a protein

None of this operations apply to a protein sequence and trying this will
raise an exception:

``` python
>>> my_protein.complement()
Traceback (most recent call last):
...
ValueError: Proteins do not have complements!
```

You can use them on `Seq` objects with a generic alphabet:

``` python
>>> my_seq
Seq('AGTACACTGGT', Alphabet())
>>> my_seq.complement()
Seq('TCATGTGACCA', Alphabet())
```
