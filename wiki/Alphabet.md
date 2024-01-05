---
title: History and replacement of Bio.Alphabet
permalink: wiki/Alphabet
layout: wiki
tags:
 - Wiki Documentation
---

# Introduction

This page is intended to help people updating existing code using Biopython to
cope with the removal of the ``Bio.Alphabet`` module in Biopython 1.78
(September 2020).

The objects from `Bio.Alphabet` had two main uses:
* Record the molecule type (DNA, RNA or protein) of a sequence,
* Declare the expected characters in a sequence, alignment, motif, etc.

# Motivation

The intended purpose of the alphabet objects was never clearly defined, and
there were drawbacks to the twenty year old design. In particular the
`AlphabetEncoder` class (used to add information about gap or stop symbols)
was overly complicated, and made even determining the molecule type difficult.
Taking the consensus of multiple alphabet objects (for example during string
addition) was also complicated. While you could specify a strict alphabet like
unambiguous IUPAC DNA for a sequence, this did not enforce that only letters
used were `A`, `C`, `G` and `T`.

# Code changes

With no clear proposal on how to improve or replace the existing system, it
was agreed to remove it. In general you can simply remove any explicit use of
``Bio.Alphabet`` in your code.

## Seq changes

The ``Seq`` objects no longer have a `.alphabet` attribute, and do not do type
checking on `Seq` operations like adding protein to DNA anymore. Start by
removing the alphabet argument:

```python
# Old style
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

my_dna = Seq("ACGTTT", generic_dna)
```

```python
# New style
from Bio.Seq import Seq

my_dna = Seq("ACGTTT")
```

See below for where the alphabet was used to set the molecule type for an
output file format.

Another case where the alphabet was used was in declaring the gap character,
by default `-` in the various Biopython sequence and alignment parsers. If
you are using a different character, you will need to pass this to the `Seq`
object `.replace()` method explicitly now:

```python
# Old style
from Bio.Alphabet import generic_dna, Gapped
from Bio.Seq import Seq

my_dna = Seq("ACGT=TT", Gapped(generic_dna, "="))
print(my_dna.ungap())
```

```python
# New style
from Bio.Seq import Seq

my_dna = Seq("ACGT=TT")
print(my_dna.replace("=", ""))
```

## SeqRecord changes

Some sequence file formats require the molecule type when writing a file,
which previously was recorded with a ``Bio.Alphabet`` object as the
`.alphabet` attribute of the `Seq` object. This is now recorded as a molecule
type string in the `SeqRecord` object annotation dictionary instead.

```python
# Old style
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

seq = Seq("ATGCGTGCAT", generic_dna)
record = SeqRecord(seq, id="test")
SeqIO.write(record, "test_write.gb", "genbank")
```

```python
# New style
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

seq = Seq("ATGCGTGCAT")
record = SeqRecord(seq, id="test", annotations={"molecule_type": "DNA"})
SeqIO.write(record, "test_write.gb", "genbank")
```

```python
# Compatible with both pre- and post Biopython 1.78:
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

if generic_dna:
    # Newer Biopython refuses second argument
    seq = Seq("ATGCGTGCAT", generic_dna)
else:
    seq = Seq("ATGCGTGCAT")
record = SeqRecord(seq, id="test", annotations={"molecule_type": "DNA"})
SeqIO.write(record, "test_write.gb", "genbank")
```

The ``Bio.SeqIO`` parsing functions used to accept an optional alphabet
argument for setting this value if could not be determined from the
file format. That is no longer possible:

```python
# Old style
from Bio.Alphabet import generic_dna
from Bio import SeqIO

# This file has a single record only
record = SeqIO.read("Tests/Fasta/wisteria.nu", "fasta", generic_dna)
rec_start = record[:20]
SeqIO.write(rec_start, "start_only.xml", "seqxml")
```

In an example like this, you must explicitly set the molecule type in the
record annotation before writing it:

```python
# New style
from Bio import SeqIO

# This file has a single record only
record = SeqIO.read("Tests/Fasta/wisteria.nu", "fasta")
rec_start = record[:20]
rec_start.annotations["molecule_type"] = "DNA"
SeqIO.write(rec_start, "start_only.xml", "seqxml")
```

Similarly, the ``Bio.SeqIO.convert`` function's optional alphabet argument has
been replaced with an optional molecule type argument:

```python
# Old style
from Bio.Alphabet import generic_dna
from Bio import SeqIO

SeqIO.convert("example.fasta", "fasta", "example.xml", "seqxml", generic_dna)
```

```python
# New style
from Bio import SeqIO

SeqIO.convert("example.fasta", "fasta", "example.xml", "seqxml", "DNA")
```

This is one way to write a backward compatible version:

```python
# Compatible with both pre- and post Biopython 1.78:
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = "DNA"
from Bio import SeqIO

SeqIO.convert("example.fasta", "fasta", "example.xml", "seqxml", generic_dna)
```

## Other changes

Code which used an alphabet to specify the expected list of letters or symbols
now generally expects the valid characters as a string (e.g. ``Bio.motifs``).
