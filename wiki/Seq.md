---
title: Seq
permalink: wiki/Seq
layout: wiki
---

In Biopython, sequences are usually held as a Seq object, which holds
the sequence string and an associated alphabet.

There is a whole chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html) on the
Seq object.

If you need to store additional information like a sequence identifer or
name, or even more details like a description or annotation, then we use
a [SeqRecord](SeqRecord "wikilink") object instead. These are the
sequence record used by the [SeqIO](SeqIO "wikilink") module for reading
and writing sequence files.
