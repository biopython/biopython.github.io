---
title: Seq
permalink: wiki/Seq
layout: wiki
---

In Biopython, sequences are usually held as Seq objects, which hold the
sequence string and an associated alphabet. There is a whole chapter in
the [Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) on the Seq
object.

If you need to store additional information like a sequence identifer or
name, or even more details like a description or annotation, then we use
a [SeqRecord](SeqRecord "wikilink") object instead. These are the
sequence records used by the [SeqIO](SeqIO "wikilink") module for
reading and writing sequence files.
