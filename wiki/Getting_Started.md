---
title: Getting Started
permalink: wiki/Getting_Started
layout: wiki
---

Download and Installation
-------------------------

For Windows we provide click-and-run installers. Most Linux
distributions will include an optional Biopython package (although this
may be out of date). Otherwise you typically download and uncompress the
archive, and install from source. See our [downloads
page](Download "wikilink") for details including the prerequisites.

You can check your installation has worked at the python prompt:

``` python
>>> import Bio
```

If that gives no error, you should be done. If you get something like
"ImportError: No module named Bio" something has gone wrong.

Tutorial
--------

The Biopython Tutorial and Cookbook
([HTML](http://biopython.org/DIST/docs/tutorial/Tutorial.html),
[PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)) contains
the bulk of our documentation. See
[Documentation](Documentation "wikilink") for more links.

Quick example
-------------

Try executing this in python:

``` python
from Bio.Seq import Seq

#create a sequence object
my_seq = Seq('CATGTAGACTAG')

#print out some details about it
print('seq %s is %i bases long' % (my_seq, len(my_seq)))
print('reverse complement is %s' % my_seq.reverse_complement())
print('protein translation is %s' % my_seq.translate())
```

You should get the following output:

    seq CATGTAGACTAG is 12 bases long
    reverse complement is CTAGTCTACATG
    protein translation is HVD*

This was a very quick demonstration of Biopython's [Seq](Seq "wikilink")
(sequence) object and some of its methods.

Reading and writing Sequence Files
----------------------------------

Use the [SeqIO](SeqIO "wikilink") module for reading or writing
sequences as [SeqRecord](SeqRecord "wikilink") objects. For multiple
sequence alignment files, you can alternatively use the
[AlignIO](AlignIO "wikilink") module.

Beginners
---------

-   Learn how to program in [Python](http://www.python.org)
    -   [A Byte of
        Python](http://swaroopch.info/text/Byte_of_Python:Main_Page)
    -   [Dive Into Python](http://www.diveintopython.org/toc/index.html)
    -   [Python Quick
        Reference](http://rgruet.free.fr/PQR25/PQR2.5.html)
-   Browse the [Biopython
    Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
-   Read this paper <biblio>Bassi2007 pmid=18052533</biblio>
-   Examine the [Class Diagram](http://biopython.org/DIST/docs/api) if
    you'd like to know more about the relationships between the modules.

Further reading
---------------

-   Use the Wiki Search tools to find more information on
    specific topics.

