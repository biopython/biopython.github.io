---
title: Getting Started
permalink: wiki/Getting_Started
layout: wiki
tags:
 - Documentation
 - TODO
---

--[Jblucks](User%3AJblucks "wikilink") 16:36, 19 January 2007 (EST):
This page modeled after the one for bioperl

Download
--------

See [Getting BioPerl](Getting_BioPerl "wikilink")

Installation
------------

See [Installing BioPerl](Installing_BioPerl "wikilink")

Quick example
-------------

Executing this:

``` python
from Bio.Seq import Seq,translate

#create a sequence object of some DNA
my_seq = Seq('CATGTAGATAG')

#print out some details about it
print 'seq is %i bases long' % len(my_seq)
print 'reverse complement is %s' % my_seq.reverse_complement().tostring()

#or see the whole record
print 'sequence record:', my_seq

#translate the sequence into a protein
my_protein = translate(my_seq)

print 'protein translation is %s' % my_protein.tostring()
print 'protein record:', my_protein
```

Produces:

    seq is 11 bases long
    reverse complement is CTATCTACATG
    sequence record: Seq('CATGTAGATAG', Alphabet())
    protein translation is HVD
    protein record: Seq('HVD', HasStopCodon(IUPACProtein(), '*'))

Beginners
---------

-   Learn how to program in [Perl](Perl "wikilink"), see several
    [Tutorials](Tutorials "wikilink")
    -   [Genome Informatics](http://stein.cshl.org/genome_informatics/)
    -   [PerlMonks
        Tutorials](http://www.perlmonks.org/index.pl?node=Tutorials)
    -   [Perl.com
        tutorials](http://www.perl.com/cs/user/query/q/6?id_topic=74)
    -   [Picking up Perl](http://www.ebb.org/PickingUpPerl/)
    -   [Learn Perl](http://learn.perl.org/)
-   Read the
-   Browse the [Bioperl Tutorial](Bptutorial "wikilink")
-   Examine the [Class Diagram](Class_Diagram "wikilink") if you'd like
    to know more about the relationships between the modules.

Further reading
---------------

-   Read the other [HOWTOs](HOWTOs "wikilink")
-   Use the Wiki Search tools to find more information on
    specific topics.

