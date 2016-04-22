---
title: Population genetics in Biopython
permalink: wiki/PopGen
layout: wiki
tags:
 - Wiki Documentation
---

--------------------------------

Biopython makes available Population Genetics functionality in the
Bio.PopGen module.

1.  Access to Genepop methods (exact tests for Hardy–Weinberg
    equilibrium, population differentiation, genotypic disequilibrium,
    F-statistics, null allele frequencies, allele size-based statistics
    for microsatellites and much more).
2.  Coalescent simulation via
    [SimCoal2](http://cmpg.unibe.ch/software/simcoal2/) and
    [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/).
3.  Selection detection via
    [fdist2](http://www.maths.bris.ac.uk/~mamab/)

Documentation on how to access Genepop using Bio.PopGen can be found
[here](PopGen_Genepop "wikilink").

Documentation for interacting with SimCoal2 and fdist2 can be found on
the PopGen Chapter in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)).

Other software
--------------

Python programmers might want to consider [PyPop](http://pypop.org)
which includes intra-population frequentist statistics and
[simuPOP](http://simupop.sourceforge.net/) for forward-time population
genetics.

Future developments
-------------------

Any notes (or discussion) about possible future developments should be
restricted to the [PopGen development](PopGen_dev "wikilink") page.
