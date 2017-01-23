---
title: Developing statistics for the Population Genetics Module.
permalink: wiki/PopGen_dev_Statistics
layout: wiki
---

Introduction
------------

A few observations about population genetics and bioinformatics

1.  Population genetics is used in a wide variety of settings, from
    cancer studies to conservation of endangered species. Even names
    might not be standardized, e.g. a human geneticist might use *Short
    Tandem Repeat* (STR) while an animal geneticist might
    use *microsatellite*.
2.  Different research communities have completely different sets of
    requirements: markers change (SNPs, microsatellites/STSs, RFLPs,
    sequences); number of individuals sampled change (for few units, to
    thousands to very large synthetic datasets); number of populations
    change; the density of markers change (from 10 million SNPs in the
    HapMap for humans and full sequences appearing for many individuals
    on parasites, to the opposite, where only 20-30 loci are available)
3.  A good implementation should support the vast majority of scenarios
    above: should support as many markers as possible, big and small
    datasets
4.  The most used format in population genetics is still the
    **Genepop** format. This format is not used much outside the
    popgen community. For an overview of the importance of the format,
    read [Excoffier and Heckel (2006)](http://dx.doi.org/10.1038/nrg1904).
    This happens to be a marker-independent format.
5.  The popgen community while producing lots of free software cares
    very little about licensing or data-standardization issues.
    Different ad-hoc formats exist (again, see the above paper).

As such it is important to characterize all types of existing statistics
and to create a framework that accommodates all of them if people decide
to implement them in the future.

**Restricted use cases should be handled with extreme care**. If
possible people from different backgrounds should contribute.

Different dimensions
--------------------

Here we characterize the different dimensions that need to be accounted
for. For a good intro see page 98 of the [Arlequin3
manual](http://cmpg.unibe.ch/software/arlequin35/man/Arlequin35.pdf).

### Marker dependent versus marker independent statistics

Some statistics require a special kind of marker. For instance Tajima D
requires a sequence or a RFLP. Allelic range requires
microsatellites/STRs. Other statistics are marker independent. For
instance Fst relies on allele counts per population.

### Intra-Population versus Inter-population statistics

As an example observed heterozygosity has no notion of population
structure.

Fst is an example of a statistic that exists to measure population
differentiation. These kind of statistics require some notion of
population differentiation.

### Data type

Haplotypic, genotypic (phase unknown), genotypic (phase known),
genoptypic dominant, frequency only.

Say, for expected heterozygosity frequencies are enough, for observed
heterozygosity genotypic (phase unknown) data is necessary.

### Single locus versus multi-loci

Single locus as in allelic richness, ExpHe, Fst.

Multi-loci as in number of polimorphic sites, LD or EHH.

### Temporal/longitudinal vs single point in time

Say temporal-Fst versus Fst.

### Population versus Landscape

This issue I suggest abandon for now.

### Example of statistic classification

-   ExpHz non-temporal, intra, single-locus, marker independent,
    genotypic - gametic unk
-   ObsHz non-temporal, intra, single-locus, independent, genotypic -
    gametic kn
-   Fst(CW) non-temporal, inter, single-locus, indep, genotypic -
    gametic unk
-   temporal-Fst temporal, intra, single-locus, indep, genotypic -
    gametic unk
-   LD(D') non-temporal, intra, multi-locus, indep, haplo/geno
-   Fk temporal, intra, single-locus, indep, geno
-   S (polimorphic sites), non-temporal, intra, multi-locus, indep,
    haplo/geno
-   Alleic range, nt, intra, single-locus, microsat, haplo/geno
-   EHH, nt, positional
-   Tajima D, nt, intra, single-locus, sequence/RFLP

Design
------

The design will have to be able to cope with all the dimensions above


Pending issues
--------------

There is still the issue of statistical tests (say Hardy-Weinberg
deviation).
