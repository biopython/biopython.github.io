---
title: Concatenating multiple alignments NEXUS files with the Bio.Nexus module.
permalink: wiki/Concatenate_nexus
layout: wiki
tags:
 - Cookbook
---

The Problem
-----------

It is common to make species-level phylogenetic inferences from multiple
genes or proteins. Demographic (and other) processes can cause single
gene trees to diverge from the species tree, so support from multiple
genes for the same tree topology is considered stronger evidence than
single gene inferences (of course, we still need to test that each gene
is telling the same story).

This is usually handled by aligning each gene separately then creating a
single "supermatrix" from the individual gene alignments, i.e. you
create a single alignment containing one row for each taxon where the
data for each row is the concatenated aligned gene sequences for the
taxon. In NEXUS files (used by the phylogenetic software PAUP\*,
MrBayes, and others) multiple genes can be explicitly represented as
different 'character partitions' or 'sets' within a data matrix that
contains one long sequence for each taxon. In this way you can create a
supermatrix but still apply different substitution models to each gene
within in it or run PAUP\*'s Partition Homogeneity Test to check for
significant difference in the rate/topology of each gene tree.

The `Bio.Nexus` module makes concatenating multiple alignments into a
supermatrix relatively straight forward.

The Solution
------------

Say we have NEXUS files for three genes,
[btCOI.nex](../examples/btCOI.nex), [btCOII.nex](../examples/btCOII.nex)
and [btITS.nex](../examples/btITS.nex), containing alignments:

    #COI
    bt1 GGGGGGGGGGGG
    bt2 GGGGGGGGGGGG
    bt3 GGGGGGGGGGGG
    #COII
    bt1 AAAAAAAAAAAA
    bt2 AAAAAAAAAAAA
    bt3 AAAAAAAAAAAA
    #ITS
    bt1 -TTTTTTT
    bt2 -TTTTTTT
    bt3 -TTTTTTT
    bt4 -TTTTTTT

We can use the `Nexus` module to make a supermatrix:

``` python
from Bio.Nexus import Nexus

# the combine function takes a list of tuples [(name, nexus instance)...],
# if we provide the file names in a list we can use a list comprehension to
# create these tuples

file_list = ["btCOI.nex", "btCOII.nex", "btITS.nex"]
nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]

combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open("btCOMBINED.nex", "w"))
```

That was easy! Let's look at our combined file

    #NEXUS
    begin data;
        dimensions ntax=4 nchar=32;
        format datatype=dna missing=? gap=-;
    matrix
    bt1 GGGGGGGGGGGG-TTTTTTTAAAAAAAAAAAA
    bt2 GGGGGGGGGGGG-TTTTTTTAAAAAAAAAAAA
    bt3 GGGGGGGGGGGG-TTTTTTTAAAAAAAAAAAA
    bt4 ????????????-TTTTTTT????????????
    ;
    end;

    begin sets;
    charset btITS.nex = 13-20;
    charset btCOI.nex = 1-12;
    charset btCOII.nex = 21-32;
    charpartition combined = btCOI.nex: 1-12, btITS.nex: 13-20, btCOII.nex: 21-32;
    end;

Ahh, it was too easy. The matrices have been combined and the character
sets and partitions set up but the ITS file had a taxon (bt4) that
wasn't in the other files. In these cases the combine function adds the
taxon with missing data (the '?'s) for the other character partitions.
Sometimes this might be the result you want but having a few taxa like
this is also a very good way to make a Partition Homogeneity Test run
for a week. Let's write a function that tests that the same taxa are
represented in a set of nexus instances and provides a useful error
message if not (i.e. what to delete from your NEXUS files if you want
them to combine nicely).

``` python
def check_taxa(matrices):
    """Verify Nexus instances have the same taxa information.

    Checks that nexus instances in a list [(name, instance)...] have
    the same taxa, provides useful error if not and returns None if
    everything matches
    """
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]
        if first_only:
            missing = ", ".join(first_only)
            msg = "%s taxa %s not in martix %s" % (nexi[0][0], missing, name)
            raise Nexus.NexusError(msg)
        elif new_only:
            missing = ", ".join(new_only)
            msg = "%s taxa %s not in all matrices" % (name, missing)
            raise Nexus.NexusError(msg)
    return None  # will only get here if it hasn't thrown an exception


def concat(file_list, same_taxa=True):
    """Combine multiple nexus data matrices in one partitioned file.

    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this
    """
    nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)
```

And now, using our new functions:


    >>> handles = [open('btCOI.nex', 'r'), open('btCOII.nex', 'r'), open('btITS.nex', 'r')]
    # If we combine them all we should get an error and the taxon/taxa that caused it
    >>> concat(handles)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "<stdin>", line 5, in concat
      File "<stdin>", line 16, in check_taxa
    Bio.Nexus.Nexus.NexusError: btITS.nex taxa bt4 not in all matrices

    # But if we use just the first two, which do have matching taxa, it should be fine
    >>> concat(handles[:2]).taxlabels
    ['bt1', 'bt2', 'bt3']

    # Ok, can we still munge them together if we want to?
    >>> concat(handle, same_taxa=False).taxlabels
    ['bt1', 'bt2', 'bt3', 'bt4']

Discussion
----------

The details of the `Nexus` class are provided in the [API
Documentation](http://www.biopython.org/DIST/docs/api/Bio.Nexus.Nexus-pysrc.html).
