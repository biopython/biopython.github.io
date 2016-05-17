---
title: Methods for Degenerated Codons
permalink: wiki/Degenerated_Codons
layout: wiki
tags:
 - Cookbook
redirect_from:
 - /wiki/Degenerated
---

Problem
-------

Sometimes statistics on nucleotide sequences are limited to X-fold
degenerated codons. The following code provides some functions to solve
this problem, by generating a subsequence.

Solution
--------

``` python
from Bio.Data.CodonTable import unambiguous_dna_by_id


def altcodons(codon, table):
    """List codons that code for the same aminonacid / are also stop.

    @param codon
    @table code table id
    @return list of codons

    """
    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons

    try:
        aa = tab.forward_table[codon]
    except:
        return []

    return [k for (k, v) in tab.forward_table.iteritems()
            if v == aa and k[0] == codon[0] and k[1] == codon[1]]


def degeneration(codon, table):
    """Determine how many codons code for the same amino acid / are also stop

    @param codon the codon
    @param table code table id
    @param the number of codons also coding for the amino acid codon codes for

    """
    return len(altcodons(codon, table))


def is_x_degenerated(x, codon, table):
    """Determine if codon is x-fold degenerated.

    @param codon the codon
    @param table code table id
    @param true if x <= the degeneration of the codon

    """
    return (x <= len(altcodons(codon, table)))


def degenerated_subseq(seq, x, table):
    """Get a subsequence consisting of the x-fold degenerated codons only."""
    data = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].tostring()
        if isXdegenerated(x, codon, table):
            data += codon
    return data

```
