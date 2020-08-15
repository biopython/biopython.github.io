---
title: Retrieve nonmatching blast queries
permalink: wiki/Retrieve_nonmatching_blast_queries
layout: wiki
tags:
 - Cookbook
redirect_from:
 - wiki/Retrieve_nonmathching_blast_queries
---

Problem
-------

The XML output of NCBI's stand alone BLAST programs does not include
information on query sequences that have 'no hits' in the target
database. [Sometimes](https://redmine.open-bio.org/issues/2821)
you want to know which sequences don't have match a database and further
analyse/anotate them accordingly. There are a number of different ways
to do this, one is to use [SeqIO's method `.index()`](SeqIO "wikilink") to
turn the query file into a dictionary then parse the results file to get
the sequences that *did* match the dictionary. You can then use Python's
`set()` arithmetic to make a list of sequences that are in the query file
and not the results which can be used as keys to retrieve the complete
[`SeqRecord`](SeqRecord "wikilink") for each of the "no hit" queries. Got
it? Well, perhaps it's easier to just do it:

Solution
--------

Let's presume you set up a BLAST run with the sequences in a file called
`queries.fasta` searched against a database, with the results saved to
`BLAST_RESULTS.xml`

``` python
from Bio import SeqIO
from Bio.Blast import NCBIXML

# Build an index, but we don't need to parse the record
q_dict = SeqIO.index("queries.fasta", "fasta")

hits = []
for record in NCBIXML.parse(open("BLAST_RESULTS.xml")):
    # As of BLAST 2.2.19 the xml output for multiple-query blast searches
    # skips queries with no hits so we could just append the ID of every blast
    # record to our 'hit list'. It's possible that the NCBI will change this
    # behaviour in the future so let's make the appending conditional on there
    # being some hits (ie, a non-empty alignments list) recorded in the blast
    # record

    if record.alignments:
        # The blast record's 'query' contains the sequences description as a
        # string. We used the ID as the key in our dictionary so we'll need to
        # split the string and get the first field to remove the right entries
        hits.append(record.query.split()[0])

misses = set(q_dict.keys()) - set(hits)
orphan_records = [q_dict[name] for name in misses]
```

We can do a little sanity check to make sure everything worked OK:

``` pycon
>>> print("found %i records in query, %i have hits, making %i misses"
...       % (len(q_dict), len(hits), len(misses)))
...
found 11955 records in query, 2802 have hits, making 9153 misses
```

Good, now you have all the 'not hits' sequence in a list
(`orphan_records`) of `SeqRecord` objects you can annotate/analyse as you
please or use [`SeqIO.write()`](SeqIO "wikilink") to make a new file of
just these sequences that can be put through another program.

Discussion
----------

As implemented above most of the time in each run is spend populating
the list of hits from the BLAST parser, would checking each record from
the results file against the dictionary one at a time be a less memory
intensive way to go in case of very large files?
