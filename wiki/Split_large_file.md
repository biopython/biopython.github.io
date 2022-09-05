---
title: Split large file
permalink: wiki/Split_large_file
layout: wiki
tags:
 - Cookbook
redirect_from:
 - wiki/Split_fasta_file
---

Problem
-------

With modern sequencing technologies it has become relatively cheap and
easy to generate very large datasets. In fact, there are times when one
can have too much data in one file, online resources like
[PLAN](http://bioinfo.noble.org/plan) or
[TMHMM](http://www.cbs.dtu.dk/services/TMHMM/) limit the size of users
queries. In such cases it useful to be able to split a sequence file
into a set of smaller files, each containing a subset of original file's
sequences.

Solution
--------

There are many possible ways to solve this general problem, this recipe
uses a generator function to avoid having all the data in memory at
once.

``` python
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
```

Here is an example using this function to divide up a large FASTQ file,
`SRR014849_1.fastq` (from this [compressed
file](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR014/SRR014849/SRR014849_1.fastq.gz)):

``` python
from Bio import SeqIO

record_iter = SeqIO.parse(open("SRR014849_1.fastq"), "fastq")
for i, batch in enumerate(batch_iterator(record_iter, 10000)):
    filename = "group_%i.fastq" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fastq")
    print("Wrote %i records to %s" % (count, filename))
```

And the output:

```
Wrote 10000 records to group_1.fastq
Wrote 10000 records to group_2.fastq
Wrote 10000 records to group_3.fastq
Wrote 10000 records to group_4.fastq
Wrote 10000 records to group_5.fastq
Wrote 10000 records to group_6.fastq
...
Wrote 7251 records to group_27.fastq
```

You can modify this recipe to use any input and output formats supported
by [`Bio.SeqIO`](SeqIO "wikilink"), for example to break up a large FASTA
file into units of 1000 sequences:

``` python
from Bio import SeqIO

record_iter = SeqIO.parse(open("large.fasta"), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, 1000)):
    filename = "group_%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
```

How it works
------------

It is possible to use [`list(SeqIO.parse(...))`](SeqIO "wikilink") to read
the entire contents of a file into memory then write slices of the list
out as smaller files. For large files (like the ones this recipe is
about) that would take up a big hunk of memory, instead we can define a
generator function, `batch_iterator()`, that loads one record at a time
then appends it to a list, repeating the process until the list
containins one file's worth of sequences.

With that function defined it's a matter of giving it an iterator (in
this case a `SeqIO.parse(...)` instance and writing out the batching of
records it produces.

Note that you can get close to this functionality with the built in
Python function `itertools.islice(record_iter, batch_size)`, but this
gives an empty file at the end if the total count is an exact multiple
of the batch size.
