---
title: From gene sequence to predicted protein with the GFF module.
permalink: wiki/Gene_predictions_to_protein_sequences
layout: wiki
tags:
 - Cookbook
---

Contributed by Brad Chapman

Note: This requires the proposed Biopython [GFF
module](http://github.com/chapmanb/bcbb/tree/master/gff/), which has not
yet been integrated in Biopython. You will need to install it first to
run the examples

Problem
-------

Starting with a
[GlimmerHMM](http://www.cbcb.umd.edu/software/GlimmerHMM/) output file
in [GFF3](http://www.sequenceontology.org/gff3.shtml) format, produce a
FASTA file of predicted protein sequences.

Solution
--------

Setting this up, we import the required modules and parse our input
FASTA file into a standard python dictionary, using
[SeqIO](SeqIO "wikilink"). SeqIO is also used for writing the output
file. We pass the output function a python generator, which SeqIO will
convert into FASTA files one at a time. This allows the script to be
generalized to arbitrarily large data sets without running into memory
issues.

``` python
from __future__ import with_statement
import sys
import os
import operator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from BCBio import GFF

def main(glimmer_file, ref_file):
    with open(ref_file) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    base, ext = os.path.splitext(glimmer_file)
    out_file = "%s-proteins.fa" % base
    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(glimmer_file, ref_recs), out_handle, "fasta")
```

Next, we provide a function to parse the GlimmerHMM output. Thanks to
the [GFF](GFF_Parsing "wikilink") library, this is very easy. We iterate
over a group of lines at a time to handle very large output files, and
provide as input the reference records we parsed earlier with SeqIO:

``` python
def glimmer_predictions(in_handle, ref_recs):
    """Parse Glimmer output, generating SeqRecord and SeqFeatures for predictions
    """
    for rec in GFF.parse(in_handle, target_lines=1000, base_dict=ref_recs):
        yield rec
```

This is the meat of the program. We loop over each of the records in the
GlimmerHMM output, then over each of the gene features in that record.
Each of these SeqFeatures contains a gene prediction with exons present
as sub\_features of the top level feature. Using the location of these
exons, we extract the nucleotide sequence for the current gene and
append it to a list. We then add together the list to get the full
sequence. If the gene is on the reverse strand, the sequence is reverse
complemented. Finally, we translate the sequence to protein. A SeqRecord
with the sequence and prediction identifier is yielded, providing the
generator we used earlier to write the output file:

``` python
def protein_recs(glimmer_file, ref_recs):
    """Generate protein records from GlimmerHMM gene predictions.
    """
    with open(glimmer_file) as in_handle:
        for rec in glimmer_predictions(in_handle, ref_recs):
            for feature in rec.features:
                seq_exons = []
                for cds in feature.sub_features:
                    seq_exons.append(rec.seq[
                        cds.location.nofuzzy_start:
                        cds.location.nofuzzy_end])
                gene_seq = reduce(operator.add, seq_exons)
                if feature.strand == -1:
                    gene_seq = gene_seq.reverse_complement()
                protein_seq = gene_seq.translate()
                yield SeqRecord(protein_seq, feature.qualifiers["ID"][0], "", "")
```

Discussion
----------

Compare the full scripts for:

-   [GFF3 to
    proteins](http://github.com/chapmanb/bcbb/blob/master/biopython/glimmergff_to_proteins.py)
-   [Custom output to
    proteins](http://github.com/chapmanb/bcbb/blob/master/biopython/glimmer_to_proteins.py)

This highlights how using standard libraries with standardized formats
can save coding versus rolling your own. The GFF3 version is more
general since it allows multiple contigs to be processed concurrently,
and requires less custom code to solve the problem. This allows you to
focus on the biological problem, and avoid bugs in the parsing code.

Example files to run the scripts are available thanks to Michael Kuhn:

-   [GFF3 GlimmerHMM
    file](http://github.com/chapmanb/bcbb/blob/master/biopython/glimmer.gff)
-   [Reference
    file](http://github.com/chapmanb/bcbb/blob/master/biopython/ref.fa)

See also:

-   [The FriendFeed post that inspired this
    entry](http://friendfeed.com/the-life-scientists/553c495e/does-anyone-know-how-to-take-gff3-file-gene)
-   [Iddo's post on parsing standard Glimmer
    output](http://bytesizebio.net/index.php/2009/10/29/short-bioinformatics-hacks-glimmer-splitter/)
