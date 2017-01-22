---
title: Retrieve and annotate Entrez Gene IDS with the Entrez module.
permalink: wiki/Annotate_Entrez_Gene_IDs
layout: wiki
tags:
 - Cookbook
---

Problem
-------

If you deal with a large quantity of gene IDs (such as the ones produced
by microarray analysis), annotating them is important if you want to
determine their potential biological meaning. However, a lot of
annotation systems are only web-based, or do not work with Python.

Solution
--------

Thanks to the Entrez module it is possible to annotate batches of Entrez
Gene IDs without trouble, using the "retrieve\_ids" function provided
below.

This example assumes you have a list of Entrez Gene IDs. *Note*: they
should be stored as strings, rather than integers, even if they are
numbers.

``` python

import sys

from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "your email here"

def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "An error occurred while retrieving the annotations."
        print "The error returned was %s" % e
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)

    print "Retrieved %d annotations for %d genes" % (len(annotations),
            len(id_list))

    return annotations
```

As the function's docstring says, the function returns a list of
dictionaries, one for each gene, that you can use in any way you want.
The following example prints out ID, Gene Symbol and Gene Name for a
retrieved annotation:

``` python

def print_data(annotation):
    for gene_data in annotation:
        gene_id = gene_data["Id"]
        gene_symbol = gene_data["NomenclatureSymbol"]
        gene_name = gene_data["Description"]
        print "ID: %s - Gene Symbol: %s - Gene Name: %s" % (gene_id, gene_symbol, gene_name)
```

More
----

Tao Liu expanded on this code with a [full example of converting between
Entrez IDs and RefSeq
IDs](https://github.com/taoliu/taolib/blob/master/Scripts/convert_gene_ids.py).

------------------------------------------------------------------------
