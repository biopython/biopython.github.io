---
title: Biopython and the Gene Ontology (GO) consortium.
permalink: wiki/Gene_Ontology
layout: wiki
---

This page describes Biopython's support for the [Gene Ontology
(GO)](http://www.geneontology.org/). Currently this is in the discussion
stages. As we implement more support, this page will evolve into
documentation for using Biopython packages to work with GO and GO
annotations.

Brief Summary of GO
===================

The GO consortium developed GO to standardize the vocabulary life
scientists use to annotate genes, and to make such annotations amenable
to computational assessment of their semantics. GO is actually three
separate but related
[ontologies](http://www.geneontology.org/GO.ontology.structure.shtml#ontology)
for three categories: biological process, cellular component, and
molecular function. The ontologies are composed of
[terms](http://www.geneontology.org/GO.ontology.structure.shtml#term).
The terms have
[relations](http://www.geneontology.org/GO.ontology.relations.shtml)
which connect them into child/ancestor relationships, where children
transitively inherit the meanings of all of their ancestors. Thus, child
nodes are more specific in their semantic meaning. Using the GO terms as
nodes, and relations as directed edges (from a child term to a parent
term), we can construct a data structure known as a [directed acyclic
graph (DAG)](http://en.wikipedia.org/wiki/Directed_acyclic_graph). We
use this DAG to traverse the ontologies and identify relationships
between GO terms.

Design of GO Support
====================

GO support should include supporting both the actual GO ontologies and
GO annotations.

Ontology Support
----------------

### Parser

We will only support the [OBO format
v.1.2](http://www.geneontology.org/GO.format.obo-1_2.shtml) in the
beginning. Parsers in Python may be around. Ed Cannon [appears to have
written an OBO to OWL
parser](http://lists.open-bio.org/pipermail/biopython-dev/2009-August/006701.html)
that is [hosted on
GitHub](http://github.com/eoc21/biopython/tree/eoc21Branch/).

### GO Directed Acyclic Graph

GO is best [represented as a directed acyclic graph
(DAG)](http://www.geneontology.org/GO.ontology.structure.shtml#go-as-a-graph).
To facilitate this data structure, we'll use
[NetworkX](https://networkx.github.io/), a popular, well-supported Python
graph library with no required dependencies other than Python. We'll use
the directed graph class
[DiGraph](https://networkx.github.io/documentation/stable/reference/classes/digraph.html) of
NetworkX to represent the ontologies.

### Sources of Code/Inspiration

#### BioPerl

[BioPerl](http://bioperl.org/) has [a package providing GO
support](http://doc.bioperl.org/bioperl-live/Bio/Ontology/toc.html). It
is not obvious yet what each component does and it seems quite heavily
engineered. Beginning GO support for Biopython will probably not reach
this level of sophistication but we can borrow from the ideas of the
BioPerl library.

#### Python Code

-   Brad Chapman's [PGML
    code](http://bioinformatics.org/cgi-bin/cvsweb.cgi/biopy-pgml/Bio/PGML/GO/),
    which apparently works with GO that's already loaded into a database

Design of GO Annotation Support
-------------------------------

GO annotations come in the form of GOA files. These are rather simple
tab-delimited flat files that contain gene identifiers, GO term IDs, and
more. The idea of providing support for GOA files is to provide a parser
that renders the file into Python data structures, and preferably as
instances of some standard Biopython.

### Parser

The GOA files are tab-delimited and quite simple, with one line per
annotation. The format is detailed at
<http://www.geneontology.org/GO.format.annotation.shtml>.

### Data Structure

We could use some sort of standard Biopython structure to store
annotations in. We could use `SeqRecord` but this would be awkward,
because it is intended to be used with a sequence. Note, however, that
`SeqRecord` does contain an `annotations` attribute, which would be
where we would like to place these GO annotations. (It's also important
to keep in mind each gene may have multiple annotations.) [Peter
Cock](User%3APeter "wikilink") points out that the QUAL parser uses a
class called `UnknownSeq` since it only contains information on the
sequence length. We should take a look at this.

One possibility is to create an `Annotation` class as a fully supported
Biopython class. More input from other Biopython developers would help
here.

As far as I can tell, this data structure really doesn't have to be
anything more than a data store, with some `__repr__()`
--[Gotgenes](User%3AGotgenes "wikilink") 05:05, 20 October 2009 (UTC)

Implementation
==============

Biopython branches for GO support development
---------------------------------------------

-   <http://github.com/gotgenes/biopython/>

