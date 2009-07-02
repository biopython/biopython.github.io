---
title: PhyloXML
permalink: wiki/PhyloXML
layout: wiki
---

This module handles the parsing, generation and manipulation of files in
the [phyloXML](http://www.phyloxml.org/) format.

This code is not yet part of Biopython, and therefore the documentation
has not been integrated into the Biopython Tutorial yet either.

Installation
------------

The source code for this module currently lives on the [phyloxml
branch](http://github.com/etal/biopython/tree/phyloxml) in GitHub. If
you're interested in testing this code before it's been merged into
Biopython, follow the instructions there to create your own fork, or
just clone the phyloxml branch onto your machine.

Requirements:

-   Biopython 1.50 or newer (older may work, but hasn't been tested)
-   Python 2.4 or newer
-   ElementTree module

The XML parser used in this module is ElementTree, new to the Python
standard library in Python 2.5. To use this module in Python 2.4, you'll
need to install a separate package that provides the ElementTree
interface. Two exist:

-   [lxml](http://codespeak.net/lxml/)
-   [elementtree](http://effbot.org/zone/element-index.htm)
    (or cElementTree)

The module attempts to import each of these compatible ElementTree
implementations until it succeeds. The given XML file handle is then
parsed incrementally to instantiate an object hierarchy containing the
relevant phylogenetic information.

Usage
-----

The most useful parts of this package are available from the top level
of the module:

``` python
from Bio import PhyloXML
```

A complete phyloXML document has a root node with the tag "phyloxml".
Directly under the root is a sequence of "phylogeny" elements
(phylogenetic trees), possibly followed by other arbitrary data not
included in the phyloXML spec. The main structural element of these
phylogenetic trees is the Clade: a tree has a clade attribute, along
with other attributes, and each clade contains a series of clades (and
other attributes), recursively.

The child nodes and attributes of each XML node are mapped onto classes
in the PhyloXML.Tree module, keeping the names the same where possible;
the XML document structure is closely mirrored in the Tree.Phyloxml
objects produced by the Bio.PhyloXML package.

### Parsing phyloXML files

This module provides two functions, read() and parse(). Both functions
accept either a file name or an open file handle, so phyloXML data can
be also loaded from compressed files, StringIO objects, and so on.

The read() function returns a single Tree.Phyloxml object representing
the entire file's data. The phylogenetic trees are in the "phylogenies"
attribute, and any other arbitrary data is stored in "other".

``` python
>>> phx = PhyloXML.read('phyloxml_examples.xml')
>>> print phx
Phyloxml
>>> len(phx.phylogenies)
13
>>> len(phx.other)
1
>>> print phx.other
[Other(tag=alignment, namespace=http://example.org/align)]
>>> print phx.other[0].children
[Other(tag=seq, namespace=http://www.phyloxml.org, value=acgtcgcggcccgtggaagtcctctcct),
Other(tag=seq, namespace=http://www.phyloxml.org, value=aggtcgcggcctgtggaagtcctctcct),
Other(tag=seq, namespace=http://www.phyloxml.org, value=taaatcgc--cccgtgg-agtccc-cct)]
```

If you aren't interested in the "other" data, you can use parse() to
iteratively construct just the phylogenetic trees contained in the file.
If there's only one tree, then the next() method on the resulting
generator will return it.

``` python
>>> for tree in PhyloXML.parse('phyloxml_examples.xml'):
...     print tree.name
```

    example from Prof. Joe Felsenstein's book "Inferring Phylogenies"
    example from Prof. Joe Felsenstein's book "Inferring Phylogenies"
    same example, with support of type "bootstrap"
    same example, with species and sequence
    same example, with gene duplication information and sequence relationships
    similar example, with more detailed sequence data
    network, node B is connected to TWO nodes: AB and C
    ...

``` python
>>> tree = PhyloXML.parse('phyloxml_examples.xml').next()
>>> tree.name
'example from Prof. Joe Felsenstein\'s book "Inferring Phylogenies"'
```

### Writing phyloXML files

The PhyloXML.Writer module exports one function to the top level of the
package: write(). It accepts a Phyloxml object (the result of read() or
to\_phyloxml()) and either a file name or a handle to an open file-like
object. Optionally, an encoding other than UTF-8 can be specified.

``` python
>>> phx = PhyloXML.read('phyloxml_examples.xml')
>>> print phx.other
[Other(tag=alignment, namespace=http://example.org/align)]
>>> phx.other = []
>>> PhyloXML.write(phx, 'ex_no_other.xml')
>>> phx_no = PhyloXML.read('ex_no_other.xml')
>>> phx_no.other
[]
```

### Using PhyloXML objects

Standard Python syntactic sugar is supported wherever it's reasonable.

-   \_\_str\_\_ makes a string of the object's class name and an
    identifier, suitable for labeling a node in generated graph
-   \_\_repr\_\_ makes a string in the usual Biopython format,
    resembling the object constructor call
-   \_\_iter\_\_ is supported by Phyloxml and Clade objects, iterating
    over the contained phylogenies and sub-clades, respectively
-   \_\_len\_\_ is supported by the same objects that support iteration,
    with expected results

Clade objects also support extended indexing:

``` python
tree = PhyloXML.parse('example.xml').next()
assert tree.clade[0] == tree.clade.clades[0]
assert tree.clade[0,1] == tree.clade.clades[0].clades[1]
```

Since valid Phylogeny objects always have a single clade attribute, this
style of indexing is a handy way to reach specific nodes buried deep in
the tree if you happen to know exactly where they are. In the future,
slicing may be supported as well to grab a range a sub-clades at once.

A couple of methods allow converting a selection to a new PhyloXML
object: Phylogeny.to\_phyloxml() and Clade.to\_phylogeny(). A few use
cases:

-   Parse a phyloXML containing multiple phylogenetic trees. Check each
    tree sequentially, and upon finding a tree with the desired
    characteristic, isolate it as a new PhyloXML object.

``` python
for tree in PhyloXML.parse('example.xml'):
    if tree.name == 'monitor lizards':
        return tree.to_phyloxml()
```

-   Extract a specific sub-clade and make it a separate phylogeny (and
    probably a new phyloXML file).

``` python
tree = PhyloXML.parse('example.xml').next()
best = None
for clade in tree.clade:
    if (clade.confidences[0].type == 'bootstrap'
            and (best is None
                or clade.confidences[0].value > best.confidences[0].value)):
        best = clade
phyloxml = best.to_phylogeny(rooted=True).to_phyloxml()
PhyloXML.write(phyloxml, 'example_best.xml')
```

### Integrating with the rest of Biopython

Since the phyloXML specification is very detailed, no existing Biopython
classes are currently used directly by the phyloXML parser. Instead,
methods are available for converting between phyloXML and standard
Biopython types.

The Tree.Sequence class contains methods for converting to and from
Biopython [SeqRecord](SeqRecord "wikilink") objects. This includes the
molecular sequence (mol\_seq) as a [Seq](Seq "wikilink") object, and the
protein domain architecture as list of
[SeqFeature](SeqFeature "wikilink") objects.

*At some point this module should be merged into the Biopython trunk,
and it would be nice to have a common interface with Bio.Nexus and
Newick. Should these three modules be reorganized to extract a common
Bio.TreeIO interface? Let's discuss it at some point.*

### Utilities

Some additional tools are located in Bio.PhyloXML.Utils.

-   dump\_tags -- Print the XML tags as they are encountered in
    the file. For testing and debugging, mainly.

<!-- -->

    >>> PhyloXML.dump_tags('phyloxml_examples.xml')
    {http://www.phyloxml.org}phyloxml
    {http://www.phyloxml.org}phylogeny
    {http://www.phyloxml.org}name
    {http://www.phyloxml.org}description
    {http://www.phyloxml.org}clade
    ...

-   pretty\_print -- Produces a plain-text representation of the
    entire tree. Uses str() to display nodes by default; for the
    longer repr() representation, add the argument show\_all=True.

<!-- -->

    >>> PhyloXML.pretty_print('phyloxml_examples.xml')
    Phyloxml
        Phylogeny example from Prof. Joe Felsenstein's book "Inferring Phylogenies"
            Clade
                Clade
                    Clade A
                    Clade B
                Clade C
    ...

    >>> PhyloXML.pretty_print('phyloxml_examples.xml', show_all=True)
    Phyloxml()
        Phylogeny(description=phyloXML allows to use either a "branch_length"
    attribute or element to indicate branch lengths., name=example from Prof. Joe
    Felsenstein's book "Inferring Phylogenies")
            Clade()
                Clade(branch_length=0.06)
                    Clade(branch_length=0.102, name=A)
                    Clade(branch_length=0.23, name=B)
                Clade(branch_length=0.4, name=C)
    ...

### Performance

This parser is meant to be able to handle large files, meaning several
thousand external nodes. (Benchmarks of relevant XML parsers for Python
are [here](http://effbot.org/zone/celementtree.htm#benchmarks).) It has
been tested with files of this size; for example, the complete NCBI
taxonomy parses in about 100 seconds and consumes about 1.3 GB of
memory. Provided enough memory is available on the system, the writer
can also rebuild phyloXML files of this size.

The read() and parse() functions process a complete file in about the
same amount of CPU time. Most of the underlying code is the same, and
the majority of the time is spent building Clade objects. For small
files (smaller than ncbi\_taxonomy\_mollusca.xml), the write() function
serializes the complete object back to an equivalent file in around
twice the CPU time required by the corresponding read() call; for large
files write() is more efficient, finishing up to 3.6 times as fast as
read() for the largest file tested.

Here are some times on a 2.00GHz Intel Xeon E5405 processor (only 1 CPU
core used) with 7.7GB memory:

| File                         | Size (uncompressed) | Read (s) | Parse (s) | Write (s) |
|------------------------------|---------------------|----------|-----------|-----------|
| phyloxml\_examples.xml       | 15 KB               | 0.0056   | 0.0053    | 0.0103    |
| apaf.xml                     | 38 KB               | 0.0112   | 0.0112    | 0.0252    |
| bcl\_2.xml                   | 105 KB              | 0.0293   | 0.0294    | 0.0455    |
| ncbi\_taxonomy\_mollusca.xml | 1.5 MB              | 0.688    | 0.666     | 0.683     |
| ncbi\_taxonomy\_metazoa.xml  | 33 MB               | 16.383   | 16.986    | 9.026     |
| ncbi\_taxonomy.xml           | 31 MB (unindented)  | 99.715   | 100.588   | 27.166    |

Summer of Code project
----------------------

This module is being developed by [Eric
Talevich](User%3AEricTalevich "wikilink") as a project for Google Summer
of Code 2009, with NESCent as the mentoring organization and Brad
Chapman as the primary mentor.

Main SoC project page: [PhyloSoC:Biopython support for parsing and
writing
phyloXML](https://www.nescent.org/wg_phyloinformatics/PhyloSoC:Biopython_support_for_parsing_and_writing_phyloXML)

Other software
--------------

[Christian Zmasek](http://monochrome-effect.net/), one of the authors of
the phyloXML specification, has released some software that uses this
format:

-   [Forester](http://www.phylosoft.org/forester/) -- a collection of
    Java and Ruby libraries for working with phylogenetic data
-   [Archaopteryx](http://www.phylosoft.org/archaeopteryx/) -- Java
    application for the visualization of annotated phylogenetic trees
    (also available in applet form)

Another list is maintained at
[phylosoft.org](http://www.phylosoft.org/).
