---
title: Phylo
permalink: wiki/Phylo
layout: wiki
tags:
 - Wiki Documentation
---

This module provides classes, functions and I/O support for working with
phylogenetic trees.

This code is not yet included with Biopython, and therefore the
documentation has not been integrated into the Biopython Tutorial yet
either.

Availability
------------

The source code for this module is being developed on GitHub, but has
not yet been included with a stable Biopython release. If you're
interested in testing this code before the official release, see
[SourceCode](SourceCode "wikilink") for instructions on getting a copy
of the development branch.

Requirements:

-   Python 2.4 or newer
-   ElementTree module

To draw trees (optional), you'll also need these packages:

-   [NetworkX](http://networkx.lanl.gov/index.html) 1.0rc1 (or 0.36 for
    snapshot at the end of GSoC 2009)
-   [PyGraphviz](http://networkx.lanl.gov/pygraphviz/) 0.99.1 or
    [pydot](http://dkbza.org/pydot.html)
-   [matplotlib](http://matplotlib.sourceforge.net/)

The I/O and tree-manipulation functionality will work without them;
they're imported on demand when the functions to\_networkx() and
draw\_graphviz() are called.

The XML parser used in the PhyloXMLIO sub-module is ElementTree, added
to the Python standard library in Python 2.5. To use this module in
Python 2.4, you'll need to install a separate package that provides the
ElementTree interface. Two exist:

-   [lxml](http://codespeak.net/lxml/)
-   [elementtree](http://effbot.org/zone/element-index.htm)
    (or cElementTree)

PhyloXMLIO attempts to import each of these compatible ElementTree
implementations until it succeeds. The given XML file handle is then
parsed incrementally to instantiate an object hierarchy containing the
relevant phylogenetic information.

I/O functions
-------------

Wrappers for supported file formats are available from the top level of
the module:

``` python
from Bio import Phylo
```

Like SeqIO and AlignIO, this module provides four I/O functions:
parse(), read(), write() and convert(). Each function accepts either a
file name or an open file handle, so data can be also loaded from
compressed files, StringIO objects, and so on. If the file name is
passed as a string, the file is automatically closed when the function
finishes; otherwise, you're responsible for closing the handle yourself.

The second argument to each function is the target format. Currently,
the following formats are supported:

-   phyloxml
-   newick
-   nexus

See the [PhyloXML](PhyloXML "wikilink") page for more examples of using
tree objects.

### parse()

Incrementally parse each tree in the given file or handle, returning an
iterator of Tree objects (i.e. some subclass of the Bio.Phylo.BaseTree
Tree class, depending on the file format).

``` python
>>> trees = Phylo.parse('phyloxml_examples.xml', 'phyloxml')
>>> for tree in trees:
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

If there's only one tree, then the next() method on the resulting
generator will return it.

``` python
>>> tree = Phylo.parse('phyloxml_examples.xml', 'phyloxml').next()
>>> tree.name
'example from Prof. Joe Felsenstein\'s book "Inferring Phylogenies"'
```

Note that this doesn't immediately reveal whether there are any
remaining trees -- if you want to verify that, use read() instead.

### read()

Parse and return exactly one tree from the given file or handle. If the
file contains zero or multiple trees, a ValueError is raised. This is
useful if you know a file contains just one tree, to load that tree
object directly rather than through parse() and next(), and as a safety
check to ensure the input file does in fact contain exactly one
phylogenetic tree at the top level.

``` python
tree = Phylo.read('example.xml', 'phyloxml')
print tree
```

### write()

Write a sequence of Tree objects to the given file or handle. Passing a
single Tree object instead of a list or iterable will also work. (See,
Phylo is friendly.)

``` python
tree1 = Phylo.read('example1.xml', 'phyloxml')
tree2 = Phylo.read('example2.xml', 'phyloxml')
Phylo.write([tree1, tree2], 'example-both.xml', 'phyloxml')
```

### convert()

Given two files (or handles) and two formats, both supported by
Bio.Phylo, convert the first file from the first format to the second
format, writing the output to the second file.

``` python
Phylo.convert('example.nhx', 'newick', 'example2.nex', 'nexus')
```

### Sub-modules

Within the Phylo module are parsers and writers for specific file
formats, conforming to the basic top-level API and sometimes adding
additional features.

**PhyloXMLIO:** Support for the [phyloXML](http://www.phyloxml.org/)
format. See the [PhyloXML](PhyloXML "wikilink") page for details.

**NewickIO:** A port of the parser in Bio.Nexus.Trees to support the
Newick (a.k.a. New Hampshire) format through the Phylo API.

**NexusIO:** Wrappers around Bio.Nexus to support the Nexus tree format.

The Nexus format actually contains several sub-formats for different
kinds of data; to represent trees, Nexus provides a block containing
some metadata and one or more Newick trees. (Another kind of Nexus block
can represent alignments; this is handled in
[AlignIO](AlignIO "wikilink").) So to parse a complete Nexus file with
all block types handled, use Bio.Nexus directly, and to extract just the
trees, use Bio.Phylo. Integration between Bio.Nexus and Bio.Phylo will
be improved in the future.

Tree and Subtree classes
------------------------

The basic objects are defined in Bio.Phylo.BaseTree.

### Format-specific extensions

To support additional information stored in specific file formats,
sub-modules within Tree offer additional classes that inherit from
BaseTree classes.

Each sub-class of BaseTree.Tree or Node has a class method to promote an
object from the basic type to the format-specific one. These sub-class
objects can generally be treated as instances of the basic type without
any explicit conversion.

**PhyloXML**: Support for the phyloXML format. See the
[PhyloXML](PhyloXML "wikilink") page for details.

**Newick**: The Newick module provides minor enhancements to the
BaseTree classes, plus several shims for compatibility with the existing
Bio.Nexus module. The API for this module is under development and
should not be relied on, other than the functionality already provided
by BaseTree.

Utilities
---------

Some additional tools are located in the Utils module under Bio.Phylo.
These functions are also loaded to the top level of the Phylo module on
import for easy access.

Where a third-party package is required, that package is imported when
the function itself is called, so these dependencies are not necessary
to install or use the rest of the Tree module.

### Exporting to other object representations

Although any phylogenetic tree can reasonably be represented by a
directed acyclic graph, the Phylo module does not attempt to provide a
generally usable graph library -- only the minimum functionality to
represent phylogenetic trees. Instead, it provides functions for
exporting tree objects to the standard graph representations, adjacency
list (dict) and adjacency matrix, using third-party libraries.

**to\_networkx** returns the given tree as a
[NetworkX](http://networkx.lanl.gov/) LabeledDiGraph or LabeledGraph
object (depending on whether the tree is rooted). You'll probably need
to import networkx directly for subsequent operations on the graph
object. From this point you can also try using one of networkx's drawing
functions to display the tree, and for simple, fully labeled trees it
may even work -- but you'll have better results with Phylo's own
draw\_graphviz function, discussed below.

``` python
import networkx, pylab
tree = Phylo.read('example.xml', 'phyloxml')
net = Phylo.to_networkx(tree)
networkx.draw(net)
pylab.show()
```

**to\_adjacency\_matrix** produces an adjacency matrix as an instance of
a NumPy 2-dimensional array, where cell values are branch lengths and
rows and columns are vertices in the graph (i.e. nodes in the tree, the
root of each clade). The returned tuple includes a list of all clade
objects in the original tree, used for determining the indexes of cells
in the matrix corresponding to clades or branches in the tree.

### Displaying trees

**pretty\_print** produces a plain-text representation of the entire
tree. Uses str() to display nodes by default; for the longer repr()
representation, add the argument show\_all=True.

Strings are automatically truncated to ensure reasonable display.

    >>> phx = Phylo.parse('phyloxml_examples.xml', 'phyloxml').next()
    >>> Phylo.pretty_print(phx)
    Phylogeny example from Prof. Joe Felsenstein's book "Inferring Phylogenies"
        Clade
            Clade
                Clade A
                Clade B
            Clade C
    ...

    >>> Phylo.pretty_print(phx, show_all=True)
    Phylogeny(description='phyloXML allows to use either a "branch_length"
    attribute or element to indicate branch lengths.', name='example from
    Prof. Joe Felsenstein's book "Inferring Phylogenies"')
        Clade()
            Clade(branch_length=0.06)
                Clade(branch_length=0.102, name='A')
                Clade(branch_length=0.23, name='B')
            Clade(branch_length=0.4, name='C')
    ...

**draw\_graphviz** mimics the networkx function of the same name, with
some tweaks to improve the display of the graph. If a file name is
given, the graph is drawn directly to that file, and options such as
image format (default PDF) may be used.

<img src="Phylo-apaf.png" title="Phylogram with colored nodes" alt="Phylogram with colored nodes" width="256" />

Prerequisites: In addition to networkx, you'll need a local installation
of Graphviz, [matplotlib](http://matplotlib.sourceforge.net/) and either
[PyGraphviz](http://networkx.lanl.gov/pygraphviz/) or
[pydot](http://dkbza.org/pydot.html).

Drawing a basic dendrogram is simple:

``` python
import pylab
tree = Phylo.read('apaf.xml', 'phyloxml')
Phylo.draw_graphviz(tree)
pylab.show()
```

<img src="Phylo-apaf-node0.png" title="fig:Phylogram with plain text nodes" alt="Phylogram with plain text nodes" width="256" />
Here's the same tree without the circles at each labelled node:

``` python
Phylo.draw_graphviz(tree, node_size=0)
```

See the function's docstring for more explanation.

*TODO: Set up a cookbook page to demonstrate the more exotic options.*

**draw\_ascii** prints an ascii-art rooted phylogram to standard output,
or another file handle if specified. Only terminal node labels are
shown; these are the result of *str(clade)* (usually clade names). The
width of the text field used for drawing is 80 characters by default,
adjustable with the *column\_width* keyword argument, and the height in
character rows is twice the number of terminals in the tree.

A simple tree with defined branch lengths looks like this:

    >>> tree = Phylo.parse('phyloxml_examples.xml', 'phyloxml').next()
    >>> Phylo.draw_ascii(tree)
              _____________ A
      _______|
    _|       |_______________________________ B
     |
     |_______________________________________________________ C

The same topology without branch lengths is drawn with equal-length
branches:

                                  ___________________________ A
      ___________________________|
    _|                           |___________________________ B
     |
     |___________________________ C

A larger tree (apaf.xml, 31 leaf nodes) drawn with the default column
width demonstrates how relatively short branches are handled:

    >>> apaf = Phylo.read('apaf.xml', 'phyloxml')
    >>> Phylo.draw_ascii(apaf)
                                       _ 22_MOUSE
                                      |
                                     _| Apaf-1_HUMAN
                                    | |
                                   ,| | 12_CANFA
                                   ||
                                  _||___ 11_CHICK
                                 | |
                                 | |___________ 16_XENLA
                          _______|
                         |       |       , 14_FUGRU
                         |       |     __|
                         |       |____|  |__ 15_TETNG
                    _____|            |
                   |     |            |____ 17_BRARE
                   |     |
                   |     |    ______ 1_BRAFL
                   |     | __|
             ______|     ||  |_________ 18_NEMVE
            |      |      |
            |      |      |____________ 23_STRPU
            |      |
           _|      |          _________ 26_STRPU
          | |      |_________|
          | |                |________ 25_STRPU
          | |
          | |                                    ___ CED4_CAEEL
          | |___________________________________|
      ____|                                     |_ 31_CAEBR
     |    |
     |    |                                ___ 28_DROPS
     |    |          _____________________|
     |    |   ______|                     |____ Dark_DROME
     |    |  |      |
     |    |__|      |_______________________ 29_AEDAE
     |       |
     |       |__________________________ 30_TRICA
     |
     |                                                           _ 34_BRAFL
     |                                 _________________________|
    _|                           _____|                         |_ 35_BRAFL
     |                          |     |
     |                        __|     |_______________ 8_BRAFL
     |                       |  |
     |                       |  |        ___________________ 20_NEMVE
     |         ______________|  |_______|
     |        |              |          |__________________________ 21_NEMVE
     |        |              |
     |     ___|              |______________________________ 9_BRAFL
     |    |   |
     |    |   |                _____________ 3_BRAFL
     |    |   |          _____|
     |    |   |_________|     |_________________ 2_BRAFL
     |____|             |
          |             |_______________ 19_NEMVE
          |
          |                                     _____ 37_BRAFL
          |            ________________________|
          |___________|                        |____ 36_BRAFL
                      |
                      |______________________ 33_BRAFL
