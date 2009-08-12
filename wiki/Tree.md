---
title: Tree
permalink: wiki/Tree
layout: wiki
---

This module provides the classes for representing phylogenetic trees.

This code is not yet part of Biopython, and therefore the documentation
has not been integrated into the Biopython Tutorial yet either.

Availability
------------

This module is being developed alongside phyloXML support as a Google
Summer of Code 2009 project. To use Tree, see the
[PhyloXML](PhyloXML "wikilink") page for instructions on installing from
a [GitHub](GitUsage "wikilink") branch.

Basic objects (BaseTree)
------------------------

See the [PhyloXML](PhyloXML "wikilink") page for examples of using Tree
objects.

Format-specific extensions
--------------------------

To support additional information stored in specific file formats,
sub-modules within Tree offer additional classes that inherit from
BaseTree classes.

Each sub-class of BaseTree.Tree or Node has a class method to promote an
object from the basic type to the format-specific one. These sub-class
objects can generally be treated as instances of the basic type without
any explicit conversion.

### PhyloXML

Support for the phyloXML format. See the [PhyloXML](PhyloXML "wikilink")
page for details.

### Newick

Porting Newick objects to Bio.Tree is in progress, but not usable yet.

Utilities
---------

Some additional tools are located in the Utils module under Bio.Tree.
These functions are loaded to the top level of the module on import for
easy access:

``` python
from Bio import Tree
```

Where a third-party package is required, that package is imported when
the function itself is called, so these dependencies are not necessary
to install or use the rest of the Tree module.

### pretty\_print()

Produces a plain-text representation of the entire tree. Uses str() to
display nodes by default; for the longer repr() representation, add the
argument show\_all=True.

Strings are automatically truncated to ensure reasonable display.

    >>> phx = PhyloXMLIO.read('phyloxml_examples.xml')
    >>> Tree.pretty_print(phx)
    Phyloxml
        Phylogeny example from Prof. Joe Felsenstein's book "Inferring Phylogenies"
            Clade
                Clade
                    Clade A
                    Clade B
                Clade C
    ...

    >>> Tree.pretty_print(phx, show_all=True)
    Phyloxml()
        Phylogeny(description='phyloXML allows to use either a "branch_length"
    attribute or element to indicate branch lengths.', name='example from Prof. Joe
    Felsenstein's book "Inferring Phylogenies"')
            Clade()
                Clade(branch_length=0.06)
                    Clade(branch_length=0.102, name='A')
                    Clade(branch_length=0.23, name='B')
                Clade(branch_length=0.4, name='C')
    ...

### Graph export

Although any phylogenetic tree can reasonably be represented by a
directed acyclic graph, the Tree module does not attempt to provide a
generally usable graph library. Instead, it provides two functions for
exporting Tree objects to the excellent
[NetworkX](http://networkx.lanl.gov/) library's graph objects and using
that library, along with
[matplotlib](http://matplotlib.sourceforge.net/) and/or
[PyGraphviz](http://networkx.lanl.gov/pygraphviz/), to display trees.

**to\_networkx** returns the given tree as a networkx LabeledDiGraph or
LabeledGraph object. You'll probably need to import networkx directly
for subsequent operations on the graph object. To display a (somewhat
haphazard-looking) dendrogram on the screen, for interactive work, use
matplotlib or pylab along with one of networkx's drawing functions.

``` python
import networkx, pylab
tree = TreeIO.read('example.xml', 'phyloxml')
net = Tree.to_networkx(tree)
networkx.draw(net)
pylab.show()
```

**draw\_graphviz** mimics the networkx function of the same name, with
some tweaks to improve the display of the graph. If a file name is
given, the graph is drawn directly to that file, and options such as
image format (default PDF) may be used.

``` python
tree = TreeIO.read('example.xml', 'phyloxml')
# Draw it a few different ways
Tree.draw_graphviz(tree, 'example.pdf')
Tree.draw_graphviz(tree, 'example.png', format='png')
Tree.draw_graphviz(tree)
```

Using with BioSQL
-----------------

BioSQL offers an extension for storing phylogenetic trees called
[PhyloDB](http://biosql.org/wiki/Extensions). The BaseTree classes were
written with PhyloDB compatibility in mind, but no BioSQL adapter for
Bio.Tree has been written yet.
