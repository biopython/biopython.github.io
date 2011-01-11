---
title: PhyloXML
permalink: wiki/PhyloXML
layout: wiki
---

Within Bio.[Phylo](Phylo "wikilink"), Biopython's module for working
with phylogenetic trees, the PhyloXML and PhyloXMLIO sub-modules handle
the parsing, generation and manipulation of files in the
[phyloXML](http://www.phyloxml.org/) format.

About the format
----------------

A complete phyloXML document has a root node with the tag "phyloxml".
Directly under the root is a sequence of "phylogeny" elements
(phylogenetic trees), possibly followed by other arbitrary data not
included in the phyloXML spec. The main structural element of these
phylogenetic trees is the Clade: a tree has a clade attribute, along
with other attributes, and each clade contains a series of clades (and
other attributes), recursively.

The child nodes and attributes of each XML node are mapped onto classes
in the PhyloXML module, keeping the names the same where possible; the
XML document structure is closely mirrored in the Phyloxml objects
produced by Bio.Phylo.PhyloXMLIO.read(), and the Phylogeny objects
produced by Bio.Phylo.read() and parse().

For example, this XML (from Tests/PhyloXML/example.xml):

    <?xml version="1.0" encoding="UTF-8"?>
    <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org">
       <phylogeny rooted="true">
          <name>An example</name>
          <clade>
             <clade branch_length="0.06">
                <clade branch_length="0.102">
                   <name>A</name>
                </clade>
                <clade branch_length="0.23">
                   <name>B</name>
                </clade>
             </clade>
             <clade branch_length="0.4">
                <name>C</name>
             </clade>
          </clade>
       </phylogeny>
    </phyloxml>

produces an object hierarchy like this:

``` python
>>> from Bio import Phylo
>>> tree = Phylo.read('example.xml','phyloxml')
>>> print tree
```

    Phylogeny(rooted='True', description='phyloXML allows to use either a "branch_length" attribute...',
              name='example from Prof. Joe Felsenstein's book "Inferring Phyl...')
        Clade()
            Clade(branch_length='0.06')
                Clade(branch_length='0.102', name='A')
                Clade(branch_length='0.23', name='B')
            Clade(branch_length='0.4', name='C')

which represents a phylogeny like this:

``` python
>>> Phylo.draw_ascii(tree)
```

                 __________________ A
      __________|
    _|          |___________________________________________ B
     |
     |___________________________________________________________________________ C

The tree objects are derived from base classes in
[Bio.Phylo](Phylo "wikilink"); see that page for more about this object
representation.

I/O functions
-------------

To start working with phyloXML files, use the [Phylo](Phylo "wikilink")
package with 'phyloxml' as the format argument:

``` python
>>> from Bio import Phylo
>>> tree = Phylo.read('some-trees.xml', 'phyloxml')
# ValueError: There are multiple trees in this file; use parse() instead.
>>> trees = Phylo.parse('some-trees.xml', 'phyloxml')
>>> Phylo.write(trees.next(), 'first-tree.xml', 'phyloxml')
1
>>> Phylo.write(trees, 'rest-trees.xml', 'phyloxml')
12
```

These functions work with Phylogeny objects (derived from BaseTree.Tree)
from the Bio.Phylo.PhyloXML module. This standard API is enough for most
use cases.

### PhyloXMLIO

Within Bio.Phylo, the I/O functions for the phyloXML format are
implemented in the PhyloXMLIO sub-module. For access to some additional
functionality beyond the basic Phylo I/O API, or to skip specifying the
'phyloxml' format argument each time, this can be imported directly:

``` python
from Bio.Phylo import PhyloXMLIO
```

The read() function returns a single Bio.Phylo.PhyloXML.Phyloxml object
representing the entire file's data. The phylogenetic trees are in the
"phylogenies" attribute, and any other arbitrary data is stored in
"other".

``` python
>>> phx = PhyloXMLIO.read('phyloxml_examples.xml')
>>> print phx
Phyloxml
>>> len(phx.phylogenies)
13
>>> len(phx.other)
1
>>> print phx.other
[Other(tag='alignment', namespace='http://example.org/align')]
>>> print phx.other[0].children
[Other(tag='seq', namespace='http://www.phyloxml.org', value='acgtcgcggcccgtggaagtcctctcct'),
Other(tag='seq', namespace='http://www.phyloxml.org', value='aggtcgcggcctgtggaagtcctctcct'),
Other(tag='seq', namespace='http://www.phyloxml.org', value='taaatcgc--cccgtgg-agtccc-cct')]
```

If you aren't interested in the "other" data, you can use parse() to
iteratively construct just the phylogenetic trees contained in the file
-- this is exactly the same as calling Phylo.parse() with the 'phyloxml'
format argument.

PhyloXMLIO.write() is similar to Phylo.write(), but also accepts a
Phyloxml object (the result of read() or to\_phyloxml()) to serialize.
Optionally, an encoding other than UTF-8 can be specified.

``` python
>>> phx = PhyloXMLIO.read('phyloxml_examples.xml')
>>> print phx.other
[Other(tag='alignment', namespace='http://example.org/align')]
>>> phx.other = []
>>> PhyloXMLIO.write(phx, 'ex_no_other.xml')
13
>>> phx_no = PhyloXMLIO.read('ex_no_other.xml')
>>> phx_no.other
[]
```

PhyloXMLIO also contains a utility called dump\_tags() for printing all
of the XML tags as they are encountered in a phyloXML file. This can be
helpful for debugging, or used along with grep or sort -u on the command
line to obtain a list of the tags a phyloXML file contains.

    >>> PhyloXMLIO.dump_tags('phyloxml_examples.xml')
    {http://www.phyloxml.org}phyloxml
    {http://www.phyloxml.org}phylogeny
    {http://www.phyloxml.org}name
    {http://www.phyloxml.org}description
    {http://www.phyloxml.org}clade
    ...

Using PhyloXML objects
----------------------

Standard Python syntactic sugar is supported wherever it's reasonable.

-   str() makes a string of the object's class name and an identifier,
    suitable for labeling a node in generated graph
-   repr() makes a string resembling the object constructor call, such
    that `eval(repr(obj))` will return `obj` for simpler PhyloXML
    objects, and at least partially rebuild more complex objects.
-   iter() is supported by Phyloxml and Clade objects, iterating over
    the contained phylogenies and sub-clades, respectively
-   len() is supported by the same objects that support iteration, with
    expected results

Clade objects also support slicing and multiple indexing:

``` python
tree = Phylo.parse('example.xml', 'phyloxml').next()
assert tree.clade[0] == tree.clade.clades[0]
assert tree.clade[0,1] == tree.clade.clades[0].clades[1]
```

Since valid Phylogeny objects always have a single clade attribute, this
style of indexing is a handy way to reach specific nodes buried deep in
the tree if you happen to know exactly where they are.

A couple of methods allow converting a selection to a new PhyloXML
object: Phylogeny.to\_phyloxml() and Clade.to\_phylogeny(). A few use
cases:

-   Parse a phyloXML containing multiple phylogenetic trees. Check each
    tree sequentially, and upon finding a tree with the desired
    characteristic, isolate it as a new PhyloXML object.

``` python
for tree in Phylo.parse('example.xml', 'phyloxml'):
    if tree.name == 'monitor lizards':
        return tree.to_phyloxml()
```

-   Extract a specific sub-clade and make it a separate phylogeny (and
    probably a new phyloXML file).

``` python
tree = Phylo.parse('example.xml', 'phyloxml').next()
best = None
for clade in tree.clade:
    if (clade.confidences[0].type == 'bootstrap'
            and (best is None
                or clade.confidences[0].value > best.confidences[0].value)):
        best = clade
phyloxml = best.to_phylogeny(rooted=True).to_phyloxml()
Phylo.write(phyloxml, 'example_best.xml', 'phyloxml')
```

### Core classes

Phyloxml

-   Container for phylogenies; not used by the top-level Bio.Phylo I/O
    functions

Phylogeny

-   Derived from Tree -- the global tree object

Clade

-   Derived from Subtree -- represents a node in the object tree, and
    local info

Other

-   Represents data included in the phyloXML file but not described by
    the phyloXML specification

### Annotation types

**(to do)**

### Integrating with the rest of Biopython

The classes used by this module inherit from the
[Phylo](Phylo "wikilink") module's generalized BaseTree classes, and
therefore have access to the methods defined on those base classes.
Since the phyloXML specification is very detailed, these subclasses are
kept in a separate module, Bio.Phylo.PhyloXML, and offer additional
methods for converting between phyloXML and standard Biopython types.

The PhyloXML.Sequence class contains methods for converting to and from
Biopython [SeqRecord](SeqRecord "wikilink") objects -- to\_seqrecord()
and from\_seqrecord(). This includes the molecular sequence (mol\_seq)
as a [Seq](Seq "wikilink") object, and the protein domain architecture
as list of [SeqFeature](SeqFeature "wikilink") objects. Likewise,
PhyloXML.ProteinDomain objects have a to\_seqfeature() method.

Performance
-----------

This parser is meant to be able to handle large files, meaning several
thousand external nodes. (Benchmarks of relevant XML parsers for Python
are [here](http://effbot.org/zone/celementtree.htm#benchmarks).) It has
been tested with files of this size; for example, the complete NCBI
taxonomy parses in about 100 seconds and consumes about 1.3 GB of
memory. Provided enough memory is available on the system, the writer
can also rebuild phyloXML files of this size.

The read() and parse() functions process a complete file in about the
same amount of CPU time. Most of the underlying code is the same, and
the majority of the time is spent building Clade objects (the most
common node type). For small files (smaller than
ncbi\_taxonomy\_mollusca.xml), the write() function serializes the
complete object back to an equivalent file slightly slower than the
corresponding read() call; for very large files, write() finishes faster
than read().

Here are some times on a 2.00GHz Intel Xeon E5405 processor (only 1 CPU
core used) with 7.7GB memory, running the standard Python 2.6.2 on
Ubuntu 9.04, choosing the best of 3 runs for each function:

| File                         | Ext. Nodes | Size (uncompressed) | Read (s) | Parse (s) | Write (s) |
|------------------------------|------------|---------------------|----------|-----------|-----------|
| apaf.xml                     |            | 38 KB               | 0.01     | 0.01      | 0.02      |
| bcl\_2.xml                   |            | 105 KB              | 0.02     | 0.02      | 0.04      |
| ncbi\_taxonomy\_mollusca.xml | 5632       | 1.5 MB              | 0.51     | 0.49      | 0.80      |
| tol\_life\_on\_earth\_1.xml  | 57124      | 46 MB               | 10.28    | 10.67     | 10.36     |
| ncbi\_taxonomy\_metazoa.xml  | 73907      | 33 MB               | 15.76    | 16.15     | 10.69     |
| ncbi\_taxonomy.xml           | 263691     | 31 MB (unindented)  | 109.70   | 109.14    | 32.39     |

On 32-bit architectures, [psyco](http://psyco.sourceforge.net/) might
improve these times significantly, at the risk of increasing memory
usage. (I haven't tested it.) For comparison, the Java-based parser used
in Forester and ATV (see below) reads the same files about 3-5 times as
quickly, or up to 15x for the largest file.

For Python 2.4, performance depends on which ElementTree implementation
is used. Using the original pure-Python elementtree, reading/parsing
takes about twice as much time for all file sizes, but writing is only
significantly slower for very large files.

Summer of Code project
----------------------

This module was developed by [Eric
Talevich](User%3AEricTalevich "wikilink") as a Google Summer of Code
2009 project to provide support for phyloXML in Biopython, with NESCent
as the mentoring organization and Brad Chapman and Christian Zmasek as
the mentors. The main page for the project is here: [PhyloSoC:Biopython
support for parsing and writing
phyloXML](https://www.nescent.org/wg_phyloinformatics/PhyloSoC:Biopython_support_for_parsing_and_writing_phyloXML)

The [Phylo](Phylo "wikilink") module was developed afterward in order to
integrate this code with the rest of Biopython.

Related software
----------------

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
