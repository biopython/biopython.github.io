---
title: TreeIO
permalink: wiki/TreeIO
layout: wiki
---

This module handles the parsing and generation of phylogenetic tree file
formats.

This code is not yet part of Biopython, and therefore the documentation
has not been integrated into the Biopython Tutorial yet either.

Installation
------------

This module is being developed alongside phyloXML support as a Google
Summer of Code 2009 project. To use TreeIO, see the
[PhyloXML](PhyloXML "wikilink") page for instructions on installing from
a [GitHub](GitUsage "wikilink") branch.

Usage
-----

Wrappers for supported file formats are available from the top level of
the module:

``` python
from Bio import TreeIO
```

This module provides three functions, read(), parse() and write(). Each
function accepts either a file name or an open file handle, so data can
be also loaded from compressed files, StringIO objects, and so on. If
the file name is passed as a string, the file is automatically closed
when the function finishes; otherwise, you're responsible for closing
the handle yourself.

The second argument to each function is the target format. Currently,
'phyloxml' and 'nexus' are supported.

See the [PhyloXML](PhyloXML "wikilink") page for more usage examples.
