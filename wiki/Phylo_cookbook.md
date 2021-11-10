---
title: Bio.Phylo Cookbook.
permalink: wiki/Phylo_cookbook
layout: wiki
tags:
 - Cookbook
---

Here are some examples of using [`Bio.Phylo`](Phylo "wikilink") for some
likely tasks. Some of these functions might be added to Biopython in a
later release, but you can use them in your own code with Biopython
1.54.

Convenience functions
---------------------

### Get the parent of a clade

The Tree data structures in `Bio.Phylo` don't store parent references for
each clade. Instead, the `get_path` method can be used to trace the path
of parent-child links from the tree root to the clade of choice:

``` python
def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]


# Select a clade
myclade = tree.find_clades("foo").next()
# Test the function
parent = get_parent(tree, myclade)
assert myclade in parent
```

Note that `get_path` has a linear run time with respect to the size of
the tree -- i.e. for best performance, don't call `get_parent` or
`get_path` inside a time-critical loop. If possible, call `get_path`
outside the loop, and look up parents in the list returned by that
function.

Alternately, if you need to repeatedly look up the parents of arbitrary
tree elements, create a dictionary mapping all nodes to their parents:

``` python
def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents


# Example
parents = all_parents(tree)
myclade = tree.find_clades("foo").next()
parent_of_myclade = parents[myclade]
assert myclade in parent_of_myclade
```

### Index clades by name

For large trees it can be useful to be able to select a clade by name,
or some other unique identifier, rather than searching the whole tree
for it during each operation.

``` python
def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names
```

Now you can retrieve a clade by name in constant time:

``` python
tree = Phylo.read("ncbi_taxonomy.xml", "phyloxml")
names = lookup_by_names(tree)
for phylum in ("Apicomplexa", "Euglenozoa", "Fungi"):
    print("Phylum size: %d" % len(names[phylum].get_terminals()))
```

A potential issue: The above implementation of `lookup_by_names` doesn't
include unnamed clades, generally internal nodes. We can fix this by
adding a unique identifier for each clade. Here, all clade names are
prefixed with a unique number (which can be useful for searching, too):

``` python
def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = "%d_%s" % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names
```

### Calculate distances between neighboring terminals

*Suggested by Joel Berendzen*

``` python
import itertools


def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""

    def generate_pairs(self):
        pairs = itertools.tee(self)
        pairs[1].next()
        return itertools.izip(pairs[0], pairs[1])

    return [self.distance(*i) for i in generate_pairs(self.find_clades(terminal=True))]
```

### Test for "semi-preterminal" clades

*Suggested by Joel Berendzen*

The existing tree method `is_preterminal` returns `True` if all of the
direct descendants are terminal. This snippet will instead return `True`
if *any* direct descendent is terminal, but still `False` if the given
clade itself is terminal.

``` python
def is_semipreterminal(clade):
    """True if any direct descendent is terminal."""
    for child in clade:
        if child.is_terminal():
            return True
    return False
```

In Python 2.5 and later, this is simplified with the built-in `any`
function:

``` python
def is_semipreterminal(clade):
    return any(child.is_terminal() for child in clade)
```

Comparing trees
---------------

*TODO:*

-   Symmetric difference / partition metric, a.k.a. topological
    distance (Robinson-Foulds)
-   Quartets distance
-   Nearest-neighbor interchange
-   Path-length-difference

Consensus methods
-----------------

*TODO:*

-   Majority-rules consensus
-   Strict consensus
-   Adams ([Adams 1972](http://dx.doi.org/10.2307/2412432),
    [pdf](http://www.faculty.biol.ttu.edu/Strauss/Phylogenetics/Readings/Adams1972.pdf))
-   Asymmetric median tree
    ([Phillips and Warnow 1996](http://dx.doi.org/10.1007/3-540-61258-0_18))

Rooting methods
---------------

The basic method on the Tree class (not TreeMixin) is
`root_with_outgroup`:

``` python
tree = Phylo.read("example.nwk", "newick")
print(tree)
# ...
tree.root_with_outgroup({"name": "A"})  # Operates in-place
print(tree)
```

Normally you'll want the outgroup to be a monophyletic group, rather
than a single taxon. This isn't automatically checked, but you can do it
yourself with the `is_monophyletic` method.

To save some typing, try keeping the query in a separate list and
reusing it:

``` python
outgroup = [{"name": taxon_name} for taxon_name in ("E", "F", "G")]
if tree.is_monophyletic(outgroup):
    tree.root_with_outgroup(*outgroup)
else:
    raise ValueError("outgroup is paraphyletic")
```

*TODO:*

-   Root at the midpoint between the two most distant nodes (or "center"
    of all tips)

Graphics
--------

*TODO:*

-   Party tricks with draw methods, covering each keyword argument

Exporting to other types
------------------------

### Convert to an 'ape' tree, via Rpy2

The **R** statistical programming environment provides support for
phylogenetics through the [**ape**](http://ape-package.ird.fr/) package and
several others that build on top of **ape**. The Python package
[rpy2](http://rpy2.bitbucket.org/) provides an interface
between R and Python, so it's possible to convert a `Bio.Phylo` tree into
an **ape** tree object:

``` python
import tempfile
from rpy2.robjects import r


def to_ape(tree):
    """Convert a tree to the type used by the R package `ape`, via rpy2.

    Requirements:
        - Python package `rpy2`
        - R package `ape`
    """
    with tempfile.NamedTemporaryFile() as tmpf:
        Phylo.write(tree, tmpf, "newick")
        tmpf.flush()
        rtree = r(
            """
            library('ape')
            read.tree('%s')
            """
            % tmpf.name
        )
    return rtree
```

See that it works:

``` pycon
>>> from StringIO import StringIO
>>> from Bio import Phylo
>>> tree = Phylo.read(StringIO("(A,(B,C),(D,E));"), "newick")
>>> rtree = to_ape(tree)
>>> len(rtree)
3
>>> print(r.summary(rtree))
Phylogenetic tree: structure(list(edge = structure(c(6, 6, 7, 7, 6, 8, 8, 1, 7,  2, 3, 8, 4, 5),
 .Dim = c(7L, 2L)), tip.label = c("A", "B", "C",  "D", "E"), Nnode = 3L),
 .Names = c("edge", "tip.label", "Nnode" ), class="phylo")

  Number of tips: 5
  Number of nodes: 3
  No branch lengths.
  No root edge.
  Tip labels: A
              B
              C
              D
              E
  No node labels.
NULL
>>> r.plot(rtree)
```

See the **rpy2** documentation for further guidance.

### Convert to a DendroPy or PyCogent tree

The tree objects used by Biopython, DendroPy and PyCogent are different.
Nonetheless, all three toolkits support the Newick file format, so
interoperability is straightforward at that level by writing to a
temporary file or `StringIO` object with one library, then reading the
same string again with another.

``` python
from Bio import Phylo
import cogent

Phylo.write(bptree, "mytree.nwk", "newick")  # Biopython tree
ctree = cogent.LoadTree("mytree.nwk")  # PyCogent tree
```

``` python
import dendropy

# Create or load a tree in DendroPy
dtree = dendropy.Tree.get_from_string("(A, (B, C), (D, E))", "newick")
dtree.write_to_path("tmp.nwk", "newick", suppress_rooting=True)
# Load the same tree in Biopython
bptree = Phylo.read("tmp.nwk", "newick")
```

### Convert to a NumPy array or matrix

**Adjacency matrix:** cells are 1 (true) if a parent-child relationship
exists, otherwise 0 (false).

``` python
import numpy


def to_adjacency_matrix(tree):
    """Create an adjacency matrix (NumPy array) from clades/branches in tree.

    Also returns a list of all clades in tree ("allclades"), where the position
    of each clade in the list corresponds to a row and column of the numpy
    array: a cell (i,j) in the array is 1 if there is a branch from allclades[i]
    to allclades[j], otherwise 0.

    Returns a tuple of (allclades, adjacency_matrix) where allclades is a list
    of clades and adjacency_matrix is a NumPy 2D array.
    """
    allclades = list(tree.find_clades(order="level"))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    adjmat = numpy.zeros((len(allclades), len(allclades)))
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            adjmat[lookup[parent], lookup[child]] = 1
    if not tree.rooted:
        # Branches can go from "child" to "parent" in unrooted trees
        adjmat = adjmat + adjmat.transpose()
    return (allclades, numpy.matrix(adjmat))
```

**Distance matrix:** cell values are branch lengths if a branch exists,
otherwise infinity (this plays well with graph algorithms).

``` python
import numpy


def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.

    A cell (i,j) in the array is the length of the branch between allclades[i]
    and allclades[j], if a branch exists, otherwise infinity.

    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a NumPy 2D array.
    """
    allclades = list(tree.find_clades(order="level"))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    distmat = numpy.repeat(numpy.inf, len(allclades) ** 2)
    distmat.shape = (len(allclades), len(allclades))
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            if child.branch_length:
                distmat[lookup[parent], lookup[child]] = child.branch_length
    if not tree.rooted:
        distmat += distmat.transpose
    return (allclades, numpy.matrix(distmat))
```

Enhancements:

-   Use an `OrderedDict` for `allclades`, so the separate dictionary
    `lookup` isn't needed. (Python 2.7+)
-   Use NumPy's [structured arrays](http://docs.scipy.org/doc/numpy-1.10.1/user/basics.rec.html) to
    assign clade names to rows and columns of the matrix, so `allclades`
    isn't needed either (this works nicely along with the
    `tabulate_names` function given earlier).

*TODO:*

-   Relationship matrix? See [Martins and Housworth
    2002](http://dx.doi.org/10.1080/10635150290102573)

