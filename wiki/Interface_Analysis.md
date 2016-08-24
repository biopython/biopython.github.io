---
title: Interface Analysis
permalink: wiki/Interface_Analysis
layout: wiki
---

### Description

The `Interface` module integrated in Biopython provides an easy and
friendly way to extract and analyze interface from PDB complexes.
Different information are calculated and provided with the extraction of
the interface: - polar/apolar/charged residues distribution - buried
surface area.

Developed during the Google Summer of Code 2011 by Mikael Trellet, this
module is still not available in the official repository of Biopython.
Nevertheless you can find the code (open-source) at the current link:
[GitHub
source](https://github.com/mtrellet/biopython/tree/interface_analysis)
In order to be relevant and complete, the `Interface` module works in
parallel with extended residues, a subclass of residue created during
the same period and also available at the previous link.

### Requirements

The main part of the `Interface` analysis module requires only a stable
installation of Biopython and Python 2.7. A more precise definition of
the interface can be done using the NACCESS module present in Biopython
but require a stable version of Naccess (available in
[NACCESS](http://www.bioinf.manchester.ac.uk/naccess/)).

### How to use

**Initialization**

Extraction of an interface is done from a complex PDB

``` python
from Bio.PDB import InterfaceBuilder

parser = PDBParser()
structure = parser.get_structure('test', '/home/directory/of/your/PDB/test.pdb')
```

Then the extraction of the interface is done in only one line

``` python
interface = InterfaceBuilder.InterfaceBuilder(structure[0]).get_interface()
```

From the `Interface` object some function and information are included

-   Get chains

``` python
chains = interface.get_chains()
for c in chains:
`   print(c)
```

-   Add a residue

``` python
interface.add(structure[0]['A'][24])
```

-   Get secondary structure distribution

``` python
ss = interface.secondary_structure
```

**Further own calculations**

Several statistics and information can be calculated from a single
interface

-   Polarity statistics

``` python
percent = interface.calculate_polarity()
```

-   Set up a dictionary of neighbors for each residues of the Interface

``` python
neighbors = interface.set_neighbors()
```

-   Calculate accessibility of Interface and the complex

``` pyhon
access = interface.calculate_accessibility()
```

**Comparison of 2 interfaces**

In order to compare 2 interfaces, so indirectly 2 complexes, few
calculations can be done:

-   Root Mean Square Deviation

``` python
rmsd = interface.rmsd(interface2)
```

-   Fraction of Common Contacts

``` python
fcc = interface.fcc(interface2)
```
