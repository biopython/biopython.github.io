---
title: The Structure Modul
permalink: wiki/Struct
layout: wiki
---

This module extends `Bio.PDB`, providing additional features that prove
useful to structural biologists using Biopython.

### Availability

This module is only available (for now) in [João's GSOC2010 branch in
Github](http://github.com/JoaoRodrigues/biopython/tree/GSOC2010).

After getting the branch, it can be accessed easily by an import
statement. It should suffice to access all of the functions. Importing
specific parts of the module is not supposed to be necessary unless
likewise specific operations are sought after.

``` python
from Bio import Struct
```

### I/O Functions

`Bio.Struct` provides a simple I/O interface to read and write structures
in PDB format. It has two methods `read` and `write` that wrap
`Bio.PDB.PDBParser` and `Bio.PDB.PDBIO` respectively. It does not intend to
replace these methods, merely being a simpler way of performing I/O
tasks.

#### read()

``` python
s = Struct.read("protein_A.pdb")
```

The name of the resulting Structure object is based on the filename.

It has an additional *id* argument that allows specification of the
Structure object name, much like the first argument in
`PDBParser.get_structure()`.

#### write()

``` python
Struct.write(s)
```

The name of the resulting output is based on the Structure object id.

It has an additional *name* argument that allows naming the output file
(much like the main argument in `PDBIO.save()`). If a file already
exists by that name, the write() function automatically renames the
output adding `_0` (or `_1`, etc).

### Protein Class

We created a Structure-based class named `Protein` to confer protein
specific methods to a given SMCRA object. A method was also added to
`Bio.PDB.Structure` that allows easy interconversion between the two
classes: `as_protein()`. We will discuss each method in detail below.

#### as\_protein()

Converts a Structure object to the `Protein` class. The conversion also
filters all residues and excludes all those that are not aminoacids
(HETATMs are also excluded). This filtering can be disabled by setting
the optional `filter_residues` argument to `False`.

``` python
s = Struct.read("protein_A.pdb")

dir(s)
# Edited for shortening purposes
# ['__doc__', ... , 'renumber_residues', 'set_parent', 'xtra']

p = s.as_protein()

dir(p)
# ['__doc__', ... , 'renumber_residues', 'search_ss_bonds', 'set_parent', 'xtra']
```

#### search\_ss\_bonds()

The function returns an iterator with tuples of pairs of Cysteine
residues in close enough proximity to be forming a SS bond. The
threshold for a S-S contact to be defined as a SS bond is 3.0A, but it
can be manually specified through the `threshold` argument.

``` python
for bond in p.search_ss_bonds():
    print(bond)
# (<Residue CYS het=  resseq=5 icode= >, <Residue CYS het=  resseq=55 icode= >)
# (<Residue CYS het=  resseq=14 icode= >, <Residue CYS het=  resseq=38 icode= >)
# (<Residue CYS het=  resseq=30 icode= >, <Residue CYS het=  resseq=51 icode= >)
```

#### coarse\_grain()

Despite coarse-graining being general to all molecules, the current
implemented methods concern proteins only. To ease the introduction of
new CG-models, a `CG_models.py` class is present that defines how each
residue should be coarse grained. As of now, three models are supported.

##### Example Usage: MARTINI

``` python
cg_martini = p.coarse_grain("MARTINI")

for residue in cg_martini.get_residues():
    print(residue.resname + " " + residue.child_list)
# ARG [<Atom BB>, <Atom S1>, <Atom S2>]
# PRO [<Atom BB>, <Atom S1>]
# ASP [<Atom BB>, <Atom S1>]
# PHE [<Atom BB>, <Atom S1>, <Atom S2>, <Atom S3>]
# ...
# GLY [<Atom BB>]
# ALA [<Atom BB>]
```

##### Supported Models

###### CA Trace

Protein residues are reduced to their alpha carbon atom. This is the
default method.

###### ENCAD 3pt Model

ENCAD is the Energy Calculation and Dynamics program, developed by
Michal Levitt since the 80s. Its 3pt coarse grained model reduces each
protein residue to 3 points (some exceptions are reduced to 4 points):
Ca, O, and a side chain bead centered on a particular pre-defined atom.

###### MARTINI Model

MARTINI is a well known CG model that, as ENCAD, reduces protein
residues to 3/4 beads: Ca, O, and a side chain bead in a particular
position.

#### check\_missing\_atoms()

Compares the residues in the `Protein` object with a pre-defined topology
(derived from AMBER) to check for missing atoms. Outputs a dictionary of
tuples for each incomplete residue. Automatically ignores Hydrogen atoms
(`ha_only` argument can be set to `False` to override this) and allows
the usage of a particular template through the `template` argument
(default `None`). Templates are dictionaries with residue names as keys
and lists with atom names as values.

``` python
p.check_missing_atoms()
```

#### find\_seq\_homologues()

Bridges `Bio.Blast.NCBIWWW qblast()` function. It allows a direct sequence
homology search through that function using the `Protein` object's
sequence. Auto-adjusts parameters for short sequences. For more complex
homology searches, use the `Bio.Blast.NCBIWWW` module directly as this is
supposed to be just a convenience function.

It returns a list ranked by Expectation Value with some informational
values (e-value, identities, positives, gaps), the PDB code of the
match, and the alignment.

Allows an argument, `raw_output`, that replaces the default parsed
results with raw XML output from the BLAST search.

``` python
seq_homologues = p.find_seq_homologues()

for homologues in seq_homologues:
    print(homologues[0] + " " + homologues[1])
    print(homologues[-1])
    print()
# 2BUO 1.82482e-31
# DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW-TETLLVQNANPDCKTILKALGPGATLEE--TACQG
# DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW TETLLVQNANPDCKTILKALGPGATLEE  TACQG
# DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG
# ...
```
