---
title: Remove PDB disordered atoms with the Bio.PDB module.
permalink: wiki/Remove_PDB_disordered_atoms
layout: wiki
tags:
 - Cookbook
---

*Contributed by Ramon Crehuet*

Problem
-------

You have a PDB with disordered atoms, i.e. different atomic positions
with occupancies that add up to 100%. From this PDB you want to create a
new one having only one set of the disordered atoms. This can be
necessary if you want to perform RMSD calculations or Molecular Dynamics
simulations.

Solution
--------

[`Bio.PDB`](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
is proficient in dealing with disordered atoms. Each disordered atom has
a property indicating its alternative positions: `atom.altloc`. Usually
there are only two alternative positions labelled *'A'* and *'B'*. The key
is to save a PDB with the optional `select` argument. This argument
needs to return a `True` value for the atoms that have to be saved. In the
following example we save all not-disordered atoms and the *'A'* positions
of the disordered ones.

``` python
from Bio.PDB import *

parser = PDBParser()
s = parser.get_structure('my_pdb', 'my_pdb.pdb')
io = PDBIO()


class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == 'A'

io = PDBIO()
io.set_structure(s)
io.save("ordered.pdb", select=NotDisordered())
```

Note that the code above does not eliminate the alternate location
identifier (*'A'* in the example above). It is the programmer's
responsibility to eliminate the identifier when necessary.

``` python
keepAltID = ...


class NMROutputSelector2(Select):  # Inherit methods from Select class
    def accept_atom(self, atom):
        if (not atom.is_disordered()) or atom.get_altloc() == keepAltID:
            atom.set_altloc(' ')  # Eliminate alt location ID before output.
            return True
        else:  # Alt location was not one to be output.
            return False
        # end of accept_atom()
# end of NMROutputSelector2()
```

Discussion
----------

It is trivial to change that to save *'B'* altloc positions. One can even
do more complicated selections based on other atom properties. The key
is to generate a class that returns `True` or `False` for a given atom. One
could also think of deleting atoms with *'B'* values in `atom.altloc`.

``` python
# Will not work!

for atom in all_atoms:  # all_atoms is a list containg all atoms
   if atom.altloc == 'B':
       del atom
```

but that **does not work**, because it only deletes the local variable
and not the PDB structure.
