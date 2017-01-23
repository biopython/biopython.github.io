---
title: Reading large PDB files with Biopython.
permalink: wiki/Reading_large_PDB_files
layout: wiki
tags:
 - Cookbook
---

If one is using Biopython to work with **PDB** files that have been
produced by molecular dynamics (MD) codes then one quickly runs into
problems with missing atoms on reading. The typical error message says
that atoms will be missing:

```
WARNING: Residue (' ', 1, ' ') redefined at line 31.
PDBConstructionException: Blank altlocs in duplicate residue SOL (' ', 1, ' ') at line 31.
Exception ignored.
Some atoms or residues will be missing in the data structure.
```

The problem is simply that these files can be large with hundreds of
thousands of atoms and residues (for instance, each water molecule is a
separate residue) and the PDB format has not enough space in the
appropriate columns of the
[ATOM](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM) or
[HETATM](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM)
record to accommodate atom numbers (*serial*) &gt;99,999 and residue
numbers (*resSeq*) &gt; 9999. Thus, these numbers are simply written
modulo 100,000 (*serial*) or modulo 10,000 (*resSeq*). This creates
duplicate entries in the chain and the `Bio.PDB.PDBParser` (or rather, the
StructureBuilder) complains. The effect is that not all atoms are read.

The [code below](#classes "wikilink") derives a new class from
`Bio.PDB.StructureBuilder.StructureBuilder` that simply increases the
*resSeq* if necessary. Because it is not very careful it is called a
`SloppyStructureBuilder`.

There's also a new class named `SloppyPDBIO` that writes pdb files with
*serial* and *resSeq* wrapped so that the resulting pdb files are [legal
PDB format](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html).

Example
-------

Load a big pdb file and write it out again. How to do interesting things
to the pdb file such as deleting some water molecules around a ligand is
left for another article.

``` python
from Bio.PDB import PDBParser
import xpdb   # this is the module described below

# read
sloppyparser = PDBParser(PERMISSIVE=True,
                         structure_builder=xpdb.SloppyStructureBuilder())
structure = sloppyparser.get_structure('MD_system', 'my_big_fat.pdb')

# ... do something here ...

# write
sloppyio = xpdb.SloppyPDBIO()
sloppyio.set_structure(structure)
sloppyio.save('new_big_fat.pdb')
```

Classes
-------

This is the Python implementation. Store it as a module `xpdb.py`
somewhere on your `PYTHONPATH`.

``` python
# xpdb.py -- extensions to Bio.PDB
# (c) 2009 Oliver Beckstein
# Relased under the same license as Biopython.
# See http://biopython.org/wiki/Reading_large_PDB_files

import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue


class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Cope with resSeq < 10,000 limitation by just incrementing internally.

    # Q: What's wrong here??
    #   Some atoms or residues will be missing in the data structure.
    #   WARNING: Residue (' ', 8954, ' ') redefined at line 74803.
    #   PDBConstructionException: Blank altlocs in duplicate residue SOL
    #   (' ', 8954, ' ') at line 74803.
    #
    # A: resSeq only goes to 9999 --> goes back to 0 (PDB format is not really
    #    good here)
    """

    # NOTE/TODO:
    # - H and W records are probably not handled yet (don't have examples
    #   to test)

    def __init__(self, verbose=False):
        Bio.PDB.StructureBuilder.StructureBuilder.__init__(self)
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        """Initiate a new Residue object.

        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for
            hetero residues, otherwise blanc.
        o resseq - int, sequence identifier
        o icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of
                # H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)

        if resseq > self.max_resseq:
            self.max_resseq = resseq

        if field == " ":
            fudged_resseq = False
            while (self.chain.has_id(res_id) or resseq == 0):
                # There already is a residue with the id (field, resseq, icode)
                # resseq == 0 catches already wrapped residue numbers which
                # do not trigger the has_id() test.
                #
                # Be sloppy and just increment...
                # (This code will not leave gaps in resids... I think)
                #
                # XXX: shouldn't we also do this for hetero atoms and water??
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)    # use max_resseq!
                fudged_resseq = True

            if fudged_resseq and self.verbose:
                sys.stderr.write("Residues are wrapping (Residue " +
                                 "('%s', %i, '%s') at line %i)."
                                 % (field, resseq, icode, self.line_counter) +
                                 ".... assigning new resid %d.\n"
                                 % self.max_resseq)
        residue = Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue = residue


class SloppyPDBIO(Bio.PDB.PDBIO):
    """PDBIO class that can deal with large pdb files as used in MD simulations

    - resSeq simply wrap and are printed modulo 10,000.
    - atom numbers wrap at 99,999 and are printed modulo 100,000

    """
    # The format string is derived from the PDB format as used in PDBIO.py
    # (has to be copied to the class because of the package layout it is not
    # externally accessible)
    _ATOM_FORMAT_STRING = "%s%5i %-4s%c%3s %c%4i%c   " + \
        "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname,
                       resseq, icode, chain_id, element="  ", charge="  "):
        """ Returns an ATOM string that is guaranteed to fit the ATOM format.

        - Resid (resseq) is wrapped (modulo 10,000) to fit into %4i (4I) format
        - Atom number (atom_number) is wrapped (modulo 100,000) to fit into
          %5i (5I) format

        """
        if hetfield != " ":
            record_type = "HETATM"
        else:
            record_type = "ATOM  "
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        x, y, z = atom.get_coord()
        bfactor = atom.get_bfactor()
        occupancy = atom.get_occupancy()
        args = (record_type, atom_number % 100000, name, altloc, resname,
                chain_id, resseq % 10000, icode, x, y, z, occupancy, bfactor,
                segid, element, charge)
        return self._ATOM_FORMAT_STRING % args


# convenience functions

sloppyparser = Bio.PDB.PDBParser(PERMISSIVE=True,
                                 structure_builder=SloppyStructureBuilder())


def get_structure(pdbfile, pdbid='system'):
    return sloppyparser.get_structure(pdbid, pdbfile)
```

See Also
--------

The `SloppyStructureBuilder()` was used as the basis for a small Python
module
[edPDB](http://sbcb.bioch.ox.ac.uk/users/oliver/software/GromacsWrapper/epydoc/edPDB-pysrc.html)
to edit PDB files in preparation for MD simulations.
