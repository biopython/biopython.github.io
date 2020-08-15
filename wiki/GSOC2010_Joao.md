---
title: GSOC2010 Joao
permalink: wiki/GSOC2010_Joao
layout: wiki
---

Author & Mentors
----------------

[João Rodrigues](User%3AJoaor "wikilink") anaryin@gmail.com

**Mentors**


Eric Talevich

Diana Jaunzeikare

Peter Cock

Abstract
--------

Biopython is a very popular library in Bioinformatics and Computational
Biology. Its `Bio.PDB` module, originally developed by Thomas Hamelryck,
is a simple yet powerful tool for structural biologists. Although it
provides a reliable PDB parser feature and it allows several
calculations (Neighbour Search, RMS) to be made on macromolecules, it
still lacks a number of features that are part of a researcher's daily
routine. Probing for disulphide bridges in a structure and adding polar
hydrogen atoms accordingly are two examples that can be incorporated in
`Bio.PDB`, given the module's clever structure and good overall
organisation. Cosmetic operations such as chain removal and residue
renaming – to account for the different existing nomenclatures – and
renumbering would also be greatly appreciated by the community.

Another aspect that can be improved for `Bio.PDB` is a smooth
integration/interaction layer for heavy-weights in macromolecule
simulation such as MODELLER, GROMACS, AutoDock, HADDOCK. It could be
argued that the easiest solution would be to code hooks to these
packages' functions and routines. However, projects such as the recently
developed edPDB or the more complete [Biskit library](http://biskit.pasteur.fr/) render,
in my opinion, such interfacing efforts redundant. Instead, I believe it
to be more advantageous to include these software' input/output formats
in Biopython's `SeqIO` and `AlignIO` modules. This, together with the
creation of interfaces for model validation/structure checking
services/software would allow Biopython to be used as a pre- and
post-simulation tool. Eventually, it would pave the way for its
inclusion in pipelines and workflows for structure modelling, molecular
dynamics, and docking simulations.

Project Schedule
----------------

The schedule below was organised to be flexible, which means that some
features will likely be done early. Also, the weeks include
documentation and unit testing efforts for the features, with extended
periods for reviewing these efforts at the two points during the project
(halfway, final week).

### Community Bonding Period

-   Getting familiar with development environment (GitHub account, Git,
    Biopython's repository, Bug tracking system, etc)

-   Gather scientific literature and discuss some of the
    to-be-implemented methods.

### Week 1 \[31st May - 6th June\]

#### Renumbering residues of a structure

-   Read SEQRES record to account for gaps
-   Alternatively read ATOM records.

#### Probe disulphide bridges in the structure

-   Via `NeighbourSearch` class
-   Also use SSBOND in header

#### Extract Biological Unit

-   REMARK350 contains rotation and translation information
-   If REMARK is absent, do nothing.

### Week 2 \[7th – 13th June\]

#### Structure Hydrogenation

-   Add all/polar hydrogens through interface with WHATIF server.
-   Optionally define a set pH
    -   [pKa values
        algorithm](http://www3.interscience.wiley.com/journal/112117957/abstract)

#### Hydrogenation Report

-   Produces a brief list of polar hydrogen atoms in the structure.
    -   Chain \| Residue \[number\] \| Atom

### Weeks 3-5 \[14th June- 4th July\]

#### Removal of disordered atoms

-   [Solution proposed in the Biopython
    wiki](Remove_PDB_disordered_atoms "wikilink")

#### Residue name normalisation

-   Build conversion table from different nomenclatures (research them
    during c.bonding period )
-   Write function to make a given structure compliant with a given
    software nomenclature:
    -   Amber
    -   CNS/HADDOCK
    -   GROMACS

#### Coarse Grain Structure

-   Implement function to reduce complexity of a structure
    -   1pt\*c-alpha
    -   2pt\*c-alpha / c-beta
    -   3pt\*c-alpha / c-beta / side-chain pseudo-centroid OR side-chain
        centroid

### Week 6 (Mid-Term) \[5th - 11th July\]


-   Testing and consolidating the features thoroughly.

-   Write documentation & examples for each feature, to be included in
    Biopython's Wiki and `Bio.PDB`'s FAQ.

-   Mid-term Evaluations. Discussing with mentors current state of project
    and adjust following schedule to comply with project's needs.

### Week 7 \[12th - 19th July\]

#### Add support for MODELLER's PIR format to Biopython

-   [Format
    Description](http://www.salilab.org/modeller/manual/node495.html#alignmentformat)
-   `SeqIO`
-   `AlignIO`

#### Allow conversion of Structure Object to Sequence Object

-   Based on `Bio.PDB.Polypeptide` function

### Weeks 8-10 \[20th July - 9th August\]

#### Add Sequence/Structure Homology functions

-   Create call to Biopython's BLAST interfaces
    -   Allow direct blast from structure object ( e.g.
        `protein.find_homoseq()` )
    -   Returns list of tuples with E-Value \*Dictionary (name, length
        of alignment, etc..)
-   Create interface with structural homology web services
    -   e.g. [Dali
        server](http://ekhidna.biocenter.helsinki.fi/dali_server/)
    -   Return list of tuples with Z-Score\*Dictionary (name, etc...)

#### Implement basic structure validation checks

-   Via NeighbourSearch class
    -   Same Charge contacts
    -   Atom Clashes
-   Via ResidueDepth Class
    -   Buried Charges
-   Interface WHATIF PDBReport web service
    -   Parse WARNING and ERROR messages

### Week 11 \[10th - 17th August\]

#### Reviewing documentation, code, write tests for new functions.

Project Code
------------

Hosted at [this GitHub
branch](http://github.com/JoaoRodrigues/biopython/tree/GSOC2010)

Project Progress
----------------

Since I'm adding some methods that are useful/logical only for proteins,
having them exposed in `Structure.py` for every molecule could be
misleading. We decided then to add a `as_protein()` method that allows
protein-specific methods to be accessed. The following example
demonstrates how this call works. Note how the `search_ss_bonds`
method is absent from `dir(s)` but not from `dir(prot)`.

``` pycon
>>> from Bio.PDB import PDBParser
>>> p = PDBParser()
>>> s = p.get_structure('example', '4PTI.pdb')
>>> dir(s)  # Cut for viewing purposes
['__doc__', ... , 'renumber_residues', 'set_parent', 'xtra']
>>> prot = s.as_protein()
>>> dir(prot)
['__doc__', ... , 'renumber_residues', 'search_ss_bonds', 'set_parent', 'xtra']
```

### Renumbering residues of a structure

Since `parse_pdb_header` is far from optimal and is likely to change in
the future, I opted to forfeit reading SEQREQ records to account for
gaps. However, ignoring this information and renumbering based on ATOM
records would make us lose information on gaps. I opted to subtract the
first residue number-1 to all residues thus making the numbering start
in 1 and still keep gaps. I also added an argument (start) to allow the
user to set which number to start the counting from.

Example:

``` python
from Bio.PDB import PDBParser

p = PDBParser()
s = p.get_structure('example', '1IHM.pdb')

print(list(s.get_residues())[0])
# <Residue ASP het=  resseq=1029 icode= >

s.renumber_residues()
print(list(s.get_residues())[0])
# <Residue ASP het=  resseq=1 icode= >
```

### Probe disulphide bridges in the structure

The same rationale from SEQRES applies for the exclusion of looking up
SSBOND. Also, instead of using `NeighborSearch` to look for pairs of
cysteins in bond distance, I instead used the minus operator since it
has been overloaded to return the distance between two atoms (Page 10 of
the [FAQ](http://biopython.org/DIST/docs/tutorial/biopdb_faq.pdf)).
The average distance cited in the literature is 2.05A but other software
packages and my own tests set 3.0A as a good threshold. Still, the user
can set his own threshold manually.

The function returns an iterator with tuples of pairs of residues.

``` python
from Bio.PDB import PDBParser

p = PDBParser()
s = p.get_structure('example', '4PTI.pdb')

prot = s.as_protein()

for bond in prot.search_ss_bonds():
    print(bond)

# (<Residue CYS het=  resseq=5 icode= >, <Residue CYS het=  resseq=55 icode= >)
# (<Residue CYS het=  resseq=14 icode= >, <Residue CYS het=  resseq=38 icode= >)
# (<Residue CYS het=  resseq=30 icode= >, <Residue CYS het=  resseq=51 icode= >)
```

### Extract Biological Unit

Added parsing for REMARK350 to `parse_pdb_header` since there was
already a bit written for another REMARK section. This extracts the
transformation matrices and the translation vector from the header, that
is then fed to the `Structure` function. Each new rotated structure is
created as a new MODEL. I chose this because crystal structures very
rarely have more than one MODEL instance and also because NMR models
don't have REMARK 350 that often (at least to my knowledge).

``` python
from Bio.PDB import PDBParser

p = PDBParser()


s1 = p.get_structure('a', '4PTI.pdb')
s1.build_biological_unit()
# 'Processed 0 transformations on the structure.' # Identity matrix is ignored.

s2 = p.get_structure('b', 'homol_1bd8.pdb') # A homology model
s2.build_biological_unit()
# 'PDB File lacks appropriate REMARK 350 entries to build Biological Unit.'

s3 = p.get_structure('c', '1IHM.pdb')
s3.build_biological_unit()
# 'Processed 59 transformations on the structure.'
```

### Hydrogenation of PDB files

Following discussion between the mentors and me, we decided that maybe
it was better to not only include a webserver for this purpose but also
a local algorithm. This would not limit the user when there he/she lacks
an internet connection.

The interface for the WHATIF Protonation service has been implemented,
although it should be regarded as **highly experimental** for now.
Interfacing this server included writing a small parser for a
PDBXML-like format, which is expected to have serious bugs in its
initial versions. I ran some simple tests and it works. It doesn't
support water molecules yet, nor any other molecules other than
proteins. Such issues will be hopefully solved later on..

For those brave enough to want to test it (and help me debug it), here's
an example usage.

``` python
from Bio.Struct.WWW import WHATIF
from Bio import Struct

server = WHATIF.WHATIF() # Performs a sort of PING to the server. Gracefully exits if the servers are down.

# Get the protein structure
structure = Struct.read('4PTI.pdb')
protein = structure.as_protein() # This excludes water molecules

# Upload the structure to the WHATIF server
# This should convert the structure from a Structure object to a string via tempfile and PDBIO
# I was having some issues uploading structures...

id = server.UploadPDB(protein)

# Protonate
# Returns a Structure Object / WARNING! Bug prone for now.

protein_h = server.PDBasXMLwithSymwithPolarH(id)
```

Regarding the local implementation, after much reading I settled on
using PyMol's algorithm. It seems to allow for protonation of any
structure, regardless of its nature (protein, DNA, etc). Its vectorial
and matrix operations can likely be optimized with Numpy and Biopython's
`Vector` module. This first implementation works for proteins only.
I'll add general molecule support later.

``` python
from Bio import Struct
from Bio.Struct import Hydrogenate as H

s = Struct.read('1ctf.pdb')
p = s.as_protein()

prot = H.Hydrogenate_Protein()
prot.add_hydrogens(p)
```

### Coarse Grain Structure

A Center of Mass function was developed first as part of a new module
`Bio.Struct.Geometry`. It allows for calculation of the center of geometry
(all masses are equal) and center of mass (taking into account elemental
masses for the atoms). The masses are a new Atom object feature derived
from [this list](http://www.chem.qmul.ac.uk/iupac/AtWt/) and from PyMol.
Essentially, all atoms of a structure now get their mass defined when
the structure is created (check `Atom.py` and [this
thread](http://lists.open-bio.org/pipermail/biopython-dev/2010-June/007880.html)
for details). This is obviously experimental.

To calculate the center of mass of any Entity (Structure, Model, Chain,
Residue) or a List of Atoms:

``` pycon
>>> from Bio.Struct.Geometry import center_of_mass
>>> from Bio import Struct
>>> s = Struct.read('4PTI.pdb')
>>> print(center_of_mass.__doc__)
Returns gravitic or geometric center of mass of an Entity.
Geometric assumes all masses are equal (geometric=True)
Defaults to Gravitic.
>>> print(center_of_mass(s))
[14.833301303933874, 21.431581746366263, 4.1218478418007134]
>>> print(center_of_mass(s, geometric=True))
[14.805324902127458, 21.365571977563405, 4.1108949403803985]
```

As of now, 3 CG models are supported.

1) CA-Trace 2) ENCAD 3-point model (CA, O, Side Chain bead) 3) MARTINI
protein model (BB, Side Chain points \[S1 to S4\])

An example, picking up the s Structure from above:

``` pycon
>>> p = s.as_protein() # To expose the CG method
>>> ca_trace = p.coarse_grain()
>>> # One atom per residue
>>> print(len(list(p.get_residues())) == len(list(ca_trace.get_atoms())))
True
>>> cg_encad = p.coarse_grain('ENCAD_3P')
>>> for residue in cg_encad.get_residues():
...     print(residue.resname, residue.child_list)
...
ARG [<Atom CA>, <Atom O>, <Atom CMA>]
PRO [<Atom CA>, <Atom O>, <Atom CMA>]
ASP [<Atom CA>, <Atom O>, <Atom CMA>]
PHE [<Atom CA>, <Atom O>, <Atom CMA>]
CYS [<Atom CA>, <Atom O>, <Atom CMA>]
LEU [<Atom CA>, <Atom O>, <Atom CMA>]
GLU [<Atom CA>, <Atom O>, <Atom CMA>]
PRO [<Atom CA>, <Atom O>, <Atom CMA>]
PRO [<Atom CA>, <Atom O>, <Atom CMA>]
TYR [<Atom CA>, <Atom O>, <Atom CMA>]
...
CYS [<Atom CA>, <Atom O>, <Atom CMA>]
GLY [<Atom CA>, <Atom O>]
GLY [<Atom CA>, <Atom O>]
ALA [<Atom CA>, <Atom O>, <Atom CMA>]
>>> cg_martini = p.coarse_grain('MARTINI')
>>> for residue in cg_martini.get_residues():
...     print(residue.resname, residue.child_list)
...
ARG [<Atom BB>, <Atom S1>, <Atom S2>]
PRO [<Atom BB>, <Atom S1>]
ASP [<Atom BB>, <Atom S1>]
PHE [<Atom BB>, <Atom S1>, <Atom S2>, <Atom S3>]
CYS [<Atom BB>, <Atom S1>]
LEU [<Atom BB>, <Atom S1>]
GLU [<Atom BB>, <Atom S1>]
PRO [<Atom BB>, <Atom S1>]
PRO [<Atom BB>, <Atom S1>]
TYR [<Atom BB>, <Atom S1>, <Atom S2>, <Atom S3>]
......
CYS [<Atom BB>, <Atom S1>]
GLY [<Atom BB>]
GLY [<Atom BB>]
ALA [<Atom BB>]
```

### Removal of disordered atoms

Implement as part of `Structure.py` and based loosely on the [contribution
of Ramon
Crehuet](http://www.biopython.org/wiki/Remove_PDB_disordered_atoms). The
`DisorderedAtom` objects are removed from the residue and a single `Atom`
object is added corresponding to the location of the user's choice
(`keep_loc` argument) which defaults to A.

An example, still keeping s from above:

``` pycon
>>> s = s.remove_disordered_atoms(verbose=True)
0 residues were modified
>>> # Now if we load a structure with disordered atoms
>>> ds = Struct.read('1MC2.pdb')
>>> ds.remove_disordered_atoms(verbose=True)
Residue TRP:1010 has 8 disordered atoms: CD1/CD2/NE1/CE2/CE3/CZ2/CZ3/CH2
Residue VAL:1018 has 3 disordered atoms: CB/CG1/CG2
Residue LEU:1024 has 4 disordered atoms: CB/CG/CD1/CD2
Residue ARG:1043 has 7 disordered atoms: CB/CG/CD/NE/CZ/NH1/NH2
Residue MET:1092 has 4 disordered atoms: CB/CG/SD/CE
Residue ARG:1107 has 7 disordered atoms: CB/CG/CD/NE/CZ/NH1/NH2
Residue GLU:1108 has 4 disordered atoms: CG/CD/OE1/OE2
Residue ASP:1111 has 4 disordered atoms: CB/CG/OD1/OD2
Residue SER:1116 has 1 disordered atoms: OG
Residue SER:1131 has 1 disordered atoms: O
10 residues were modified
```

### Sequence Homologues from Structures

Biopython supports BLAST (local and remote through NCBI servers). We
bridged both `Bio.PDB` and `Bio.Blast` modules to allow an easier search for
sequence homologues. For now, it supports remote BLAST through
`Bio.Blast.NCBIWWW` and functions as a blackbox - i.e. users cannot change
any search parameter. If one wants to fully use BLAST he/she should use
the regular `BLAST` module. This is just a convenience function.

It is accessible only to `Protein` objects. It queries the PDB subset
database of NCBI BLAST servers with the Structure object's sequence,
auto-adjusting parameters for short sequences (less than 15 residues).

It returns a list ranked by Expectation Value with some informational
values (e-value, identities, positives, gaps), the PDB code of the
match, and the alignment.

``` pycon
>>> from Bio import Struct
>>> s = Struct.read('1A8O.pdb')
>>> p = s.as_protein()
>>> seq_homologues = p.find_seq_homologues()
>>> for homologues in seq_homologues:
...    print(homologues[0], homologues[1])
...    print(homologues[-1])
...    print()
...
2BUO 1.82482e-31
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW-TETLLVQNANPDCKTILKALGPGATLEE--TACQG
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW TETLLVQNANPDCKTILKALGPGATLEE  TACQG
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG

1AUM 1.82482e-31
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW-TETLLVQNANPDCKTILKALGPGATLEE--TACQG
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW TETLLVQNANPDCKTILKALGPGATLEE  TACQG
DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG
...
```

### Support for MODELLER PIR format in SeqIO

MODELLER PIR format support was added to `SeqIO` as 'pir-modeller'.
Currently, the format can be read but not written. An example of the
format follows, as well as an example of the parser's usage.

    >P1;5fd1
    structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19
    AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
    EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*

``` pycon
>>> from Bio import SeqIO
>>> for i in SeqIO.parse("test_pir.txt", "pir-modeller"):
...     print(i)
...
ID: 5fd1
Name: 5fd1
Description: ferredoxin
Number of features: 0
/r_factor= 0.19
/end_residue=106
/initial_chain=a
/end_chain=a
/record_type=X-Ray Structure
/initial_residue=1
/resolution= 1.90
/source_organism=Azotobacter vinelandii
Seq('AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAI...LER')
```
