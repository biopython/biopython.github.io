---
title: GSoC2011 mtrellet
permalink: wiki/GSoC2011_mtrellet
layout: wiki
redirect_from:
 - /wiki/User:Mtrellet
---

Author & Mentors
----------------

[Mikael Trellet](User%3AMtrellet "wikilink") mikael.trellet@gmail.com

**Mentors**

  
João Rodrigues

Eric Talevich

Interface analysis wikipage
---------------------------

You can find a descriptive and a "how to use" part of the Interface
analysis module at the following link:

[Interface analysis](http://biopython.org/wiki/Interface_Analysis)

Abstract
--------

Analysis of protein-protein complexes interfaces at a residue level
yields significant information on the overall binding process. Such
information can be broadly used for example in binding affinity studies,
interface design, and enzymology. To tap into it, there is a need for
tools that systematically and automatically analyze protein structures,
or that provide means to this end. Protorop
(http://www.bioinformatics.sussex.ac.uk/protorp/) is an example of such
a tool and the elevated number of citations the server has had since its
publication acknowledge its importance. However, being a webserver,
Protorop is not suited for large-scale analysis and it leaves the
community dependent on its maintainers to keep the service available. On
the other hand, Biopython’s structural biology module, Bio.PDB, provides
the ideal parsing machinery and programmatic structures for the
development of an offline, open-source library for interface analysis.
Such a library could be easily used in large-scale analysis of
protein-protein interfaces, for example in the CAPRI experiment
evaluation or in benchmark statistics. It would be also reasonable, if
time permits, to extend this module to deal with protein-DNA or
protein-RNA complexes, as Biopython supports nucleic acids already.

Project Schedule
----------------

### Week 1 \[23rd May - 31st June\]

#### Add the new module backbone in current Bio.PDB code base

-   Evaluate possible code reuse and call it into the new module
-   Try simple calculations to be sure that there is stability between
    the different modules (parsing for example) and functions

#### Define a stable benchmark

-   Select few PDB files among interface size and proteins size would be
    different
-   Add some basics unit tests

### Weeks 2-3 \[1st - 13th June\]

#### Extend IUPAC.Data module with residue information

-   Deduce residues weight from Atom instead of direct dictionnary
    storage
-   Polar/charge character (dictionary or influenced by pH)
-   Hydrophobicity scale(s)
    -   [pKa values
        algorithm](http://www3.interscience.wiley.com/journal/112117957/abstract)

### Weeks 4 \[14th - 21st June\]

#### Implement Extended Residue class as a subclass of Residue

-   Build Extended Residue on the fly or have it hard-coded (?)
-   Allow regular operations on Residue to be performed seamlessly in
    Extended Residue (should come with inheritance)
-   Unit tests on pdb files containing particular residues

### Weeks 5-6-7 \[22nd June - 11th July\]

#### Implement InterfaceAnalysis module

-   Develop Interface class as a subclass of Model
-   Develop method to automatically extract Interface from parsed
    structure upon class instantiation
    -   e.g. I = Interface(Structure)
    -   Allow threshold for distance
    -   Allow chain pairs to ignore (to avoid intra-molecule contacts)
-   Unit tests with results from usual scripts, broadly used by
    scientists

### Mid term evaluation

### Weeks 7-8 \[12th July - 25th July\]

#### Develop functions for interface analysis

-   Calculation of interface polar character statistics (% of polar
    residues, apolar, etc)
-   Calculation of BSA calling MSMS or HSA
-   Calculation of SS element statistics in the interface through DSSP
-   ...
-   Unit tests and use of results as input for further calculations by
    other tools and scripts

### Weeks 9-10 \[26th July - 8th August\]

#### Develop functions for Interface comparison

-   Perhaps adapt current RMSD functions to allow usage of Interface
    Residues
-   Otherwise, should be called through something like Ia.rmsd\_to(Ib)
    where Ia and IB are interface objects
-   Calculation of iRMSD
-   Calculation of FCC (Fraction of Common Contacts)
-   Rough Identity and Similarity percentage
-   ...
-   Unit tests, comparison with specific tools as Profit

### Weeks 11 \[9th July - 8th August\]

#### Code organization and final testing

Unit tests
----------

Unit tests will be perfomed along the project, allowing to do only a
larger test at the end gathering every tests already performed.

Then the aim will be to optimized, if possible, some parts of the code
in efficiency and rapidity without changes at algorithmic level. Several
days will be booked to package code and be sure that everything can
communicate with Biopython.

------------------------------------------------------------------------

Project Progress
----------------

### Implementation of Interface object backbone

-   Theory

We began to think of an easy way to add the Interface as a new part of
the SMCRA scheme. The idea was to have this new scheme = SM-I-CRA.
Unfortunately the Interface object is not as well defined as just a
child of model and a parent of chains. Indeed, the main part of the
interface is residues, and even residues pairs. We want to keep the
information of the chain but we can't keep them as they are defined
actually, since we will get some overlaps, duplication and
miscompatibility between the chains of our model and the chains of our
interface. In the same way, our try to link the creation of the
interface with existing modules as StructureBuilder and Model wasn't
successful. So, we decided to simplify a bit the concept in adding the
classes related to the Interface in an independent way. Obviously links
will exist between the different levels of SMCRA but Interface would be
considered now as a parallel entity, not integrated completely in the
SMCRA scheme.

-   Coding

Interface.py is the definition of the Interface object inherited from
Entity with the following methods : **\_\_init\_\_**(self, id),
**add**(self, entity) and **get\_chains**(self).

The add module overrides the add method of Entity in order to have an
easy way to class residues according to their respective chains. The
get\_chains modules returns the chains involved in the interface defined
by the Interface object.

The second class created is InterfaceBuilder.py which deals directly
with the interface building (hard to guess..!) We find these different
modules : **\_\_init\_\_**(self, model, id=None, threshold=5.0,
include\_waters=False, \*chains), **\_unpack\_chains**(self,
list\_of\_tuples), **get\_interface**(self), **\_add\_residue**(self,
residue), **\_build\_interface**(self, model, id, threshold,
include\_waters=False, \*chains)

**\_\_init\_\_** : In order to initialize an interface you need to
provide the model for which you want to calculate the interface, that's
the only mandatory argument.

**\_unpack\_chains**: Method used by \_\_init\_\_ so as to create
self.chain\_list, variable read in many parts of the class. It
transforms a list of tuples (given by the user) in a list of characters
representing the chains which will be involved in the definition of the
interface.

**get\_interface**: Returns simply the interface

**\_add\_residue**: Allows the user to add some specific residues to his
interface

**\_build\_interface**: The machinery to build the interface, it uses
NeighborSearch and Selection in order to define the interface depending
on the arguments given by the user.

-   Github repository

[Interface.py](https://github.com/mtrellet/biopython/commit/4cfa4359d0f927609c076ed7b66f37add5aabdfb)
[InterfaceBuilder.py](https://github.com/mtrellet/biopython/commit/194efe37ac8f88d688e0cf528f1fb896c8441866)

### Extension of residues

-   Theory

In order to have several useful information about residues and in the
aim to use them for further calculations inside interfaces, we want to
implement a subclass of residues. This new class would be integrate as
an inherited class of Residues (as DisorderedResidues) and would
calculate by default few information about residues as polar charge,
hydrophobicity or weight. We first thought to create a copy of an
existing residue to create an extended residue but, with regards to the
memory consumming, we preferred change directly the residue type. Thus,
we change the class of the residue and perform calculation during the
initialization of the new ExtendedResidue class. So, changing of a
residue is done by Residue.extend()

-   Coding

ExtendedResidue is an inherited class from Residue, these are the
different methods inside :

**\_\_init\_\_**: Sets charge, mass and hydrophobicity of a residue.

**\_\_repr\_\_**: Defines a new id for the residue, based on the residue
id.

**set\_hydrophobicity**: Sets the residue hydrophobicoty depending on
the scale choice, KD {Kyle and Doolittle}, OHM {Sweet and Eisenberg} or
consensus (by default= consensus). Returns the scale chosen.

**set\_charge**: Sets the relative charge of a residue with regards to
its pka and the pH defined by the user (by default 7.0)

**set\_mass**: Sets the mass of a residue according to its atoms mass

### Interface information

-   Theory

In order to work with complexes interface, few important information are
needed. Some of them are calculated during the initialization of the
interface and others are calculated via functions.

-   Coding

Private functions in InterfaceBuilder (used during the Interface
initialization):

**\_secondary\_structure**: Calculates secondary structure of each
residue of the Interface using DSSP.

**\_rsa\_calculation**: Used to define the interface according to the
residues accessibility.

Public functions in Interface (used by the user):

**set\_neighbors**:Creates a dictionary of neighbors for each residue of
the interface.

**calculate\_accessibility**:Calculates the buried surface area of the
complex and indirectly the accessibility of the interface. Total
accessibility and accessibility of the chains are also calculated in
this function. This calculation is done via the NACCESS module.

### Interfaces comparison

-   Theory

Several ways are used by scientists to measure differences between two
complexes. In particular in structural biology, an important part of the
work is the analyse of the models generated by a scientific method.
Reliability and accurately of a docking for instance. To know if a model
can be considered as acceptable, a reference is taken and from this
reference, some calculations are done, comparing models with reference.
We can cite RMSD or Fraction of native contacts for example. These are
two calculations we added in the Interface analysis module.

-   Coding

**interface\_ref.rmsd(interface\_model)**: This function used
Superimposer module in order to calculate Root Mean Square Deviation
between residues of 2 interfaces. It gives the interface RMSD as well.

**interface\_ref.fnat(interface\_model)**: This function calculates the
fraction of common contacts between 2 interfaces.
