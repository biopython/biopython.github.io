---
title: Google Summer of Code
permalink: wiki/Google_Summer_of_Code
layout: wiki
---

Introduction
------------

For the past several years, Biopython-related GSoC projects have been
successfully run under the mentorship of the [Open Bioinformatics
Foundation (OBF)](http://www.open-bio.org/wiki/Google_Summer_of_Code)
and the [National Evolutionary Synthesis Center
(NESCent)](http://informatics.nescent.org/wiki/Phyloinformatics_Summer_of_Code_2013).

Please read those organizations' GSoC pages and the [main Google Summer
of Code page](http://code.google.com/soc) for more details about the
program.

Current project proposals
-------------------------

See the [Open Bioinformatics Foundation (OBF) GSoC wiki
page](http://www.open-bio.org/wiki/Google_Summer_of_Code) and [OBF GSoC
page](http://obf.github.io/GSoC/) as our usual mentoring organization.

Any project ideas for Biopython can be posted at
<http://obf.github.io/GSoC/ideas/> (update
<https://github.com/OBF/GSoC/blob/gh-pages/00_ideas.md> via a GitHub
pull request).

We encourage potential students and mentors to join the [BioPython
mailing lists](http://biopython.org/wiki/Mailing_lists) and actively
participate in developing these project ideas, either by submitting
their own ideas or contributing to improving existing ones.

Past Mentors
------------

Usually, each BioPython proposal has one or more mentors assigned to it.
Nevertheless, we encourage potential students/mentors to contact the
[mailing lists](http://biopython.org/wiki/Mailing_lists) with their own
ideas for proposals. There is therefore not a set list of 'available'
mentors, since it highly depends on which projects are proposed every
year.

Past mentors include:

-   [James Casbon](http://casbon.me/)
-   [Brad Chapman](https://github.com/chapmanb)
-   [Peter Cock](http://www.hutton.ac.uk/staff/peter-cock)
-   [Thomas Hamelryck](http://wiki.binf.ku.dk/User:Thomas_Hamelryck)
-   [Reece Hart](http://www.linkedin.com/in/reece)
-   [João Rodrigues](http://nmr.chem.uu.nl/~joao)
-   [Eric Talevich](http://etal.myweb.uga.edu/)

Past Accepted Projects
----------------------

### 2014 (OBF)

#### Indexing & Lazy-loading Sequence Parsers

Student  
Evan Parker ([blog](http://blog.evanaparker.com/))

Rationale  
[Bio.SeqIO](SeqIO "wikilink")'s indexing offers parsing on demand access
to any sequence in a large file (or collection of files on disk) as a
[SeqRecord](SeqRecord "wikilink") object. This works well when you have
many small to medium sized sequences/genomes. However, this is not ideal
for large genomes or chromosomes where only a sub-region may be needed.
A lazy-loading parser would delay reading the record until requested.
For example, if region *record\[3000:4000\]* is requested, then only
those 1000 bases need to be loaded from disk into memory, plus any
features in that region. This is how Biopython's
[BioSQL](BioSQL "wikilink") interface works. Tools like tabix and
samtools have demonstrated efficient co-ordinate indexing which could be
useful here.

Aside from being used via an index for random access, lazy-loading
parsers could be used when iterating over a file as well. This can
*potentially* offer speed ups for tasks where only a fraction of the
data is used. For example, if calculating the GC content of a collection
of genomes from GenBank, using Bio.SeqIO.parse(...) would currently
needlessly load and parse all the annotation and features. A lazy-parser
would only parse the sequence information.

Approach & Goals  
Useful features include:

-   Internal indexing of multiple file formats, including FASTA and
    richly annotated sequence formats like GenBank/EMBL
    and GTF/GFF/GFF3.
-   Full compatibility with existing SeqIO parsers which load everything
    into memory as a \`SeqRecord\` object.

Difficulty and needed skills  
Hard. Familiarity with the Biopython's existing sequence
parsing essential. Understanding of indexing large files will be vital.

Possible Mentors  
[Wibowo Arindrarto](https://github.com/bow), [Peter
Cock](https://github.com/peterjc/), others welcome

### 2013 (NESCent)

This year the Open Bioinformatics Foundation was not accepted on the
very competitive GSoC programme. Biopython instead participated under
[NEScent](http://informatics.nescent.org/wiki/Phyloinformatics_Summer_of_Code_2013).

#### Codon alignment and analysis

Student  
Zheng Ruan ([blog](http://zruanweb.com/))

Rationale  
A codon alignment is an alignment of nucleotide sequences in which the
trinucleotides correspond directly to amino acids in the translated
protein product. This carries important information which can be used
for several analysis, notably analysis of selection pressures. This
project extends Biopython to support this data type and these analyses.

Approach & Goals  
Useful features include:

-   Conversion of a set of unaligned nucleic acid sequences and a
    corresponding protein sequence alignment to a codon alignment.
-   Calculation of selection pressure from the ratio of nonsynonomous to
    synonomous site replacements, and related functions.
-   Model selection.
-   Possible extension of [AlignIO](AlignIO "wikilink") and the
    MultipleSeqAlignment class to take full advantage of codon
    alignments, including validation (testing for frame shifts, etc.)

Difficulty and needed skills  
Medium, depending on ambition. Familiarity with the Biopython's existing
alignment classes and functions, or equivalents in BioPerl, BioJava or
BioRuby (e.g.), will be helpful. Understanding of the practical uses of
codon alignments, or at least a basic understanding of molecular
biology, is important. Some basic math is involved, essentially reading
a few equations and converting them to code. One useful book to have on
hand is *[Computational Molecular
Evolution](http://abacus.gene.ucl.ac.uk/CME2006/)* by Ziheng Yang.

Mentors  
[Eric Talevich](http://etal.myweb.uga.edu/), [Peter
Cock](https://github.com/peterjc/)

#### Bio.Phylo: filling in the gaps

Student  
Yanbo Ye ([blog](http://blog.yeyanbo.com/tag/gsoc.html))

Rationale  
While the [Phylo](Phylo "wikilink") module in Biopython supports I/O and
basic tree operations, there are some important components that remain
to be implemented to better support phylogenetic workflows.

Approach & Goals  
This "idea" is intentially left open-ended -- some potentially useful
features are:

-   A "Phylo.consensus" module with functions for the consensus of
    multiple trees. E.g.: Strict consensus, as Bio.Nexus already
    implements; perhaps other methods like Adams consensus.
-   A function to calculate bootstrap support for branches given a
    target or "master" tree and a series of bootstrap replicate trees
    (usually read incrementally from a Newick file).
-   Functions for comparison of two or more trees.
-   Simple algorithms for tree inference, i.e. neighbor-joining and
    parsimony tree estimation. For small alignments (and perhaps
    medium-sized ones with PyPy), it would be nice to run these without
    an external program, e.g. to construct a guide tree for another
    algorithm or quickly view a phylogenetic clustering of sequences.
-   Tree visualizations: A proper draw\_unrooted function to perform
    radial layout, with an optional "iterations" argument to use
    Felsenstein's Equal Daylight algorithm. Circular radial diagrams are
    also popular these days. Any new function should support the same
    arguments as the existing Phylo.draw.

Difficulty and needed skills  
Medium-to-advanced programming skill in Python -- it's important for
these implementations to be reasonably efficient, though we don't aim to
compete with the fastest stand-alone implementations of
these algorithms. Knowledge of phylogenetic methods is critical; for
reference, you might like to have a copy of Joe Felsenstein's
*[Inferring Phylogenies](http://www.sinauer.com/detail.php?id=1775)*.
Tree visualizations are done with matplotlib.

Mentors  
Mark Holder, Jeet Sukumaran, [Eric Talevich](http://etal.myweb.uga.edu/)

### 2012 (OBF)

#### [SearchIO](http://biopython.org/wiki/SearchIO)

Student  
Wibowo Arindrarto ([blog](http://bow.web.id/blog/2012/08/back-on-the-main-branch/))

Rationale  
Biopython has general APIs for parsing and writing assorted sequence
file formats (SeqIO), multiple sequence alignments (AlignIO),
phylogenetic trees (Phylo) and motifs (Bio.Motif). An obvious omission
is something equivalent to BioPerl's SearchIO. The goal of this proposal
is to develop an easy-to-use Python interface in the same style as
SeqIO, AlignIO, etc but for pairwise search results. This would aim to
cover EMBOSS muscle & water, BLAST XML, BLAST tabular, HMMER, Bill
Pearson's FASTA alignments, and so on.

Approach  
Much of the low level parsing code to handle these file formats already
exists in Biopython, and much as the SeqIO and AlignIO modules are
linked and share code, similar links apply to the proposed SearchIO
module when using pairwise alignment file formats. However, SearchIO
will also support pairwise search results where the pairwise sequence
alignment itself is not available (e.g. the default BLAST
tabular output). A crucial aspect of this work will be to design a
pairwise-search-result object heirachy that reflects this, probably with
a subclass inheriting from both the pairwise-search-result and the
existing MultipleSequenceAlignment object. Beyond the initial challenge
of an iterator based parsing and writing framework, random access akin
to the Bio.SeqIO.index and index\_db functionality would be most
desirable for working with large datasets.

Challenges  
The project will cover a range of important file formats from major
Bioinformatics tools, thus will require familiarity with running these
tools, and understanding their output and its meaning. Inter-converting
file formats is part of this.

Difficulty and needed skills  
Medium/Hard depending on how many objectives are attempted. The student
needs to be fluent in Python and have knowledge of the
BioPython codebase. Experience with all of the command line tools listed
would be clear advantages, as would first hand experience using
BioPerl's SearchIO. You will also need to know or learn the git version
control system.

Mentors  
[Peter Cock](http://www.hutton.ac.uk/staff/peter-cock)

#### [Representation and manipulation of genomic variants](http://arklenna.tumblr.com/tagged/gsoc2012)

Student  
Lenna Peterson ([blog](http://arklenna.tumblr.com/tagged/gsoc2012))

Rationale  
Computational analysis of genomic variation requires the ability to
reliably communicate and manipulate variants. The goal of this project
is to provide facilities within BioPython to represent sequence
variation objects, convert them to and from common human and file
representations, and provide common manipulations on them.

Approach & Goals  

-   Object representation
    -   identify variation types to be represented (SNV, CNV, repeats,
        inversions, etc)
    -   develop internal machine representation for variation types
    -   ensure coverage of essential standards, including HGVS, GFF, VCF
-   External representations
    -   write parser and generators between objects and external string
        and file formats
-   Manipulations
    -   canonicalize variations with more than one valid representation
        (e.g., ins versus dup and left shifting repeats).
    -   develop coordinate mapping between genomic, cDNA, and protein
        sequences (HGVS)
-   Other
    -   release code to appropriate community efforts and write short
        manuscript
    -   implement web service for HGVS conversion

Difficulty and needed skills  
Easy-to-Medium depending on how many objectives are attempted. The
student will need have skills in most or all of: basic molecular biology
(genomes, transcripts, proteins), genomic variation, Python, BioPython,
Perl, BioPerl, NCBI Eutilities and/or Ensembl API. Experience with
computer grammars is highly desirable. You will also need to know or
learn the git version control system.

Mentors  
[Reece Hart](http://www.linkedin.com/in/reece)

[Brad Chapman](https://github.com/chapmanb)

[James Casbon](http://casbon.me/)

### 2011 (OBF)

#### [Biomolecular Interface Analysis](http://biopython.org/wiki/GSoC2011_mtrellet)

Student  
Mikael Trellet

Rationale  
Analysis of protein-protein complexes interfaces at a residue level
yields significant information on the overall binding process. Such
information can be broadly used for example in binding affinity studies,
interface design, and enzymology. To tap into it, there is a need for
tools that systematically and automatically analyze protein structures,
or that provide means to this end.
Protorop (http://www.bioinformatics.sussex.ac.uk/protorp/) is an example
of such a tool and the elevated number of citations the server has had
since its publication acknowledge its importance. However, being a
webserver, Protorop is not suited for large-scale analysis and it leaves
the community dependent on its maintainers to keep the
service available. On the other hand, Biopython’s structural biology
module, Bio.PDB, provides the ideal parsing machinery and programmatic
structures for the development of an offline, open-source library for
interface analysis. Such a library could be easily used in large-scale
analysis of protein-protein interfaces, for example in the CAPRI
experiment evaluation or in benchmark statistics. It would be also
reasonable, if time permits, to extend this module to deal with
protein-DNA or protein-RNA complexes, as Biopython supports nucleic
acids already.

Approach & Goals  

-   Add the new module backbone in current Bio.PDB code base
    -   Evaluate possible code reuse and call it into the new module
    -   Try simple calculations to be sure that there is stability
        between the different modules (parsing for example) and
        functions
-   Define a stable benchmark
    -   Select few PDB files among interface size and proteins size
        would be different
-   Extend IUPAC.Data module with residue information
    -   Deduce residues weight from Atom instead of direct dictionary
        storage
    -   Polar/charge character (dictionary or influenced by pH)
    -   Hydrophobicity scale(s)
-   Implement Extended Residue class as a subclass of Residue
-   Implement Interface object and InterfaceAnalysis module
-   Develop functions for interface analysis
    -   Calculation of interface polar character statistics (% of polar
        residues, apolar, etc)
    -   Calculation of BSA calling MSMS or HSA
    -   Calculation of SS element statistics in the interface through
        DSSP
    -   Unit tests and use of results as input for further calculations
        by other tools and scripts
-   Develop functions for Interface comparison
-   Code organization and final testing

Difficulty and needed skills  
Easy/Medium. Working knowledge of the Bio.PDB module of BioPython.
Knowledge of structural biology in general and associated file
formats (PDB).

Mentors  
[João Rodrigues](http://nmr.chem.uu.nl/~joao)

[Eric Talevich](http://etal.myweb.uga.edu/)

#### [A Python bridge for Mocapy++](http://biopython.org/wiki/GSOC2011_Mocapy)

Student  
Michele Silva

Rationale  
Discovering the structure of biomolecules is one of the biggest problems
in biology. Given an amino acid or base sequence, what is the three
dimensional structure? One approach to biomolecular structure prediction
is the construction of probabilistic models. A Bayesian network is a
probabilistic model composed of a set of variables and their joint
probability distribution, represented as a directed acyclic graph. A
dynamic Bayesian network is a Bayesian network that represents sequences
of variables. These sequences can be time-series or sequences of
symbols, such as protein sequences. Directional statistics is concerned
mainly with observations which are unit vectors in the plane or in
three-dimensional space. The sample space is typically a circle or
a sphere. There must be special directional methods which take into
account the structure of the sample spaces. The union of graphical
models and directional statistics allows the development of
probabilistic models of biomolecular structures. Through the use of
dynamic Bayesian networks with directional output it becomes possible to
construct a joint probability distribution over sequence and structure.
Biomolecular structures can be represented in a geometrically natural,
continuous space. Mocapy++ is an open source toolkit for inference and
learning using dynamic Bayesian networks that provides support for
directional statistics. Mocapy++ is excellent for constructing
probabilistic models of biomolecular structures; it has been used to
develop models of protein and RNA structure in atomic detail. Mocapy++
is used in several high-impact publications, and will form the core of
the molecular modeling package Phaistos, which will be released soon.
The goal of this project is to develop a highly useful Python interface
to Mocapy++, and to integrate that interface with the Biopython project.
Through the Bio.PDB module, Biopython provides excellent functionality
for data mining biomolecular structure databases. Integrating Mocapy++
and Biopython will allow training a probabilistic model using data
extracted from a database. Integrating Mocapy++ with Biopython will
create a powerful toolkit for researchers to quickly implement and test
new ideas, try a variety of approaches and refine their methods. It will
provide strong support for the field of biomolecular structure
prediction, design, and simulation.

Approach & Goals  
Mocapy++ is a machine learning toolkit for training and using
Bayesian networks. It has been used to develop probabilistic models of
biomolecular structures. The goal of this project is to develop a Python
interface to Mocapy++ and integrate it with Biopython. This will allow
the training of a probabilistic model using data extracted from
a database. The integration of Mocapy++ with Biopython will provide a
strong support for the field of protein structure prediction, design
and simulation.

Mentors  
[Eric Talevich](http://etal.myweb.uga.edu/)

[Thomas Hamelryck](http://wiki.binf.ku.dk/User:Thomas_Hamelryck)

#### [MocapyExt](http://biopython.org/wiki/GSOC2011_MocapyExt)

Student  
Justinas V. Daugmaudis

Rationale  
BioPython is a very popular library in Bioinformatics and
Computational Biology. Mocapy++ is a machine learning toolkit for
parameter learning and inference in dynamic Bayesian networks (DBNs),
which encode probabilistic relationships among random variables in
a domain. Mocapy++ is freely available under the GNU General Public
Licence (GPL) from SourceForge. The library supports a wide spectrum of
DBN architectures and probability distributions, including distributions
from directional statistics. Notably, Kent distribution on the sphere
and the bivariate von Mises distribution on the torus, which have proven
to be useful in formulating probabilistic models of protein and
RNA structure. Such a highly useful and powerful library, which has been
used in such projects as TorusDBN, Basilisk, FB5HMM with great success,
is the result of the long-term effort. The original Mocapy
implementation dates back to 2004, and since then the library has been
rewritten in C++. However, C++ is a statically typed and compiled
programming language, which does not facilitate rapid prototyping. As a
result, currently Mocapy++ has no provisions for dynamic loading of
custom node types, and a mechanism to plug-in new node types that would
not require to modify and recompile the library is of interest. Such a
plug-in interface would assist rapid prototyping by allowing to quickly
implement and test new probability distributions, which, in turn, could
substantially reduce development time and effort; the user would be
empowered to extend Mocapy++ without modifications and
subsequent recompilations. Recognizing this need, the project (herein
referred as MocapyEXT), with the aim to improve the current Mocapy++
node type extension mechanism, has been proposed by T. Hamelryck.

Approach & Goals  
The MocapyEXT project is largely an engineering effort to bring a
transparent Python plug-in interface to Mocapy++, where built-in and
dynamically loaded node types could be used in a uniform manner. Also,
externally implemented and dynamically loaded nodes could be modified by
a user and these changes will not necessitate the recompilation of the
client program, nor the accompanying Mocapy++ library. This will
facilitate rapid prototyping, ease the adaptation of currently existing
code, and improve the software interoperability whilst introducing
minimal changes to the existing Mocapy++ interface, thus facilitating a
smooth acceptance of the changes introduced by MocapyEXT.

Mentors  
[Eric Talevich](http://etal.myweb.uga.edu/)

[Thomas Hamelryck](http://wiki.binf.ku.dk/User:Thomas_Hamelryck)

### 2010 (OBF)

#### [Improving Bio.PDB](http://biopython.org/wiki/GSOC2010_Joao)

Student  
[João Rodrigues](http://nmr.chem.uu.nl/~joaor)

Rationale  
Biopython is a very popular library in Bioinformatics and
Computational Biology. Its Bio.PDB module, originally developed by
Thomas Hamelryck, is a simple yet powerful tool for
structural biologists. Although it provides a reliable PDB parser
feature and it allows several calculations (Neighbour Search, RMS) to be
made on macromolecules, it still lacks a number of features that are
part of a researcher's daily routine. Probing for disulphide bridges in
a structure and adding polar hydrogen atoms accordingly are two examples
that can be incorporated in Bio.PDB, given the module's clever structure
and good overall organisation. Cosmetic operations such as chain removal
and residue renaming – to account for the different existing
nomenclatures – and renumbering would also be greatly appreciated by
the community. Another aspect that can be improved for Bio.PDB is a
smooth integration/interaction layer for heavy-weights in macromolecule
simulation such as MODELLER, GROMACS, AutoDock, HADDOCK. It could be
argued that the easiest solution would be to code hooks to these
packages' functions and routines. However, projects such as the recently
developed edPDB or the more complete Biskit library render, in my
opinion, such interfacing efforts redundant. Instead, I believe it to be
more advantageous to include these software' input/output formats in
Biopython's SeqIO and AlignIO modules. This, together with the creation
of interfaces for model validation/structure checking services/software
would allow Biopython to be used as a pre- and post-simulation tool.
Eventually, it would pave the way for its inclusion in pipelines and
workflows for structure modelling, molecular dynamics, and
docking simulations.

Mentors  
[Eric Talevich](http://etal.myweb.uga.edu/)

[Peter Cock](http://www.hutton.ac.uk/staff/peter-cock)

Diana Jaunzeikare

### 2009 (NESCent)

This was the first year Biopython took part in GSoC, and we did so under
the banner of [NESCent's GSoC 2013
program](http://informatics.nescent.org/wiki/Phyloinformatics_Summer_of_Code_2009).

#### [PhyloXML](http://biopython.org/wiki/PhyloXML)

Rationale  
PhyloXML is an XML format for phylogenetic trees, designed to allow
storing information about the trees themselves (such as branch lengths
and multiple support values) along with data such as taxonomic and
genomic annotations. Connecting these pieces of evolutionary information
in a standard format is key for comparative genomics.

A Bioperl driver for phyloXML was created during the 2008 Summer of
Code; this project aims to build a similar module for the popular
Biopython package.

Mentors  
[Brad Chapman](https://github.com/chapmanb)

Christian Zmasek

#### [Biogeographical Phylogenetics for BioPython](http://biopython.org/wiki/BioGeography)

Rationale  
I developed Bio.Geography, a new module for the bioinformatics
programming toolkit Biopython. Bio.Geography expands upon Biopython's
traditional capabilities for accessing gene and protein sequences from
online databases by allowing automated searching, downloading, and
parsing of geographic location records from GBIF, the authoritative
aggregator of specimen information from natural history
collections worldwide. This will enable analyses of evolutionary
biogeography that require the areas inhabited by the species at the tips
of the phylogeny, particularly for large-scale analyses where it is
necessary to process thousands of specimen occurrence records. The
module will also facilitate applications such as species mapping, niche
modeling, error-checking of museum records, and monitoring
range changes.

Mentors  
[Brad Chapman](https://github.com/chapmanb)

Stephen Smith

David Kidd


