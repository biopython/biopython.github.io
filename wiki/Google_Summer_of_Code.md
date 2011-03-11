---
title: Google Summer of Code
permalink: wiki/Google_Summer_of_Code
layout: wiki
---

As part of the Open Bioinformatics Foundation, Biopython is
participating in Google Summer of Code (GSoC) again in 2010. We are
supporting João Rodrigues in his project, "[Extending Bio.PDB:
broadening the usefulness of BioPython's Structural Biology
module](GSOC2010_Joao "wikilink")."

In 2009, Biopython was involved with GSoC in collaboration with our
friends at
[NESCent](https://www.nescent.org/wg_phyloinformatics/Main_Page), and
had two projects funded:

-   Nick Matzke worked on [Biogeographical
    Phylogenetics](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009#Biogeographical_Phylogenetics_for_BioPython).
-   Eric Talevich added support for [parsing and writing
    phyloXML](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009#Biopython_support_for_parsing_and_writing_phyloXML).

In 2010, another project was funded:

-   João Rodrigues worked on [the Structural Biology module
    Bio.PDB](http://www.biopython.org/wiki/GSOC2010_Joao) adding several
    features used in everyday structural bioinformatics.

Please read the [GSoC page at the Open Bioinformatics
Foundation](http://www.open-bio.org/wiki/Google_Summer_of_Code) and the
main [Google Summer of Code](http://code.google.com/soc) page for more
details about the program. If you are interested in contributing as a
mentor or student next year, please introduce yourself on the [mailing
list](http://biopython.org/wiki/Mailing_lists).

2011 Project ideas
------------------

### Biopython and PyCogent interoperability

Rationale : [PyCogent](http://pycogent.sourceforge.net/) and [Biopython](http://biopython.org/wiki/Main_Page) are two widely used toolkits for performing computational biology and bioinformatics work in Python. The libraries have had traditionally different focuses: with Biopython focusing on sequence parsing and retrieval and PyCogent on evolutionary and phylogenetic processing. Both user communities would benefit from increased interoperability between the code bases, easing the developing of complex workflows.  

<!-- -->

Approach : The student would focus on soliciting use case scenarios from developers and the larger communities associated with both projects, and use these as the basis for adding glue code and documentation to both libraries. Some use cases of immediate interest as a starting point are:  

:\* Allow round-trip conversion between biopython and pycogent core
objects (sequence, alignment, tree, etc.).

:\* Building workflows using Codon Usage analyses in PyCogent with
clustering code in Biopython.

:\* Connecting Biopython acquired sequences to PyCogent's alignment,
phylogenetic tree preparation and tree visualization code.

:\* Integrate Biopython's [phyloXML
support](http://biopython.org/wiki/Phylo), developed during GSoC 2009,
with PyCogent.

:\* Develop a standardised controller architecture for interrogation of
genome databases by extending PyCogent's Ensembl code, including export
to Biopython objects.

Challenges : This project provides the student with a lot of freedom to create useful interoperability between two feature rich libraries. As opposed to projects which might require churning out more lines of code, the major challenge here will be defining useful APIs and interfaces for existing code. High level inventiveness and coding skill will be required for generating glue code; we feel library integration is an extremely beneficial skill. We also value clear use case based documentation to support the new interfaces.  

<!-- -->

Involved toolkits or projects :  

:\* [Biopython](http://biopython.org/wiki/Main_Page)

:\* [PyCogent](http://pycogent.sourceforge.net/)

Degree of difficulty and needed skills : Medium to Hard. At a minimum, the student will need to be highly competent in Python and become familiar with core objects in PyCogent and Biopython. Sub-projects will require additional expertise, for instance: familiarity with concepts in phylogenetics and genome biology; understanding SQL dialects.  

<!-- -->

Mentors : [Gavin Huttley](http://jcsmr.anu.edu.au/org/dmb/compgen/), [Rob Knight](http://chem.colorado.edu/index.php?option=com_content&view=article&id=263:rob-knight), [Brad Chapman](http://bcbio.wordpress.com), [Eric Talevich](User%3AEricTalevich "wikilink")  

### Accessing R phylogenetic tools from Python

Rationale : The [R statistical language](http://www.r-project.org/) is a powerful open-source environment for statistical computation and visualization. [Python](http://www.python.org/) serves as an excellent complement to R since it has a wide variety of available libraries to make data processing, analysis, and web presentation easier. The two can be smoothly interfaced using [Rpy2](http://bitbucket.org/lgautier/rpy2/), allowing programmers to leverage the best features of each language. Here we propose to build Rpy2 library components to help ease access to phylogenetic and biogeographical libraries in R.  

<!-- -->

Approach : Rpy2 contains higher level interfaces to popular R libraries. For instance, the [ggplot2 interface](http://rpy.sourceforge.net/rpy2/doc-2.1/html/graphics.html#package-ggplot2) allows python users to access powerful plotting functionality in R with an intuitive API. Providing similar high level APIs for biological toolkits available in R would help expose these toolkits to a wider audience of Python programmers. A nice introduction to phylogenetic analysis in R is available from Rich Glor at the [Bodega Bay Marine Lab wiki](http://bodegaphylo.wikispot.org/Phylogenetics_and_Comparative_Methods_in_R). Some examples of R libraries for which integration would be welcomed are:  

:\* [ape (Analysis of Phylogenetics and
Evolution)](http://ape.mpl.ird.fr/) -- an interactive library
environment for phylogenetic and evolutionary analyses

:\* [ade4](http://pbil.univ-lyon1.fr/ADE-4/home.php?lang=eng) -- Data
Analysis functions to analyse Ecological and Environmental data in the
framework of Euclidean Exploratory methods

:\* [geiger](http://cran.r-project.org/web/packages/geiger/index.html)
-- Running macroevolutionary simulation, and estimating parameters
related to diversification from comparative phylogenetic data.

:\* [picante](http://picante.r-forge.r-project.org/) -- R tools for
integrating phylogenies and ecology

:\* [mefa](http://mefa.r-forge.r-project.org/) -- multivariate data
handling for ecological and biogeographical data

Challenges : The student would have the opportunity to learn an available R toolkit, and then code in Python and R to make this available via an intuitive API. This will involve digging into the R code examples to discover the most useful parts for analysis, and then projecting this into a library that is intuitive to Python coders. Beyond the coding and design aspects, the student should feel comfortable writing up use case documentation to support the API and encourage its adoption.  

<!-- -->

Involved toolkits or projects :  

:\* [ape (Analysis of Phylogenetics and
Evolution)](http://ape.mpl.ird.fr/)

:\* [Rpy2](http://bitbucket.org/lgautier/rpy2/)

:\* [Biopython](http://biopython.org/wiki/Main_Page)

Degree of difficulty and needed skills : Moderate. The project requires familiarity with coding in Python and R, and knowledge of phylogeny or biogeography. The student has plenty of flexibility to define the project based on their biological interests (e.g. [microarrays and heatmaps](http://www.warwick.ac.uk/go/peter_cock/python/heatmap/)); there is also the possibility to venture far into data visualization once access to analysis methods is made. [GenGIS](http://kiwi.cs.dal.ca/GenGIS/Main_Page) and can give ideas about what is possible.  

<!-- -->

Mentors : [Laurent Gautier](http://dk.linkedin.com/pub/laurent-gautier/8/81/869), [Brad Chapman](http://bcbio.wordpress.com), [Peter Cock](http://www.scri.ac.uk/staff/petercock)  

### Mocapy++Biopython: from data to probabilistic models of biomolecules

Rationale : [Mocapy++](http://sourceforge.net/projects/mocapy/) is a machine learning toolkit for training and using [Bayesian networks](http://en.wikipedia.org/wiki/Bayesian_network). Mocapy++ supports the use of [directional statistics](http://en.wikipedia.org/wiki/Directional_statistics); the statistics of angles, orientations and directions. This unique feature of Mocapy++ makes the toolkit especially suited for the formulation of probabilistic models of biomolecular structure. The toolkit has already been used to develop (published and peer reviewed) models of [protein](http://www.pnas.org/content/105/26/8932.abstract?etoc) and [RNA](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000406) structure in atomic detail. Mocapy++ is implemented in C++, and does not provide any Python bindings. The goal of this proposal is to develop an easy-to-use Python interface to Mocapy++, and to integrate this interface with the Biopython project. Through its [Bio.PDB](http://biopython.org/DIST/docs/cookbook/biopdb_faq.pdf) module (initially implemented by the mentor of this proposal, [T. Hamelryck](http://www.binf.ku.dk/research/structural_bioinformatics/)), Biopython provides excellent functionality for data mining of biomolecular structure databases. Integrating Mocapy++ and Biopython would create strong synergy, as it would become quite easy to extract data from the databases, and subsequently use this data to train a probabilistic model. As such, it would provide a strong impulse to the field of protein structure prediction, design and simulation. Possible applications beyond bioinformatics are obvious, and include probabilistic models of human or animal movement, or any other application that involves directional data.  

<!-- -->

Approach : Ideally, the student would first gain some understanding of the theoretical background of the algorithms that are used in Mocapy++, such as parameter learning of Bayesian networks using [Stochastic Expectation Maximization (S-EM)](http://en.wikipedia.org/wiki/Expectation-maximization_algorithm). Next, the student would study some of the use cases of the toolkit, making use of some of the published articles that involve Mocapy++. After becoming familiar with the internals of Mocapy++, Python bindings will then be implemented using the [Boost C++ library](http://www.boost.org). Based on the use cases, the student would finally implement some example applications that involve data mining of biomolecular structure using Biopython, the subsequent formulation of probabilistic models using Python-Mocapy++, and its application to some biologically relevant problem. Schematically, the following steps are involved for the student:  

:\* Gaining some understanding of S-EM and directional statistics

:\* Study of Mocapy++ use cases

:\* Study of Mocapy++ internals and code

:\* Design of interface strategy

:\* Implementing Python bindings using Boost

:\* Example applications, involving Bio.PDB data mining

Challenges : The project is highly interdisciplinary, and ideally requires skills in programming (C++, Python, wrapping C++ libraries in Python, Boost), machine learning, knowledge of biomolecular structure and statistics. The project could be extended (for example, by implementing additional functionality in Mocapy++) or limited (for example by limiting the time spent on understanding the theory behind Mocapy++).  

<!-- -->

Involved toolkits or projects :  

:\* [Biopython](http://biopython.org/wiki/Main_Page)

:\* [Mocapy++](http://sourceforge.net/projects/mocapy/)

Degree of difficulty and needed skills : Hard. The student needs to be fluent in C++, Python and the [C++ Boost library](http://www.boost.org). Experience with machine learning, Bayesian statistics and biomolecular structure would be clear advantages.  

<!-- -->

Mentors : [Thomas Hamelryck](http://www.binf.ku.dk/research/structural_bioinformatics/)  


