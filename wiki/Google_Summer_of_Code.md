---
title: Google Summer of Code
permalink: wiki/Google_Summer_of_Code
layout: wiki
---

As part of the Open Bioinformatics Foundation, Biopython is
participating in Google Summer of Code (GSoC) again in 2011. This page
contains a list of project ideas for the upcoming summer; potential GSoC
students can base an application on any of these ideas, or propose
something new.

In 2009, Biopython was involved with GSoC in collaboration with our
friends at
[NESCent](https://www.nescent.org/wg_phyloinformatics/Main_Page), and
had two projects funded:

-   Nick Matzke worked on [Biogeographical
    Phylogenetics](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009#Biogeographical_Phylogenetics_for_BioPython).
-   Eric Talevich added support for [parsing and writing
    phyloXML](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009#Biopython_support_for_parsing_and_writing_phyloXML).

In 2010, another project was funded:

-   João Rodrigues [worked on the Structural Biology module
    Bio.PDB](GSOC2010_Joao "wikilink"), adding several features used in
    everyday structural bioinformatics. These features are now gradually
    being merged into the mainline with João's help.

Please read the [GSoC page at the Open Bioinformatics
Foundation](http://www.open-bio.org/wiki/Google_Summer_of_Code) and the
main [Google Summer of Code](http://code.google.com/soc) page for more
details about the program. If you are interested in contributing as a
mentor or student next year, please introduce yourself on the [mailing
list](http://biopython.org/wiki/Mailing_lists).

2011 Project ideas
------------------

### Mocapy++Biopython: from data to probabilistic models of biomolecules

Rationale : [Mocapy++](http://sourceforge.net/projects/mocapy/) is a machine learning toolkit for training and using [Bayesian networks](http://en.wikipedia.org/wiki/Bayesian_network). Mocapy++ supports the use of [directional statistics](http://en.wikipedia.org/wiki/Directional_statistics); the statistics of angles, orientations and directions. This unique feature of Mocapy++ makes the toolkit especially suited for the formulation of probabilistic models of biomolecular structure. The toolkit has already been used to develop (published and peer reviewed) models of [protein](http://www.pnas.org/content/105/26/8932.abstract?etoc) and [RNA](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000406) structure in atomic detail. Mocapy++ is implemented in C++, and does not provide any Python bindings. The goal of this proposal is to develop an easy-to-use Python interface to Mocapy++, and to integrate this interface with the Biopython project. Through its [Bio.PDB](http://biopython.org/DIST/docs/cookbook/biopdb_faq.pdf) module (initially implemented by the mentor of this proposal, [T. Hamelryck](http://www.binf.ku.dk/research/structural_bioinformatics/)), Biopython provides excellent functionality for data mining of biomolecular structure databases. Integrating Mocapy++ and Biopython would create strong synergy, as it would become quite easy to extract data from the databases, and subsequently use this data to train a probabilistic model. As such, it would provide a strong impulse to the field of protein structure prediction, design and simulation. Possible applications beyond bioinformatics are obvious, and include probabilistic models of human or animal movement, or any other application that involves directional data.  

<!-- -->

Approach : Ideally, the student (or several students) would first gain some understanding of the theoretical background of the algorithms that are used in Mocapy++, such as parameter learning of Bayesian networks using [Stochastic Expectation Maximization (S-EM)](http://en.wikipedia.org/wiki/Expectation-maximization_algorithm). Next, the student would study some of the use cases of the toolkit, making use of some of the published articles that involve Mocapy++. After becoming familiar with the internals of Mocapy++, Python bindings will then be implemented using the [Boost C++ library](http://www.boost.org). Based on the use cases, the student would finally implement some example applications that involve data mining of biomolecular structure using Biopython, the subsequent formulation of probabilistic models using Python-Mocapy++, and its application to some biologically relevant problem. Schematically, the following steps are involved for the student:  

:\* Gaining some understanding of S-EM and directional statistics

:\* Study of Mocapy++ use cases

:\* Study of Mocapy++ internals and code

:\* Design of interface strategy

:\* Implementing Python bindings using Boost

:\* Example applications, involving Bio.PDB data mining

Challenges : The project is highly interdisciplinary, and ideally requires skills in programming (C++, Python, wrapping C++ libraries in Python, Boost), machine learning, knowledge of biomolecular structure and statistics. The project could be extended (for example, by implementing additional functionality in Mocapy++) or limited (for example, by limiting the time spent on understanding the theory behind Mocapy++). The project would certainly benefit from several students with complementary skills.  

<!-- -->

Involved toolkits or projects :  

:\* [Biopython](http://biopython.org/wiki/Main_Page)

:\* [Mocapy++](http://sourceforge.net/projects/mocapy/)

Degree of difficulty and needed skills : Hard. The student needs to be fluent in C++, Python and the [C++ Boost library](http://www.boost.org). Experience with machine learning, Bayesian statistics and biomolecular structure would be clear advantages.  

<!-- -->

Mentors : [Thomas Hamelryck](http://www.binf.ku.dk/research/structural_bioinformatics/)  

### Variant representation, parser, generator, and coordinate converter

Rationale : Computational analysis of genomic variation requires the ability to reliably translate between human and computer representations of genomic variants. While several standards for human variation syntax have been proposed, community support is limited because of the technical complexity of the proposals and the lack of software libraries that implement them. The goal of this project is to initiate freely-available, language-neutral tools to parse, generate, and convert between representations of genomic variation.  

<!-- -->

Approach and Goals :  

:\* identify variation types to be represented (SNV, CNV, repeats,
inversions, etc)

:\* develop internal machine representation for variation types in
Python

:\* develop language-neutral grammar for the (reasonably) supportable
subset of the Human Genome Variation Society nomeclature guidelines

:\* write a Python library to convert between machine and human
representations of variation (i.e., parsing and generating)

:\* develop coordinate mapping between genomic, cDNA, and protein
sequences (at least)

:\* release code to appropriate community efforts and write short
manuscript

:\* as time permits:

:\*\* build Perl modules or Java libraries with identical functionality

:\*\* develop syntactic and semantic validation

:\*\* implement web service for coordinate conversion using NCBI
Eutilities

:\*\* develop a new variant syntax that is representation-complete

Challenges : The major challenge in this project is to design an API which cleanly separates internal representations of variation from the multiple external representations. For example, coordinate conversion per se does not require any sequence information, but validating a variant does. Ideally, the libraries developed in this project will provide low-level functionality of coordinate conversion and parsing, and high-level functionality for the most common use cases. This aim requires analyzing the proposals to determine which aspects may be impossible or difficult to represent with a simple grammar.  

<!-- -->

Involved toolkits or projects :  

:\* BioPython

:\* Related: <http://www.mutalyzer.nl/2.0/>,
<http://www.hgvs.org/mutnomen/>

Degree of difficulty and needed skills : Easy-to-Medium depending on how many objectives are attempted. The student will need have skills in most or all of: basic molecular biology (genomes, transcripts, proteins), genomic variation, Python, BioPython, Perl, BioPerl, NCBI Eutilities and/or Ensembl API. Experience with computer grammars is highly desirable.  

<!-- -->

Mentors : Reece Hart (Locus Development, San Francisco); [Brad Chapman](http://bcbio.wordpress.com)  


