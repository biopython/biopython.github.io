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
-   [Eric Talevich](User%3AEricTalevich "wikilink") added support for
    [parsing and writing
    phyloXML](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009#Biopython_support_for_parsing_and_writing_phyloXML).

In 2010, another project was funded:

-   João Rodrigues [worked on the Structural Biology module
    Bio.PDB](GSOC2010_Joao "wikilink"), adding several features used in
    everyday structural bioinformatics. These features are now gradually
    being merged into the mainline with João's help.

In 2011, three projects were funded in Biopython via OBF:

-   [Mikael Trellet](User%3AMtrellet "wikilink") added [support for
    biomolecular interface analysis](GSoC2011_mtrellet "wikilink") to
    the Bio.PDB module.
-   Michele Silva wrote a [Python bridge for
    Mocapy++](GSOC2011_Mocapy "wikilink") and linked it to Bio.PDB to
    enable statistical analysis of protein structures.
-   Justinas Daugmaudis also enhanced Mocapy++ in a complementary way,
    developing a [plugin system for
    Mocapy++](GSOC2011_MocapyExt "wikilink") allowing users to easily
    write new nodes (probability distribution functions) in Python.

Please read the [GSoC page at the Open Bioinformatics
Foundation](http://www.open-bio.org/wiki/Google_Summer_of_Code) and the
main [Google Summer of Code](http://code.google.com/soc) page for more
details about the program. If you are interested in contributing as a
mentor or student next year, please introduce yourself on the [mailing
list](http://biopython.org/wiki/Mailing_lists).

2012 Project ideas
------------------

### SearchIO

Rationale : Biopython has general APIs for parsing and writing assorted sequence file formats ([SeqIO](SeqIO "wikilink")), multiple sequence alignments ([AlignIO](AlignIO "wikilink")), phylogenetic trees ([Phylo](Phylo "wikilink")) and motifs (Bio.Motif). An obvious omission is something equivalent to [BioPerl's SearchIO](bp:HOWTO:SearchIO "wikilink"). The goal of this proposal is to develop an easy-to-use Python interface in the same style as [SeqIO](SeqIO "wikilink"), [AlignIO](AlignIO "wikilink"), etc but for pairwise search results. This would aim to cover EMBOSS muscle & water, BLAST XML, BLAST tabular, HMMER, Bill Pearson's FASTA alignments, and so on.  

Much of the low level parsing code to handle these file formats already
exists in Biopython, and much as the [SeqIO](SeqIO "wikilink") and
[AlignIO](AlignIO "wikilink") modules are linked and share code, similar
links apply to the proposed SearchIO module when using pairwise
alignment file formats. However, SearchIO will also support pairwise
search results where the pairwise sequence alignment itself is not
available (e.g. the default BLAST tabular output). A crucial aspect of
this work will be to design a pairwise-search-result object heirachy
that reflects this, probably with a subclass inheriting from both the
pairwise-search-result and the existing MultipleSequenceAlignment
object.

Beyond the initial challenge of an iterator based parsing and writing
framework, random access akin to the Bio.SeqIO.index and index\_db
functionality would be most desirable for working with large datasets.

Challenges : The project will cover a range of important file formats from major Bioinformatics tools, thus will require familiarity with running these tools, and understanding their output and its meaning. Inter-converting file formats is part of this.  

<!-- -->

Involved toolkits or projects :  

:\* Biopython

Degree of difficulty and needed skills : Medium/Hard depending on how many objectives are attempted. The student needs to be fluent in Python. Experience with all of the command line tools listed would be clear advantages, as would first hand experience using [BioPerl's SearchIO](bp:HOWTO:SearchIO "wikilink"). You will also need to know or learn the git version control system.  

<!-- -->

Mentors : Peter Cock  

### Representation and manipulation of genomic variants

Rationale : Computational analysis of genomic variation requires the ability to reliably communicate and manipulate variants. The goal of this project is to provide facilities within BioPython to represent sequence variation objects, convert them to and from common human and file representations, and provide common manipulations on them.  

<!-- -->

Approach and Goals  

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

Challenges : The major challenge in this project is to design an API that separates internal representations of variation from the multiple external representations. Ideally, the libraries developed in this project will provide low-level functionality of coordinate conversion and parsing, and high-level functionality for the most common use cases. This aim requires analyzing the proposals to determine which aspects may be impossible or difficult to represent with a simple grammar.  

<!-- -->

Resources  

-   [BioPython](http://biopython.org)
-   [HGVS "nomenclature"](http://www.mutalyzer.nl/2.0/)
-   [VCF](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41),
    \[VCFtools <http://vcftools.sourceforge.net/>\]
-   [GFF3](http://www.sequenceontology.org/gff3.shtml),
    [GVF](http://www.sequenceontology.org/resources/gvf.html)
-   [GATK](http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit)

Degree of difficulty and needed skills : Easy-to-Medium depending on how many objectives are attempted. The student will need have skills in most or all of: basic molecular biology (genomes, transcripts, proteins), genomic variation, Python, BioPython, Perl, BioPerl, NCBI Eutilities and/or Ensembl API. Experience with computer grammars is highly desirable. You will also need to know or learn the git version control system.  

<!-- -->

Mentors: [Reece Hart](http://linkedin.com/in/reece) ([Locus Development](http://locusdevelopmentinc.com), San Francisco); [Brad Chapman](http://bcbio.wordpress.com); [James Casbon](http://casbon.me)  


