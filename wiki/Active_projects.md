---
title: Active projects
permalink: wiki/Active_projects
layout: wiki
---

This page provides a central location to collect references to active
projects. This is a good place to start if you are interested in
contributing to Biopython and want to find larger projects in progress.
For developers, use this to reference git branches or other projects
which you will be working on for an extended period of time. Please keep
it up to date as projects are finished and integrated into Biopython.

Current projects
----------------

### Population Genetics development

Giovanni and Tiago are working on expanding population genetics code in
Biopython. See the [PopGen development page](PopGen_dev "wikilink") for
more details.

### GFF parser

Brad is working on a Biopython GFF parser. Source code is available from
[git hub](http://github.com/chapmanb/bcbb/tree/master/gff). See blog
posts on the [initial
implementation](http://bcbio.wordpress.com/2009/03/08/initial-gff-parser-for-biopython/)
and [MapReduce parallel
version](http://bcbio.wordpress.com/2009/03/22/mapreduce-implementation-of-gff-parsing-for-biopython/).

### PhyloXML driver (GSoC)

Eric is working on supporting the [PhyloXML](http://www.phyloxml.org/)
format, as a
[project](http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798969)
for Google Summer of Code 2009. Brad is mentoring this project. The code
lives on a branch in
[GitHub](http://github.com/etal/biopython/tree/phyloxml), and you can
see a timeline and other info about ongoing development
[here](http://github.com/etal/biopython/tree/phyloxml/Bio/PhyloXML/).
The new module is being documented on this wiki as
[PhyloXML](PhyloXML "wikilink").

### Biogeography (GSoC)

[Nick](Matzke "wikilink") is working on developing a Biogeography module
for BioPython. This work is funded by [Google Summer of Code
2009](http://socghop.appspot.com/program/home/google/gsoc2009) through
NESCENT's [Phyloinformatics Summer of Code
2009](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009).
See the project proposal at: [Biogeographical Phylogenetics for
BioPython](http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798250).
The mentors are [Stephen Smith](http://blackrim.org/) (primary), [Brad
Chapman](http://bcbio.wordpress.com/), and [David
Kidd](http://evoviz.nescent.org/). The code currently lives at the
nmatzke branch on
[GitHub](http://github.com/nmatzke/biopython/tree/master), and you can
see a timeline and other info about ongoing development
[here](http://github.com/nmatzke/biopython/tree/master). The new module
is being documented on this wiki as
[BioGeography](BioGeography "wikilink").

### Open Enhancement Bugs

This [Bugzilla
Search](http://bugzilla.open-bio.org/buglist.cgi?product=Biopython&bug_status=NEW&bug_status=ASSIGNED&bug_status=REOPENED&bug_severity=enhancement)
will list all open enhancement bugs (any filed by core developers are
fairly likely to be integrated, some are just wish list entries).

Project ideas
-------------

Please add any ideas or proposals for new additions to Biopython. Bugs
and enhancements for current code should be discussed though our
bugzilla interface.

-   Use SQLAlchemy, an object relational mapper, for BioSQL internals.
    This would add an additional external dependency to Biopython, but
    provides ready support for additional databases like SQLite. It also
    would provide a raw object interface to BioSQL databases when the
    SeqRecord-like interface is not sufficient. Brad has some initial
    code for this.

<!-- -->

-   Revamp the GEO SOFT parser, drawing on the ideas used in [Sean
    Davis' GEOquery parser in
    R/Bioconductor](http://www.bioconductor.org/packages/bioc/html/GEOquery.html).
    See also [this page](http://www.warwick.ac.uk/go/peter_cock/r/geo/).

<!-- -->

-   Roche 454 SFF support in Bio.SeqIO is being discussed on the
    mailing lists.

Enhancement list
----------------

Maintaining software involves incremental improvements for new format
changes and removal of bugs. Please see our
[bugzilla](http://bugzilla.open-bio.org/) page for a current list. Post
to the developer mailing list if you are interested in tackling any open
issues.
