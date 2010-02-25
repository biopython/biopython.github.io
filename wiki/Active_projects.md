---
title: Active projects
permalink: wiki/Active_projects
layout: wiki
---

This page provides a central location to collect references to active
projects. This is a good place to start if you are interested in
contributing to Biopython and want to find larger projects in progress.
For developers, use this to reference [git](git "wikilink") branches or
other projects which you will be working on for an extended period of
time. Please keep it up to date as projects are finished and integrated
into Biopython.

Current projects
----------------

### Population Genetics development

Giovanni and Tiago are working on expanding population genetics code in
Biopython. See the [PopGen development page](PopGen_dev "wikilink") for
more details.

### GFF parser

Brad is working on a Biopython GFF parser. Source code is available from
[git hub](http://github.com/chapmanb/bcbb/tree/master/gff).
Documentation is in progress at [GFF Parsing](GFF_Parsing "wikilink").
See blog posts on the [initial
implementation](http://bcbio.wordpress.com/2009/03/08/initial-gff-parser-for-biopython/)
and [MapReduce parallel
version](http://bcbio.wordpress.com/2009/03/22/mapreduce-implementation-of-gff-parsing-for-biopython/).

### Phylo

[Eric](User%3AEricTalevich "wikilink") is working on a new module for
phylogenetics, [Bio.Phylo](Phylo "wikilink"). It grew out of a Google
Summer of Code 2009
[project](http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798969),
mentored by Brad, to add support for
[phyloXML](http://www.phyloxml.org/) to Biopython; it also refactors
part of Bio.Nexus. Most of the code has been pushed to the main
development branch on GitHub already, but new features appear first on
Eric's [phyloxml
branch](http://github.com/etal/biopython/tree/phyloxml).

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
Bio/Geography directory of the [Geography fork of the nmatzke branch on
GitHub](http://github.com/nmatzke/biopython/tree/Geography), and you can
see a timeline and other info about ongoing development
[here](BioGeography "wikilink"). The new module is being documented on
this wiki as [BioGeography](BioGeography "wikilink").

### Roche 454 SFF parsing in Bio.SeqIO

See [Bug 2837](http://bugzilla.open-bio.org/show_bug.cgi?id=2837), based
on code from Jose Blanca. Recently merged into the trunk and should be
included in Biopython 1.54 onwards.

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

-   Build a general tool to filter sequences containing ambiguous or low
    quality bases. Chris Fields from BioPerl is interested in
    coordinating the BioPerl/Biopython implementations. See these
    threads on the mailing lists for discussion:
    <http://lists.open-bio.org/pipermail/biopython/2009-July/005355.html>,
    <http://lists.open-bio.org/pipermail/biopython/2009-July/005342.html>

<!-- -->

-   Use SQLAlchemy, an object relational mapper, for BioSQL internals.
    This would add an additional external dependency to Biopython, but
    provides ready support for additional databases like SQLite. It also
    would provide a raw object interface to BioSQL databases when the
    SeqRecord-like interface is not sufficient. Brad and Kyle have some
    initial code for this.

<!-- -->

-   Revamp the GEO SOFT parser, drawing on the ideas used in [Sean
    Davis' GEOquery parser in
    R/Bioconductor](http://www.bioconductor.org/packages/bioc/html/GEOquery.html).
    See also [this page](http://www.warwick.ac.uk/go/peter_cock/r/geo/).

Enhancement list
----------------

Maintaining software involves incremental improvements for new format
changes and removal of bugs. Please see our
[bugzilla](http://bugzilla.open-bio.org/) page for a current list. Post
to the developer mailing list if you are interested in tackling any open
issues.
