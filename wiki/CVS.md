---
title: CVS
permalink: wiki/CVS
layout: wiki
---

### About CVS

The Biopython source code is kept under a version control system which
allows multiple users from around the world to work on the same code
base at the same time. We currently use CVS (Concurrent Versions System)
for this purpose.

An hourly updated CVS checkout of Biopython is available at
<http://biopython.open-bio.org/SRC/biopython>

You can also [track CVS changes](Tracking_CVS_commits "wikilink") via
[RSS](wp:RSS_(file_format) "wikilink").

### Viewing CVS Sources

Our current development CVS sources are available for viewing via
[ViewCVS](http://viewcvs.sourceforge.net/):

[Biopython CVS
Home](http://cvs.biopython.org/cgi-bin/viewcvs/viewcvs.cgi/?cvsroot=biopython)

The main repository is the biopython directory. You can download these
sources by clicking the Download tarball link at the bottom of the page.

### Anonymous CVS Access

You can also access a read-only version of the CVS sources using our
anonymous CVS server. Directions are available at:

[<http://cvs.biopython.org/>](http://cvs.biopython.org/)

To summarize, first login (password is 'cvs'), then checkout the current
code:

`cvs -d :pserver:cvs@code.open-bio.org:/home/repository/biopython login`  
`cvs -d :pserver:cvs@code.open-bio.org:/home/repository/biopython checkout biopython`

Once this is done, you can at a later date update your local copy with
one line from within the biopython directory:

`cvs update`

### Write CVS Access

Write CVS access is available for Biopython developers. This is normally
given on a case by case basis, and the best place to discuss getting
write access is on the [Biopython Development mailing
list](mailto:biopython-dev@biopython.org). Once you have access, the
instructions on [BioPerl's CVS wiki
page](http://www.bioperl.org/wiki/Using_CVS) are very helpful. To
summarize, checkout the code with this command, substituting your actual
username:

`cvs -d:ext:USERNAME@code.open-bio.org:/home/repository/biopython checkout biopython`

### Migration to Subversion

Most of the other [Open Bioinformatics Foundation](http://open-bio.org)
projects have already migrated from CVS to [Subversion
(SVN)](SVN "wikilink"). Biopython has been considering [moving from CVS
to SVN](Subversion_migration "wikilink") for a while, but we are
currently [evaluating moving to Git instead](GitMigration "wikilink").
