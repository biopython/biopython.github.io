---
title: Git
permalink: wiki/Git
layout: wiki
---

### About git

The Biopython source code is kept under a version control system which
allows multiple users from around the world to work on the same code
base at the same time. We currently use
[git](http://en.wikipedia.org/wiki/Git_%28software%29) (developed by
Linus Torvalds for Linux kernel development). Prior to September 2009,
we used to use [CVS](CVS "wikilink").

An hourly updated CVS checkout of Biopython was available at
<http://biopython.open-bio.org/SRC/biopython> and we intend to have this
updated from git shortly.

You can also [track
changes](http://github.com/feeds/biopython/commits/biopython/master) via
[RSS](wp:RSS_(file_format) "wikilink").

### Viewing git sources

Our current development git sources are available for viewing via
[GitHub](http://github.com/):

[Biopython GitHub Home](http://github.com/biopython/biopython)

This is main repository. You can download these sources by clicking the
Download link near the top of the page.

### Anonymous git Access

You can also access (clone) a copy using the git command line tool. See
[GitUsage](GitUsage "wikilink").

### Write git Access

In order to make changes to the official repository, you will need a
github account with collaborator status. Write CVS access is available
for Biopython developers. This is normally given on a case by case
basis, and the best place to discuss getting write access is on the
[Biopython Development mailing
list](mailto:biopython-dev@biopython.org).

Write access will of course be given to existing Biopython developers
who had CVS commit rights.

Once you have access, see the instructions on
[GitUsage](GitUsage "wikilink")

### Migration from CVS

Most of the other [Open Bioinformatics Foundation](http://open-bio.org)
projects migrated from CVS to [Subversion (SVN)](SVN "wikilink").
Biopython had been considering [moving from CVS to
SVN](Subversion_migration "wikilink") for a while, but instead [moved to
git](GitMigration "wikilink"). [BioRuby](http://bioruby.org) also uses
github.
