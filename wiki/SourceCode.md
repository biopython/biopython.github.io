---
title: SourceCode
permalink: wiki/SourceCode
layout: wiki
---

The source code from Biopython is freely available for your use and
contribution. The core developers maintain a stable trunk under revision
control from which we will roll releases as new functionality is
integrated and bugs are fixed.

Until September 2009, we used [CVS](CVS "wikilink") as our revision
control system, but have now [migrated to git](GitMigration "wikilink").

### About git

The Biopython source code is kept under a distributed version control
system which allows multiple users from around the world to work on the
same code base at the same time. We currently use
[git](http://en.wikipedia.org/wiki/Git_%28software%29) (developed by
Linus Torvalds for Linux kernel development).

An hourly updated CVS checkout of Biopython was available at
<http://biopython.open-bio.org/SRC/biopython> and we intend to have this
updated from git shortly.

You can [track
changes](http://github.com/feeds/biopython/commits/biopython/master) via
[RSS](wp:RSS_(file_format) "wikilink").

### Viewing git sources

Our current development git sources are available for viewing via
[GitHub](http://github.com/):

[Biopython GitHub Home](http://github.com/biopython/biopython)

This is the main repository. You can download these sources by clicking
the Download link near the top of the page.

### Anonymous Access

Getting a copy of the repository (called "cloning" in Git terminology)
without a GitHub account is very simple using the git command line tool:

`git clone `[`git://github.com/biopython/biopython.git`](git://github.com/biopython/biopython.git)

This command creates a local copy of the entire Biopython repository on
your machine (your own personal copy of the official repository with its
complete history). You can update this local copy at the command line
(from within the Biopython repository directory) with:

`git pull origin`

You can even make *local* changes and commit them to this local copy,
see [GitUsage](GitUsage "wikilink") or the git documentation for further
information.

### Write Access

In order to make changes to the official repository, you will need a
github account with collaborator status. Write access is available for
Biopython developers (including all those who previously had CVS commit
rights).

This is normally given on a case by case basis, and the best place to
discuss getting write access is on the [Biopython Development mailing
list](mailto:biopython-dev@biopython.org).

Once you have access, see the instructions on
[GitUsage](GitUsage "wikilink")

### Migration from CVS

Most of the other [Open Bioinformatics Foundation](http://open-bio.org)
projects migrated from CVS to [Subversion (SVN)](SVN "wikilink").
Biopython had been considering [moving from CVS to
SVN](Subversion_migration "wikilink") for a while, but instead [moved to
git](GitMigration "wikilink"). [BioRuby](http://bioruby.org) also uses
github.
