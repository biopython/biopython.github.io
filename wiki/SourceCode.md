---
title: SourceCode
permalink: wiki/SourceCode
layout: page
redirect_from:
 - /wiki/Git
---

### Introduction

The source code from Biopython is freely available for your use and
contribution under our [liberal
license](http://www.biopython.org/DIST/LICENSE).

The [Biopython source code](http://github.com/biopython/biopython) is
kept under a distributed version control system which allows multiple
users from around the world to work on the same code base at the same
time. We currently use
[git](http://en.wikipedia.org/wiki/Git_%28software%29) (developed by
Linus Torvalds for Linux kernel development) hosted on
[GitHub](http://github.com).

Our core developers maintain a stable trunk from which we will roll
releases as new functionality is integrated and bugs are fixed.

### Viewing the source code

You can [browse our latest source code on
github](http://github.com/biopython/biopython).

### Track changes

You can [track changes](Tracking_commits "wikilink") via
[RSS](wp:RSS_(file_format) "wikilink").

### Downloading the latest source

You can download the latest source code by clicking the Download link
near the top of the [Biopython GitHub
page](http://github.com/biopython/biopython) (this will offer you a [tar
ball](http://github.com/biopython/biopython/tarball/master) or
[zip](http://github.com/biopython/biopython/zipball/master) file).

An hourly updated copy of the code is also available at
<http://biopython.open-bio.org/SRC/biopython> (just a snapshot - no
history etc).

### Anonymous Access

Getting a copy of the repository (called "cloning" in git terminology)
is very simple using the git command line tool, you don't need an
account or password:

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
projects migrated from [CVS](CVS "wikilink") to [Subversion
(SVN)](SVN "wikilink"). Biopython had been considering [moving from CVS
to SVN](Subversion_migration "wikilink") for a while, but instead [moved
to git](GitMigration "wikilink") in September 2009. Note that
[BioRuby](http://bioruby.org) also uses github.
