---
Title: Source Code
permalink: wiki/SourceCode
layout: page
redirect_from:
 - wiki/CVS
 - wiki/Git
 - wiki/SVN
 - wiki/Subversion
 - wiki/Tracking_commits
 - wiki/Tracking_CVS_commits
---

### Introduction

Biopython is currently released under the liberal "Biopython License
Agreement", but as part of a plan to switch to the more commonly used
"3-Clause BSD License", some of the code is explicitly dual licensed
under your choice of these two options. For details, see our [license
file](https://github.com/biopython/biopython/blob/master/LICENSE.rst).

The [Biopython source code](http://github.com/biopython/biopython) is
kept under a distributed version control system which allows multiple
users from around the world to work on the same code base at the same
time. We currently use
[git](http://en.wikipedia.org/wiki/Git_%28software%29) (developed by
Linus Torvalds for Linux kernel development) hosted on
[GitHub](http://github.com).

Our core developers maintain a stable trunk from which we will roll
releases as new functionality is integrated and bugs are fixed.

### Viewing the Source Code

You can [browse our latest source code on
github](http://github.com/biopython/biopython).

### Track Changes

You can track code development by [RSS feed](https://github.com/biopython/biopython/commits/master.atom)
or the [Biopython mailing list](mailto:biopython@biopython.org).
See also our other [mailing lists](Mailing_lists "wikilink").

### Downloading the Latest Source

You can download the latest source code by clicking the Download link
near the top of the [Biopython GitHub
page](http://github.com/biopython/biopython) (this will offer you a [tar
ball](http://github.com/biopython/biopython/tarball/master) or
[zip](http://github.com/biopython/biopython/zipball/master) file).

### Anonymous Access

Getting a copy of the repository (called "cloning" in git terminology)
is very simple using the git command line tool, you don't need an
account or password:

``` bash
git clone git://github.com/biopython/biopython.git
```

This command creates a local copy of the entire Biopython repository on
your machine (your own personal copy of the official repository with its
complete history). You can update this local copy at the command line
(from within the Biopython repository directory) with:

``` bash
git pull origin
```

You can even make *local* changes and commit them to this local copy,
see [GitUsage](GitUsage "wikilink") or the git documentation for further
information.

### Write Access

Most changes are submitted as pull requests via GitHub.

In order to directly make changes to the official repository, you will need a
GitHub account with collaborator status. Write access is available for
Biopython developers (including all those who previously had CVS commit
rights).

This is normally given on a case by case basis, and the best place to
discuss getting write access is on the [Biopython mailing
list](Mailing_lists).

Once you have access, see the instructions on
[GitUsage](GitUsage "wikilink")

### Migration from CVS

Most of the other [Open Bioinformatics Foundation](http://open-bio.org)
projects migrated from CVS to SVN (Subversion), and later to git hosted
at GitHub.

While Biopython did considering moving from CVS to SVN, instead we
migrated directly from CVS to git in September 2009.
