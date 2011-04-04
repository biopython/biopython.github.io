---
title: Download
permalink: wiki/Download
layout: wiki
---

Current Release - 1.57 - 2 April 2011
=====================================

See also [What's
new](https://github.com/biopython/biopython/raw/master/NEWS).

### Files

-   [biopython-1.57.tar.gz](http://biopython.org/DIST/biopython-1.57.tar.gz)
    6,783 Kb -- Source Tarball
-   [biopython-1.57.zip](http://biopython.org/DIST/biopython-1.57.zip)
    7,446 Kb -- Source Zip File
-   [biopython-1.57.win32-py2.4-unsupported.exe](http://biopython.org/DIST/biopython-1.57.win32-py2.4-unsupported.exe)
    1,405 Kb -- 32 bit Windows Installer for Python 2.4 (which we no
    longer officially support) and NumPy 1.1
-   [biopython-1.57.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.57.win32-py2.5.exe)
    1,405 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.57.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.57.win32-py2.6.exe)
    1,428 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.57.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.57.win32-py2.7.exe)
    1,555 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5

Please report any issues on our [mailing
lists](mailing_lists "wikilink") or [bug
tracker](http://redmine.open-bio.org/projects/biopython).

Note we don't (yet) have official 64 bit Windows Installers - however,
Christoph Gohlke has kindly made [Windows 64bit
installers](http://www.lfd.uci.edu/~gohlke/pythonlibs/) for NumPy and
Biopython (and other tools) available for testing.

### Installation Instructions

For Windows we provide click-and-run installers (specific to your
version of python), but you will first need to install some prerequisite
software (listed below, in particular, NumPy).

Most Linux distributions will include an optional Biopython package
(described below), and will take care of any prerequisite software
automatically.

For Mac OS X, we recommend installing from source (see below). You will
need to have installed Apple's XCode tools *including* the optional 10.4
SDK (check the option for 10.4 support when installing Xcode tools).

Otherwise you typically install from source by downloading and
uncompressing the archive, then running the commands:

`python setup.py build`  
`python setup.py test`  
`sudo python setup.py install`

If you have trouble, see the full installation instructions:

-   [HTML Full Installation
    Instructions](http://biopython.org/DIST/docs/install/Installation.html)
-   [PDF Full Installation
    Instructions](http://biopython.org/DIST/docs/install/Installation.pdf)

If you are still stuck, sign up to the [Biopython mailing
list](Mailing_lists "wikilink") and ask for help there.

### Required Software

-   [Python 2.4, 2.5, 2.6 or 2.7](http://www.python.org) (Biopython 1.50
    was the last release of Biopython to support Python 2.3)
-   [C compiler (if compiling
    from source)](http://www.python.org/doc/current/inst/inst.html) You
    need a C compiler supported by distutils, gcc will work fine on
    UNIX-like platforms. This is not needed on Windows if using the
    install programs provided above. On Mac OS, we recommend you install
    Apple's XCode *including* the 10.4 SDK.
-   [NumPy (Numerical Python)](http://numpy.scipy.org/). Note that until
    BioPython 1.49, Biopython used the older Numeric library. We have
    tested NumPy 1.0, 1.1, 1.2, and 1.3 with Biopython.

### Optional Software

For compiling Biopython:

-   [flex: The Fast Lexical Analyzer](http://flex.sourceforge.net/) --
    for building Bio.PDB.mmCIF.MMCIFlex which is used to parse
    macromolecular Crystallographic Information Files (mmCIF)

Some parts of Biopython use the following additional python libraries:

-   [ReportLab](http://www.reportlab.org/downloads.html) -- used for pdf
    graphics code
-   [MySQLdb](http://sourceforge.net/projects/mysql-python) -- used for
    [BioSQL](BioSQL "wikilink") with a MySQL database

In addition Biopython includes wrapper code for calling a number of
third party command line tools including:

-   [Wise2](http://www.ebi.ac.uk/Wise2/) -- for command line tool dnal
-   [NCBI Standalone
    BLAST](http://www.ncbi.nlm.nih.gov/blast/download.shtml) -- command
    line tool for running BLAST on your local machine
-   [Clustalw](ftp://ftp.ebi.ac.uk/pub/software/unix/clustalw/) --
    command line tool for building sequence alignments
-   [SIMCOAL2](http://cmpg.unibe.ch/software/simcoal2/) and
    [FDist](http://www.rubic.rdg.ac.uk/~mab/software.html) -- command
    line tools for population genetics
-   [EMBOSS](http://emboss.sourceforge.net/) -- lots of useful command
    line tools.

Easy Install
============

We don't officially sanction this option, but it has been reported to
work fine. If you have
[easy\_install](http://peak.telecommunity.com/DevCenter/EasyInstall)
installed on your computer, you can download and install the latest
Biopython distribution by simply executing this command:

`easy_install -f `[`http://biopython.org/DIST/`](http://biopython.org/DIST/)` biopython`

You will have to have administrator's rights to do this. On a Unix style
system this is normally done by:

`sudo easy_install -f `[`http://biopython.org/DIST/`](http://biopython.org/DIST/)` biopython`

Packages
========

For those of you using Linux, the easiest way to install Biopython is
through your distribution's package management system. However, unless
you are running a recent release of your Linux Distribution, you may
find that the Biopython packages available to be a little out of date.
You might want to see if there is a backport available, otherwise you
will have to install Biopython from source.

### Ubuntu or Debian

You should be able to install Biopython and its dependencies using the
Synaptic GUI tool (on the main menu under System / Administration /
Synaptic Package Manager), or at the command line using:

`sudo apt-get install python-biopython`

If you want the documentation and unit tests,

`sudo apt-get install python-biopython-doc`

And if you want to use [BioSQL](BioSQL "wikilink"),

`sudo apt-get install python-biopython-sql`

However, this will probably not be the latest release (see [Ubuntu
listing here](http://packages.ubuntu.com/python-biopython), and [Debian
listing
here](http://packages.debian.org/search?searchon=sourcenames&keywords=biopython)).
If you want the latest version of Biopython, you will need to install it
from source. However, you should be able to automatically install the
build dependencies with the following command:

`sudo apt-get build-dep python-biopython`

Note: You may need to additionally install the NumPy package by hand, as
a very out of date repository may still expect Biopython to use Numeric
instead.

### Fedora

Biopython is an official Fedora package (since Fedora 5). The package is
named
[python-biopython](https://admin.fedoraproject.org/community/?package=python-biopython#package_maintenance),
and can be installed using yum as root:

`yum install python-biopython`

or via one of the GUI package management systems such as pirut and
PackageKit (available in F-9 and later).

### Gentoo Linux

Gentoo's portage tree contains an ebuild (sci-biology/biopython) which
builds from source. To install it, open a terminal as root and run:

`emerge -va biopython `

[Here](http://www.gentoo-portage.com/sci-biology/biopython) is a link to
Biopython at [Gentoo-Portage](http://www.gentoo-portage.com) which shows
the latest versions in Gentoo's Portage tree.

Ports
=====

### FreeBSD

The most easy way of installing Biopython in
[FreeBSD](http://www.freebsd.org/) is through the [Ports
Collection](http://www.freebsd.org/ports/). If you're new to this
procedure please take a look at [this
document](http://www.freebsd.org/doc/en_US.ISO8859-1/books/handbook/ports-using.html).
Supposing that you're familiar with this method and that you have an
up-to-date ports tree, all you need to do is to execute the following
commands as root:

<bash> cd /usr/ports/biology/py-biopython make install clean </bash>

Due to the great architecture of the ports system, this simple commands
will automatically fetch and install Biopython (as well as its necessary
dependencies).

Old Releases
============

Recent releases of Biopython require NumPy (and not Numeric):

-   [biopython-1.56.tar.gz](http://biopython.org/DIST/biopython-1.56.tar.gz)
    6,778 Kb -- Source Tarball (26 November 2010)
-   [biopython-1.56.zip](http://biopython.org/DIST/biopython-1.56.zip)
    7,347 Kb -- Source Zip File
-   [biopython-1.56.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.56.win32-py2.4.exe)
    1,429 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.56.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.56.win32-py2.5.exe)
    1,429 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.56.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.56.win32-py2.6.exe)
    1,451 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.56.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.56.win32-py2.7.exe)
    1,578 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5
-   [biopython-1.55.tar.gz](http://biopython.org/DIST/biopython-1.55.tar.gz)
    6,493 Kb -- Source Tarball (31 August 2010)
-   [biopython-1.55.zip](http://biopython.org/DIST/biopython-1.55.zip)
    7,058 Kb -- Source Zip File
-   [biopython-1.55.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.55.win32-py2.4.exe)
    1,448 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.55.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.55.win32-py2.5.exe)
    1,449 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.55.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.55.win32-py2.6.exe)
    1,471 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.55.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.55.win32-py2.7.exe)
    1,598 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5

<!-- -->

-   [biopython-1.55b.tar.gz](http://biopython.org/DIST/biopython-1.55b.tar.gz)
    6,428 Kb -- Source Tarball (August 18, 2010)
-   [biopython-1.55b.zip](http://biopython.org/DIST/biopython-1.55b.zip)
    6,996 Kb -- Source Zip File
-   [biopython-1.55b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.55b.win32-py2.4.exe)
    1,451 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.55b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.55b.win32-py2.5.exe)
    1,451 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.55b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.55b.win32-py2.6.exe)
    1,474 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.54.tar.gz](http://biopython.org/DIST/biopython-1.54.tar.gz)
    6,295 Kb -- Source Tarball (May 20, 2010)
-   [biopython-1.54.zip](http://biopython.org/DIST/biopython-1.54.zip)
    6,859 Kb -- Source Zip File
-   [biopython-1.54.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.54.win32-py2.4.exe)
    1,434 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.54.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.54.win32-py2.5.exe)
    1,434 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.54.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.54.win32-py2.6.exe)
    1,457 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.54b.tar.gz](http://biopython.org/DIST/biopython-1.54b.tar.gz)
    6,554 Kb -- Source Tarball (April 2, 2010)
-   [biopython-1.54b.zip](http://biopython.org/DIST/biopython-1.54b.zip)
    7,118 Kb -- Source Zip File
-   [biopython-1.54b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.54b.win32-py2.4.exe)
    1,426 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.54b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.54b.win32-py2.5.exe)
    1,427 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.54b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.54b.win32-py2.6.exe)
    1,456 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.53.tar.gz](http://biopython.org/DIST/biopython-1.53.tar.gz)
    4,185 Kb -- Source Tarball (December 15, 2009)
-   [biopython-1.53.zip](http://biopython.org/DIST/biopython-1.53.zip)
    4,652 Kb -- Source Zip File
-   [biopython-1.53.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.53.win32-py2.4.exe)
    1,129 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.53.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.53.win32-py2.5.exe)
    1,130 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.53.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.53.win32-py2.6.exe)
    1,155 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.52.tar.gz](http://biopython.org/DIST/biopython-1.52.tar.gz)
    5,486 Kb -- Source Tarball (September 22, 2009)
-   [biopython-1.52.zip](http://biopython.org/DIST/biopython-1.52.zip)
    5,930 Kb -- Source Zip File
-   [biopython-1.52.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.52.win32-py2.4.exe)
    1,107 Kb -- 32 bit Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.52.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.52.win32-py2.5.exe)
    1,108 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.52.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.52.win32-py2.6.exe)
    1,147 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.51.tar.gz](http://biopython.org/DIST/biopython-1.51.tar.gz)
    5,428 Kb -- Source Tarball (August 17, 2009)
-   [biopython-1.51.zip](http://biopython.org/DIST/biopython-1.51.zip)
    5,922 Kb -- Source Zip File
-   [biopython-1.51.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.51.win32-py2.4.exe)
    1,166 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.51.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.51.win32-py2.5.exe)
    1,167 Kb -- Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.51.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.51.win32-py2.6.exe)
    1,206 Kb -- Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.51b.tar.gz](http://biopython.org/DIST/biopython-1.51b.tar.gz)
    5,172 Kb -- Source Tarball (June 23, 2009)
-   [biopython-1.51b.zip](http://biopython.org/DIST/biopython-1.51b.zip)
    5,605 Kb -- Source Zip File
-   [biopython-1.51b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.51b.win32-py2.4.exe)
    1,161 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.51b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.51b.win32-py2.5.exe)
    1,161 Kb -- Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.51b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.51b.win32-py2.6.exe)
    1,199 Kb -- Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.50.tar.gz](http://biopython.org/DIST/biopython-1.50.tar.gz)
    4,550 Kb -- Source Tarball (April 20, 2009)
-   [biopython-1.50.zip](http://biopython.org/DIST/biopython-1.50.zip)
    4,988 Kb -- Source Zip File
-   [biopython-1.50.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.50.win32-py2.3.exe)
    1,228 Kb -- Windows Installer for Python 2.3 and NumPy 1.1
-   [biopython-1.50.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.50.win32-py2.4.exe)
    1,232 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.50.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.50.win32-py2.5.exe)
    1,232 Kb -- Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.50.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.50.win32-py2.6.exe)
    1,270 Kb -- Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.50b.tar.gz](http://biopython.org/DIST/biopython-1.50b.tar.gz)
    4,788 Kb (April 3, 2009)
-   [biopython-1.50b.zip](http://biopython.org/DIST/biopython-1.50b.zip)
    5,250 Kb
-   [biopython-1.50b.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.50b.win32-py2.3.exe)
    1,226 Kb -- Windows Installer for Python 2.3 and NumPy 1.1
-   [biopython-1.50b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.50b.win32-py2.4.exe)
    1,230 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.50b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.50b.win32-py2.5.exe)
    1,230 Kb -- Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.50b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.50b.win32-py2.6.exe)
    1,268 Kb -- Windows Installer for Python 2.6 and NumPy 1.3

<!-- -->

-   [biopython-1.49.tar.gz](http://biopython.org/DIST/biopython-1.49.tar.gz)
    4,052 Kb (November 21, 2008)
-   [biopython-1.49.zip](http://biopython.org/DIST/biopython-1.49.zip)
    4,498 Kb
-   [biopython-1.49.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.49.win32-py2.3.exe)
    1,111 Kb -- Windows Installer for Python 2.3 and NumPy 1.1
-   [biopython-1.49.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.49.win32-py2.4.exe)
    1,115 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.49.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.49.win32-py2.5.exe)
    1,115 Kb -- Windows Installer for Python 2.5 and NumPy 1.1

<!-- -->

-   [biopython-1.49b.tar.gz](http://biopython.org/DIST/biopython-1.49b.tar.gz)
    4,331 Kb (November 7, 2008)
-   [biopython-1.49b.zip](http://biopython.org/DIST/biopython-1.49b.zip)
    4,780 Kb
-   [biopython-1.49b.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.49b.win32-py2.3.exe)
    1,109 Kb -- Windows Installer for Python 2.3 and NumPy 1.1
-   [biopython-1.49b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.49b.win32-py2.4.exe)
    1,113 Kb -- Windows Installer for Python 2.4 and NumPy 1.1
-   [biopython-1.49b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.49b.win32-py2.5.exe)
    1,114 Kb -- Windows Installer for Python 2.5 and NumPy 1.1

Please note that Biopython 1.48 and older require the Numeric library,
not its replacement NumPy. Windows installers for Python 2.4 and older
are available from the [Numerical
Python](http://numpy.scipy.org/#older_array) website. A Windows
installer for Numeric 24.2 for Python 2.5 is available here:
[Numeric-24.2.win32-py2.5.exe](http://biopython.org/DIST/Numeric-24.2.win32-py2.5.exe)
446 Kb

Please note that Biopython 1.48 and older used [mxTextTools
2.0](http://www.egenix.com/www2002/python/eGenix-mx-Extensions-v2.x.html/)
in some of the parsers. There were a few niggles with mxTextTools 3.0,
so ideally install the older mxTextTools 2.0.

-   [biopython-1.48.tar.gz](http://biopython.org/DIST/biopython-1.48.tar.gz)
    4,051 Kb (September 8, 2008)
-   [biopython-1.48.zip](http://biopython.org/DIST/biopython-1.48.zip)
    4,542 Kb
-   [biopython-1.48.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.48.win32-py2.3.exe)
    1,226 Kb
-   [biopython-1.48.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.48.win32-py2.4.exe)
    1,254 Kb
-   [biopython-1.48.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.48.win32-py2.5.exe)
    1,254 Kb
-   [biopython-1.47.tar.gz](http://biopython.org/DIST/biopython-1.47.tar.gz)
    4,018 Kb (July 5, 2008)
-   [biopython-1.47.zip](http://biopython.org/DIST/biopython-1.47.zip)
    4,528 Kb
-   [biopython-1.47.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.47.win32-py2.3.exe)
    1,207 Kb
-   [biopython-1.47.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.47.win32-py2.4.exe)
    1,236 Kb
-   [biopython-1.47.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.47.win32-py2.5.exe)
    1,236 Kb
-   [biopython-1.46.tar.gz](http://biopython.org/DIST/biopython-1.46.tar.gz)
    3,926 Kb (June 29, 2008)
-   [biopython-1.46.zip](http://biopython.org/DIST/biopython-1.46.zip)
    4,426 Kb
-   [biopython-1.46.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.46.win32-py2.3.exe)
    1,206 Kb
-   [biopython-1.46.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.46.win32-py2.4.exe)
    1,235 Kb
-   [biopython-1.46.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.46.win32-py2.5.exe)
    1,235 Kb
-   [biopython-1.45.tar.gz](http://biopython.org/DIST/biopython-1.45.tar.gz)
    3,886 Kb (March 22, 2008)
-   [biopython-1.45.zip](http://biopython.org/DIST/biopython-1.45.zip)
    4,395 Kb
-   [biopython-1.45.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.45.win32-py2.3.exe)
    1,113 Kb
-   [biopython-1.45.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.45.win32-py2.4.exe)
    1,141 Kb
-   [biopython-1.45.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.45.win32-py2.5.exe)
    1,142 Kb
-   [biopython-1.44.tar.gz](http://biopython.org/DIST/biopython-1.44.tar.gz)
    3,750 Kb (October 28, 2007)
-   [biopython-1.44.zip](http://biopython.org/DIST/biopython-1.44.zip)
    4,243 Kb
-   [biopython-1.44.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.3.exe)
    1,091 Kb
-   [biopython-1.44.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.4.exe)
    1,116 Kb
-   [biopython-1.44.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.5.exe)
    1,116 Kb
-   [biopython-1.43.tar.gz](http://biopython.org/DIST/biopython-1.43.tar.gz)
    3,778 Kb (March 17, 2007)
-   [biopython-1.43.zip](http://biopython.org/DIST/biopython-1.43.zip)
    4,271 Kb
-   [biopython-1.43.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.3.exe)
    1,104 Kb
-   [biopython-1.43.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.4.exe)
    1,108 Kb
-   [biopython-1.43.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.5.exe)
    1,109 Kb
-   [biopython-1.42.tar.gz](http://biopython.org/DIST/biopython-1.42.tar.gz)
    3,841 Kb (July 16, 2006)
-   [biopython-1.42.zip](http://biopython.org/DIST/biopython-1.42.zip)
    4,399 Kb
-   [biopython-1.42.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.3.exe)
    1,070 Kb
-   [biopython-1.42.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.4.exe)
    1,074 Kb
-   [biopython-1.42.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.5.exe)
    1,075 Kb
-   [biopython-1.41.tar.gz](http://biopython.org/DIST/biopython-1.41.tar.gz)
    3,719 Kb (October 28, 2005)
-   [biopython-1.41.zip](http://biopython.org/DIST/biopython-1.41.zip)
    4,241 Kb
-   [biopython-1.41.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.3.exe)
    1,038 Kb
-   [biopython-1.41.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.4.exe)
    1,042 Kb
-   [biopython-1.40b.tar.gz](http://biopython.org/DIST/biopython-1.40b.tar.gz)
    3,437 Kb (February 18, 2005)
-   [biopython-1.40b.zip](http://biopython.org/DIST/biopython-1.40b.zip)
    3,267 Kb
-   [biopython-1.40b.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.3.exe)
    1,019 Kb
-   [biopython-1.40b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.4.exe)
    1,023 Kb
-   [biopython-1.30.tar.gz](http://biopython.org/DIST/biopython-1.30.tar.gz)
    3,186 Kb (May 14, 2004)
-   [biopython-1.24.tar.gz](http://biopython.org/DIST/biopython-1.24.tar.gz)
    3,081 Kb (February 16, 2004)
-   [biopython-1.24.zip](http://biopython.org/DIST/biopython-1.24.zip)
    3,623 Kb
-   [biopython-1.24.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.2.exe)
    892 Kb
-   [biopython-1.24.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.3.exe)
    894 Kb
-   [biopython-1.23.tar.gz](http://biopython.org/DIST/biopython-1.23.tar.gz)
    2,241 Kb (October 18, 2003)
-   [biopython-1.23.zip](http://biopython.org/DIST/biopython-1.23.zip)
    2,719 Kb
-   [biopython-1.23.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.2.exe)
    833 Kb
-   [biopython-1.23.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.3.exe)
    842 Kb
-   [biopython-1.22.tar.gz](http://biopython.org/DIST/biopython-1.22.tar.gz)
    2,214 Kb (October 9, 2003)
-   [biopython-1.22.zip](http://biopython.org/DIST/biopython-1.22.zip)
    2,691 Kb
-   [biopython-1.21.tar.gz](http://biopython.org/DIST/biopython-1.21.tar.gz)
    2,214 Kb
-   [biopython-1.21.zip](http://biopython.org/DIST/biopython-1.21.zip)
    2,897 Kb
-   [biopython-1.21.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.21.win32-py2.2.exe)
    770 Kb
-   [biopython-1.21.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.21.win32-py2.3.exe)
    832 Kb
-   [biopython-1.20.tar.gz](http://biopython.org/DIST/biopython-1.20.tar.gz)
    2,101 Kb (July 28, 2003)
-   [biopython-1.20.zip](http://biopython.org/DIST/biopython-1.20.zip)
    2,602 Kb
-   [biopython-1.10.tar.gz](http://biopython.org/DIST/biopython-1.10.tar.gz)
    1,811 Kb (December 17, 2002)
-   [biopython-1.10.zip](http://biopython.org/DIST/biopython-1.10.zip)
    2,300 Kb
-   [biopython-1.10.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.10.win32-py2.2.exe)
    1,199 Kb
-   [biopython-1.00a4.tar.gz](http://biopython.org/DIST/biopython-1.00a4.tar.gz)
    1,739Kb (December 18, 2001)
-   [biopython-1.00a4.zip](http://biopython.org/DIST/biopython-1.00a4.zip)
    2,121Kb
-   [biopython-1.00a4.win32-py2.0.exe](http://biopython.org/DIST/biopython-1.00a4.win32-py2.0.exe)
    835Kb
-   [biopython-1.00a4.win32-py2.1.exe](http://biopython.org/DIST/biopython-1.00a4.win32-py2.1.exe)
    837Kb
-   [biopython-1.00a4.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.00a4.win32-py2.2.exe)
    838Kb
-   [MacBiopython-1.00a4.sit](http://biopython.org/DIST/MacBiopython-1.00a4.sit)
    2.2Mb
-   [biopython-1.00a3.tar.gz](http://biopython.org/DIST/biopython-1.00a3.tar.gz)
    1,816Kb (September 3, 2001)
-   [biopython-1.00a3.zip](http://biopython.org/DIST/biopython-1.00a3.zip)
    2,165Kb
-   [biopython-1.00a3.win32-py2.0.exe](http://biopython.org/DIST/biopython-1.00a3.win32-py2.0.exe)
    583Kb
-   [biopython-1.00a3.win32-py2.1.exe](http://biopython.org/DIST/biopython-1.00a3.win32-py2.1.exe)
    585Kb
-   [Macbiopython-1.00a3.sit.bin](http://biopython.org/DIST/Macbiopython-1.00a3.sit.bin)
    1926Kb

