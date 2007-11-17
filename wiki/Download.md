---
title: Download
permalink: wiki/Download
layout: wiki
---

Current Release -- 1.44 -- 28 October 2007
==========================================

See also [What's new](http://biopython.open-bio.org/SRC/biopython/NEWS).

### Installation Instructions

-   [HTML Installation
    Instructions](http://biopython.org/DIST/docs/install/Installation.html)
-   [PDF Installation
    Instructions](http://biopython.org/DIST/docs/install/Installation.pdf)

### Files

-   Source Tarball --
    [biopython-1.44.tar.gz](http://biopython.org/DIST/biopython-1.44.tar.gz)
    3,750 Kb
-   Source Zip File --
    [biopython-1.44.zip](http://biopython.org/DIST/biopython-1.44.zip)
    4,243 Kb
-   Windows Installer for Python 2.3 --
    [biopython-1.44.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.3.exe)
    1,091 Kb
-   Windows Installer for Python 2.4 --
    [biopython-1.44.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.4.exe)
    1,116 Kb
-   Windows Installer for Python 2.5 --
    [biopython-1.44.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.44.win32-py2.5.exe)
    1,116 Kb

### Required Software

-   [Python 2.3 or above](http://www.python.org)
-   [A C compiler (if compiling
    from source)](http://www.python.org/doc/current/inst/inst.html) You
    need a C compiler supported by distutils, gcc will work fine on
    UNIX-like platforms. This is not needed on Windows if using the
    install programs provided above.
-   [mxTextTools
    2.0](http://www.egenix.com/www2002/python/eGenix-mx-Extensions-v2.x.html/)
    There are a few niggles with mxTextTools 3.0, so ideally install the
    older mxTextTools 2.0.
-   [Numerical Python](http://numpy.scipy.org/#older_array) Note that
    BioPython has not (yet) switched to the 'new' numpy library. You
    need the 'old' Numeric library, version 24.2 is recommended. Windows
    installers for Python 2.4 and older are available from the
    [Numerical Python](http://numpy.scipy.org/#older_array) website. A
    Windows installer for Numeric 24.2 for Python 2.5 is available here:
    [Numeric-24.2.win32-py2.5.exe](http://biopython.org/DIST/Numeric-24.2.win32-py2.5.exe)
    446 Kb

For Ubuntu Linux, install the following packages to cover these
requirements:

-   python-egenix-mxtexttools
-   python-numeric
-   python-dev
-   build-essential

### Optional Software

-   [ReportLab](http://www.reportlab.org/downloads.html) -- used for pdf
    graphics code
-   [MySQLdb](http://sourceforge.net/projects/mysql-python) -- used for
    [BioSQL](BioSQL "wikilink")
-   [flex: The Fast Lexical Analyzer](http://flex.sourceforge.net/) --
    for building Bio.PDB.mmCIF.MMCIFlex
-   [Wise2](http://www.ebi.ac.uk/Wise2/) -- for command line tool dnal
-   [NCBI Standalone
    BLAST](http://www.ncbi.nlm.nih.gov/blast/download.shtml) -- for
    running BLAST on your local machine
-   [Clustalw](ftp://ftp.ebi.ac.uk/pub/software/unix/clustalw/) --
    command line tool for building sequence alignments

Biopython-corba
===============

### Files

-   Source tarball --
    [biopython-corba-0.3.0.tar.gz](http://biopython.org/DIST/biopython-corba-0.3.0.tar.gz)
    344Kb
-   Source Zip File --
    [biopython-corba-0.3.0.zip](http://biopython.org/DIST/biopython-corba-0.3.0.zip)
    378Kb

### Required Software

-   [omniORBpy](http://omniorb.sourceforge.net)
-   [Fnorb](http://fnorb.sourceforge.net)
-   [orbit-python](http://sourceforge.net/projects/orbit-python)

Packages
========

For those of you using Linux, the easiest way to install Biopython is
through your distribution's package management system. However, unless
you are running a recent release of your Linux Distribution, you may
find that the Biopython packages available to be a little out of date.
You might want to see if there is a backport available, otherwise you
will have to install Biopython from source.

### Ubuntu

You should be able to install Biopython and its dependencies using the
Synaptic GUI tool (on the main menu under System / Administration /
Synaptic Package Manager), or at the command line using:

`sudo apt-get install python-biopython`

However, this will probably not be the latest release.

Biopython 1.41 packages for Ubuntu (Dapper):

-   [python-biopython](http://packages.ubuntulinux.org/dapper/source/python-biopython)

Ubuntu Edgy doesn't seem to have working AMD64 packages. Fiesty
biopython 1.42 AMD64 packages backported for Edgy here:

-   [python-biopython-doc\_1.42-2\_all.deb](http://students.ee.sun.ac.za/~nmarais/files/python-biopython-doc_1.42-2_all.deb)
-   [python-biopython-martel\_1.42-2\_all.deb](http://students.ee.sun.ac.za/~nmarais/files/python-biopython-martel_1.42-2_all.deb)
-   [python-biopython-sql\_1.42-2\_all.deb](http://students.ee.sun.ac.za/~nmarais/files/python-biopython-sql_1.42-2_all.deb)
-   [python-biopython\_1.42-2\_amd64.deb](http://students.ee.sun.ac.za/~nmarais/files/python-biopython_1.42-2_amd64.deb)

### Fedora

[Fedora Core
5](http://download.fedora.redhat.com/pub/fedora/linux/extras/5/i386/repoview/python-biopython.html),
[Fedora Core
6](http://download.fedora.redhat.com/pub/fedora/linux/extras/6/i386/repoview/python-biopython.html)
and
[rawhide](http://download.fedora.redhat.com/pub/fedora/linux/extras/development/i386/repoview/python-biopython.html)
(soon-to-be Fedora 7) have official packages for biopython 1.43 called
python-biopython, it can be installed using yum:

`yum install python-biopython`

Ports
=====

### FreeBSD

The most easy way of installing BioPython in
[FreeBSD](http://www.freebsd.org/) is through the [Ports
Collection](http://www.freebsd.org/ports/). If you're new to this
procedure please take a look at [this
document](http://www.freebsd.org/doc/en_US.ISO8859-1/books/handbook/ports-using.html).
Supposing that you're familiar with this method and that you have an
up-to-date ports tree, all you need to do is to execute the following
commands as root:

<bash>

1.  cd /usr/ports/biology/py-biopython
2.  make install clean

</bash>

Due to the great architecture of the ports system, this simple commands
will automatically fetch and install BioPython (as well as its necessary
dependencies).

Old Releases
============

-   [biopython-1.43.tar.gz](http://biopython.org/DIST/biopython-1.43.tar.gz)
    3,778 Kb
-   [biopython-1.43.zip](http://biopython.org/DIST/biopython-1.43.zip)
    4,271 Kb
-   [biopython-1.43.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.3.exe)
    1,104 Kb
-   [biopython-1.43.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.4.exe)
    1,108 Kb
-   [biopython-1.43.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.43.win32-py2.5.exe)
    1,109 Kb
-   [biopython-1.42.tar.gz](http://biopython.org/DIST/biopython-1.42.tar.gz)
    3,841 Kb
-   [biopython-1.42.zip](http://biopython.org/DIST/biopython-1.42.zip)
    4,399 Kb
-   [biopython-1.42.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.3.exe)
    1,070 Kb
-   [biopython-1.42.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.4.exe)
    1,074 Kb
-   [biopython-1.42.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.42.win32-py2.5.exe)
    1,075 Kb
-   [biopython-1.41.tar.gz](http://biopython.org/DIST/biopython-1.41.tar.gz)
    3,719 Kb
-   [biopython-1.41.zip](http://biopython.org/DIST/biopython-1.41.zip)
    4,241 Kb
-   [biopython-1.41.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.3.exe)
    1,038 Kb
-   [biopython-1.41.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.4.exe)
    1,042 Kb
-   [biopython-1.40b.tar.gz](http://biopython.org/DIST/biopython-1.40b.tar.gz)
    3,437 Kb
-   [biopython-1.40b.zip](http://biopython.org/DIST/biopython-1.40b.zip)
    3,267 Kb
-   [biopython-1.40b.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.3.exe)
    1,019 Kb
-   [biopython-1.40b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.4.exe)
    1,023 Kb
-   [biopython-1.30.tar.gz](http://biopython.org/DIST/biopython-1.30.tar.gz)
    3,186 Kb
-   [biopython-1.24.tar.gz](http://biopython.org/DIST/biopython-1.24.tar.gz)
    3,081 Kb
-   [biopython-1.24.zip](http://biopython.org/DIST/biopython-1.24.zip)
    3,623 Kb
-   [biopython-1.24.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.2.exe)
    892 Kb
-   [biopython-1.24.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.3.exe)
    894 Kb
-   [biopython-1.23.tar.gz](http://biopython.org/DIST/biopython-1.23.tar.gz)
    2,241 Kb
-   [biopython-1.23.zip](http://biopython.org/DIST/biopython-1.23.zip)
    2,719 Kb
-   [biopython-1.23.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.2.exe)
    833 Kb
-   [biopython-1.23.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.3.exe)
    842 Kb
-   [biopython-1.22.tar.gz](http://biopython.org/DIST/biopython-1.22.tar.gz)
    2,214 Kb
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
    2,101 Kb
-   [biopython-1.20.zip](http://biopython.org/DIST/biopython-1.20.zip)
    2,602 Kb
-   [biopython-1.10.tar.gz](http://biopython.org/DIST/biopython-1.10.tar.gz)
    1,811 Kb
-   [biopython-1.10.zip](http://biopython.org/DIST/biopython-1.10.zip)
    2,300 Kb
-   [biopython-1.10.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.10.win32-py2.2.exe)
    1,199 Kb
-   [biopython-1.00a4.tar.gz](http://biopython.org/DIST/biopython-1.00a4.tar.gz)
    1,739Kb
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
    1,816Kb
-   [biopython-1.00a3.zip](http://biopython.org/DIST/biopython-1.00a3.zip)
    2,165Kb
-   [biopython-1.00a3.win32-py2.0.exe](http://biopython.org/DIST/biopython-1.00a3.win32-py2.0.exe)
    583Kb
-   [biopython-1.00a3.win32-py2.1.exe](http://biopython.org/DIST/biopython-1.00a3.win32-py2.1.exe)
    585Kb
-   [Macbiopython-1.00a3.sit.bin](http://biopython.org/DIST/Macbiopython-1.00a3.sit.bin)
    1926Kb

