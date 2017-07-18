---
title: Download
permalink: wiki/Download
layout: page
redirect_from:
 - /download/
---

Current Release - 1.70 - 10 July 2017
=====================================

See also [What's
new](https://github.com/biopython/biopython/blob/master/NEWS.rst).

### Files

#### Biopython 1.70

-   [biopython-1.70.tar.gz](http://biopython.org/DIST/biopython-1.70.tar.gz)
    15Mb -- Source Tarball
-   [biopython-1.70.zip](http://biopython.org/DIST/biopython-1.70.zip)
    16Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.70)

### Installation Instructions

Recent versions of Python (starting with Python 2.7.9 and Python 3.4) include
the Python package management tool ``pip``, which allows an easy installation
from the command line on all platforms. Try:

``` bash
pip install biopython
```

If pip is not already installed, try:

``` bash
python -m ensurepip
```

If you need to install under a specific version of Python, try something
like this:

``` bash
python2.7 -m pip install biopython
python3.6 -m pip install biopython
pypy -m pip install biopython
```

On **Windows**, by default ``python`` and ``pip`` are not on the ``PATH``.
You can re-install Python and tick this option, or give the full path instead.
Try something like this, depending on where your copy of Python is installed:

```
C:\Python27\Scripts\pip install biopython
```

### Installation from Source

Installation from source requires an appropriate C compiler, for example
GCC on **Linux**, and MSVC on **Windows**.
For **Mac OS X**, or as it is now branded, **macOS**, if you want to
compile Biopython from source you will need to have installed Apple's
command line tools, which can de done with the terminal command:

``` bash
xcode-select --install
```

This will offer to install Apple's XCode development suite - you can, but
it is not needed and takes a lot of disk space.

You can then download and unzip a Biopython source code release, or get
our code from GitHub. Then run:

``` bash
pip install .
```

If you are still stuck, sign up to the [Biopython mailing
list](Mailing_lists "wikilink") and ask for help there.

### Required Software

-   [Python 2.7, 3.4, 3.5, or 3.6](http://www.python.org)
-   [C compiler (if compiling from
    source)](https://docs.python.org/3/using/index.html) You
    need a C compiler supported by ``setuptools``, **gcc** will work fine on
    UNIX-like platforms. This is not needed on Windows if using the
    install programs provided above. On Mac OS, you should install
    Apple's the compiler tools as described above.
-   [NumPy (Numerical Python)](http://numpy.scipy.org/).

### Optional Software

Some parts of Biopython use the following additional python libraries:

-   [ReportLab](http://www.reportlab.com/software/downloads/) -- used for pdf
    graphics code
-   [psycopg](http://initd.org/psycopg/) -- used for
    [BioSQL](BioSQL "wikilink") with a PostgreSQL database
-   [mysql-connector](http://dev.mysql.com/downloads/connector/python/)
    -- used for [BioSQL](BioSQL "wikilink") with a MySQL database
-   [MySQLdb](http://sourceforge.net/projects/mysql-python) -- An
    alternative MySQL library used by [BioSQL](BioSQL "wikilink")

In addition Biopython includes wrapper code for calling a number of
third party command line tools including:

-   [Wise2](http://www.ebi.ac.uk/~birney/wise2/) -- for command line tool dnal
-   [NCBI Standalone
    BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) -- command
    line tool for running BLAST on your local machine
-   [Clustalw](ftp://ftp.ebi.ac.uk/pub/software/unix/clustalw/) --
    command line tool for building sequence alignments
-   [SIMCOAL2](http://cmpg.unibe.ch/software/simcoal2/) and
    [FDist](http://www.maths.bris.ac.uk/~mamab/software/) -- command
    line tools for population genetics
-   [EMBOSS](http://emboss.sourceforge.net/) -- lots of useful command
    line tools.

Packages
========

We would now recommend ``pip``, however for those of you using Linux,
you can alternatively install Biopython through your Linux
distribution's package management system. However, unless
you are running a recent release of your Linux Distribution, you may
find that the Biopython packages available to be a little out of date.
You might want to see if there is a backport available, otherwise you
will have to install Biopython using ``pip`` or from compiled from source.

### Ubuntu or Debian

You should be able to install Biopython and its dependencies using the
Synaptic GUI tool (on the main menu under System / Administration /
Synaptic Package Manager), or at the command line using:

``` bash
sudo apt-get install python-biopython
```

If you want the documentation and unit tests,

``` bash
sudo apt-get install python-biopython-doc
```

And if you want to use [BioSQL](BioSQL "wikilink"),

``` bash
sudo apt-get install python-biopython-sql
```

However, this will probably not be the latest release (see [Ubuntu
listing here](http://packages.ubuntu.com/python-biopython), and [Debian
listing
here](http://packages.debian.org/search?searchon=sourcenames&keywords=biopython)).
If you want the latest version of Biopython, you will need to install it
from source. However, you should be able to automatically install the
build dependencies with the following command:

``` bash
sudo apt-get build-dep python-biopython
```

### Archlinux

Biopython is in the [official Archlinux
repository](https://www.archlinux.org/packages/?q=biopython) as
python-biopython (for Python 3) or python2-biopython (for Python 2) and
can be installed using pacman:

``` bash
pacman -S python2-biopython
```

Or, for Python 3:

``` bash
pacman -S python-biopython
```

### Fedora

Biopython is an official Fedora package (since Fedora 5). The package is
named
[python-biopython](https://apps.fedoraproject.org/packages/python-biopython) for Python 2, or
[python3-biopython](https://apps.fedoraproject.org/packages/python3-biopython) for Python 3,
and can be installed using yum as root:

``` bash
yum install python-biopython
```

or

``` bash
yum install python3-biopython
```

or via one of the GUI package management systems such as pirut and
PackageKit (available in F-9 and later).

### Gentoo Linux

Gentoo's portage tree contains an ebuild (sci-biology/biopython) which
builds from source. To install it, open a terminal as root and run:

``` bash
emerge -va biopython
```

[Here](https://packages.gentoo.org/packages/sci-biology/biopython) is a link to
Biopython at [Gentoo](https://packages.gentoo.org/) which shows
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

``` bash
cd /usr/ports/biology/py-biopython make install clean
```

Due to the great architecture of the ports system, this simple commands
will automatically fetch and install Biopython (as well as its necessary
dependencies).

Old Releases
============

Recent releases of Biopython require NumPy (and not Numeric):

-   [biopython-1.70.tar.gz](http://biopython.org/DIST/biopython-1.70.tar.gz)
    15Mb -- Source Tarball
-   [biopython-1.70.zip](http://biopython.org/DIST/biopython-1.70.zip)
    16Mb -- Source Zip File

<!-- -->

-   [biopython-1.69.tar.gz](http://biopython.org/DIST/biopython-1.69.tar.gz)
    15Mb -- Source Tarball
-   [biopython-1.69.zip](http://biopython.org/DIST/biopython-1.69.zip)
    16Mb -- Source Zip File
-   [biopython-1.69.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.69.win32-py2.7.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.7 and NumPy 1.11.0
-   [biopython-1.69.win32-py2.7.msi](http://biopython.org/DIST/biopython-1.69.win32-py2.7.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 2.7 and NumPy 1.11.0
-   [biopython-1.69.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.69.win32-py3.3.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.3 and NumPy 1.10.2
-   [biopython-1.69.win32-py3.3.msi](http://biopython.org/DIST/biopython-1.69.win32-py3.3.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.3 and NumPy 1.10.2
-   [biopython-1.69.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.69.win32-py3.4.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.4 and NumPy 1.11.0
-   [biopython-1.69.win32-py3.4.msi](http://biopython.org/DIST/biopython-1.69.win32-py3.4.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.4 and NumPy 1.11.0
-   [biopython-1.69.win32-py3.5.exe](http://biopython.org/DIST/biopython-1.69.win32-py3.5.exe)
    3Mb -- 32 bit Windows .exe Installer for Python 3.5 and NumPy 1.11.1
-   [biopython-1.69.win32-py3.5.msi](http://biopython.org/DIST/biopython-1.69.win32-py3.5.msi)
    3Mb -- 32 bit Windows .msi Installer for Python 3.5 and NumPy 1.11.1
-   [biopython-1.69.win32-py3.6.exe](http://biopython.org/DIST/biopython-1.69.win32-py3.6.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.6 and NumPy 1.11.3
-   [biopython-1.69.win32-py3.6.msi](http://biopython.org/DIST/biopython-1.69.win32-py3.6.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.6 and NumPy 1.11.3

<!-- -->

-   [biopython-1.68.tar.gz](http://biopython.org/DIST/biopython-1.68.tar.gz)
    14Mb -- Source Tarball
-   [biopython-1.68.zip](http://biopython.org/DIST/biopython-1.68.zip)
    15Mb -- Source Zip File
-   [biopython-1.68.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.68.win32-py2.6.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.6 and NumPy 1.8.2
-   [biopython-1.68.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.68.win32-py2.7.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.7 and NumPy 1.11.0
-   [biopython-1.68.win32-py2.7.msi](http://biopython.org/DIST/biopython-1.68.win32-py2.7.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 2.7 and NumPy 1.11.0
-   [biopython-1.68.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.68.win32-py3.3.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.3 and NumPy 1.10.2
-   [biopython-1.68.win32-py3.3.msi](http://biopython.org/DIST/biopython-1.68.win32-py3.3.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.3 and NumPy 1.10.2
-   [biopython-1.68.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.68.win32-py3.4.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.4 and NumPy 1.11.0
-   [biopython-1.68.win32-py3.4.msi](http://biopython.org/DIST/biopython-1.68.win32-py3.4.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.4 and NumPy 1.11.0
-   [biopython-1.68.win32-py3.5.exe](http://biopython.org/DIST/biopython-1.68.win32-py3.5.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.5 and NumPy 1.11.1
-   [biopython-1.68.win32-py3.5.msi](http://biopython.org/DIST/biopython-1.68.win32-py3.5.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.5 and NumPy 1.11.1

<!-- -->

-   [biopython-1.67.tar.gz](http://biopython.org/DIST/biopython-1.67.tar.gz)
    14Mb -- Source Tarball
-   [biopython-1.67.zip](http://biopython.org/DIST/biopython-1.67.zip)
    15Mb -- Source Zip File
-   [biopython-1.67.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.67.win32-py2.6.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.6 and NumPy 1.8.2
-   [biopython-1.67.win32-py2.6.msi](http://biopython.org/DIST/biopython-1.67.win32-py2.6.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 2.6 and NumPy 1.8.2
-   [biopython-1.67.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.67.win32-py2.7.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.7 and NumPy 1.9.1
-   [biopython-1.67.win32-py2.7.msi](http://biopython.org/DIST/biopython-1.67.win32-py2.7.msi)
    2Mb -- 32 bit Windows .exe Installer for Python 2.7 and NumPy 1.9.1
-   [biopython-1.67.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.67.win32-py3.3.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.3 and NumPy 1.9.1
-   [biopython-1.67.win32-py3.3.msi](http://biopython.org/DIST/biopython-1.67.win32-py3.3.msi)
    2Mb -- 32 bit Windows .exe Installer for Python 3.3 and NumPy 1.9.1
-   [biopython-1.67.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.67.win32-py3.4.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.4 and NumPy 1.9.1
-   [biopython-1.67.win32-py3.4.msi](http://biopython.org/DIST/biopython-1.67.win32-py3.4.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.4 and NumPy 1.9.1
-   [biopython-1.67.win32-py3.5.exe](http://biopython.org/DIST/biopython-1.67.win32-py3.5.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.4 and NumPy 1.9.3
-   [biopython-1.67.win32-py3.5.msi](http://biopython.org/DIST/biopython-1.67.win32-py3.5.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.4 and NumPy 1.9.3

<!-- -->

-   [biopython-1.66.tar.gz](http://biopython.org/DIST/biopython-1.66.tar.gz)
    14Mb -- Source Tarball
-   [biopython-1.66.zip](http://biopython.org/DIST/biopython-1.66.zip)
    15Mb -- Source Zip File
-   [biopython-1.66.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.66.win32-py2.6.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.6 and NumPy 1.8.2
-   [biopython-1.66.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.66.win32-py2.7.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 2.7 and NumPy 1.9.1
-   [biopython-1.66.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.66.win32-py3.3.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.3 and NumPy 1.9.1
-   [biopython-1.66.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.66.win32-py3.4.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.4 and NumPy 1.9.1
-   [biopython-1.66.win32-py3.5.exe](http://biopython.org/DIST/biopython-1.66.win32-py3.5.exe)
    2Mb -- 32 bit Windows .exe Installer for Python 3.5
-   [biopython-1.66.win32-py3.5.msi](http://biopython.org/DIST/biopython-1.66.win32-py3.5.msi)
    2Mb -- 32 bit Windows .msi Installer for Python 3.5

<!-- -->

-   [biopython-1.65.tar.gz](http://biopython.org/DIST/biopython-1.65.tar.gz)
    13Mb -- Source Tarball
-   [biopython-1.65.zip](http://biopython.org/DIST/biopython-1.65.zip)
    14Mb -- Source Zip File
-   [biopython-1.65.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.65.win32-py2.6.exe)
    2Mb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.8.2
-   [biopython-1.65.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.65.win32-py2.7.exe)
    2Mb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.9.1
-   [biopython-1.65.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.65.win32-py3.3.exe)
    2Mb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.9.1
-   [biopython-1.65.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.65.win32-py3.4.exe)
    2Mb -- 32 bit Windows Installer for Python 3.4 and NumPy 1.9.1

<!-- -->

-   [biopython-1.64.tar.gz](http://biopython.org/DIST/biopython-1.64.tar.gz)
    12Mb -- Source Tarball
-   [biopython-1.64.zip](http://biopython.org/DIST/biopython-1.64.zip)
    13Mb -- Source Zip File
-   [biopython-1.64.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.64.win32-py2.6.exe)
    2Mb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.8.1
-   [biopython-1.64.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.64.win32-py2.7.exe)
    2Mb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.8.1
-   [biopython-1.64.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.64.win32-py3.3.exe)
    2Mb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.8.1
-   [biopython-1.64.win32-py3.4.exe](http://biopython.org/DIST/biopython-1.64.win32-py3.4.exe)
    2Mb -- 32 bit Windows Installer for Python 3.4 and NumPy 1.8.1

<!-- -->

-   [biopython-1.63.tar.gz](http://biopython.org/DIST/biopython-1.63.tar.gz)
    11Mb -- Source Tarball
-   [biopython-1.63.zip](http://biopython.org/DIST/biopython-1.63.zip)
    12Mb -- Source Zip File
-   [biopython-1.63.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.63.win32-py2.6.exe)
    2Mb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.7
-   [biopython-1.63.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.63.win32-py2.7.exe)
    2Mb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.7
-   [biopython-1.63.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.63.win32-py3.3.exe)
    2Mb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7

<!-- -->

-   [biopython-1.63b.tar.gz](http://biopython.org/DIST/biopython-1.63b.tar.gz)
    11,123 Kb -- Source Tarball
-   [biopython-1.63b.zip](http://biopython.org/DIST/biopython-1.63b.zip)
    12,111 Kb -- Source Zip File
-   [biopython-1.63b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.63b.win32-py2.6.exe)
    1,877 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.7
-   [biopython-1.63b.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.63b.win32-py2.7.exe)
    2,003 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.7
-   [biopython-1.63b.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.63b.win32-py3.3.exe)
    2,005 Kb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7

<!-- -->

-   [biopython-1.62.tar.gz](http://biopython.org/DIST/biopython-1.62.tar.gz)
    11,123 Kb -- Source Tarball
-   [biopython-1.62.zip](http://biopython.org/DIST/biopython-1.62.zip)
    12,111 Kb -- Source Zip File
-   [biopython-1.62.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.62.win32-py2.5.exe)
    1,852 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.62.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.62.win32-py2.6.exe)
    1,877 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.62.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.62.win32-py2.7.exe)
    2,003 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5
-   [biopython-1.62.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.62.win32-py3.3.exe)
    2,005 Kb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7

<!-- -->

-   [biopython-1.62b.tar.gz](http://biopython.org/DIST/biopython-1.62b.tar.gz)
    10,658 Kb -- Source Tarball (*beta release*, 15 July 2013)
-   [biopython-1.62b.zip](http://biopython.org/DIST/biopython-1.62b.zip)
    11,607 Kb -- Source Zip File
-   [biopython-1.62b.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.62b.win32-py2.5.exe)
    1,661 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.62b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.62b.win32-py2.6.exe)
    1,686 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.62b.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.62b.win32-py2.7.exe)
    1,813 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5
-   [biopython-1.62b.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.62b.win32-py3.3.exe)
    1,814 Kb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7

<!-- -->

-   [biopython-1.61.tar.gz](http://biopython.org/DIST/biopython-1.61.tar.gz)
    10,311 Kb -- Source Tarball (5 February 2013)
-   [biopython-1.61.zip](http://biopython.org/DIST/biopython-1.61.zip)
    11,198 Kb -- Source Zip File
-   [biopython-1.61.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.61.win32-py2.5.exe)
    1,612 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.61.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.61.win32-py2.6.exe)
    1,637 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.61.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.61.win32-py2.7.exe)
    1,764 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5
-   [biopython-1.61.win32-py3.2-beta.exe](http://biopython.org/DIST/biopython-1.61.win32-py3.2-beta.exe)
    1,757 Kb -- 32 bit Windows Installer for Python 3.2 and NumPy 1.5
    (*beta* status for testing)
-   [biopython-1.61.win32-py3.3-beta.exe](http://biopython.org/DIST/biopython-1.61.win32-py3.3-beta.exe)
    1,750 Kb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7
    (*beta* status for testing)

<!-- -->

-   [biopython-1.60.tar.gz](http://biopython.org/DIST/biopython-1.60.tar.gz)
    9,280 Kb -- Source Tarball (25 June 2012)
-   [biopython-1.60.zip](http://biopython.org/DIST/biopython-1.60.zip)
    10,051 Kb -- Source Zip File
-   [biopython-1.60.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.60.win32-py2.5.exe)
    1,469 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.60.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.60.win32-py2.6.exe)
    1,492 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.60.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.60.win32-py2.7.exe)
    1,1618 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5
-   [biopython-1.60.win32-py3.2-beta.exe](http://biopython.org/DIST/biopython-1.60.win32-py3.2-beta.exe)
    1,611 Kb -- 32 bit Windows Installer for Python 3.2 and NumPy 1.5
    (*beta* status for testing)

<!-- -->

-   [biopython-1.59.tar.gz](http://biopython.org/DIST/biopython-1.59.tar.gz)
    8,377 Kb -- Source Tarball (24 February 2012)
-   [biopython-1.59.zip](http://biopython.org/DIST/biopython-1.59.zip)
    9,127 Kb -- Source Zip File
-   [biopython-1.59.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.59.win32-py2.5.exe)
    1,440 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.59.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.59.win32-py2.6.exe)
    1,463 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.59.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.59.win32-py2.7.exe)
    1,590 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5

<!-- -->

-   [biopython-1.58.tar.gz](http://biopython.org/DIST/biopython-1.58.tar.gz)
    7,847 Kb -- Source Tarball (18 August 2011)
-   [biopython-1.58.zip](http://biopython.org/DIST/biopython-1.58.zip)
    8,474 Kb -- Source Zip File
-   [biopython-1.58.win32-py2.4-unsupported.exe](http://biopython.org/DIST/biopython-1.58.win32-py2.4-unsupported.exe)
    1,427 Kb -- 32 bit Windows Installer for Python 2.4 (which we no
    longer officially support) and NumPy 1.1
-   [biopython-1.58.win32-py2.5.exe](http://biopython.org/DIST/biopython-1.58.win32-py2.5.exe)
    1,428 Kb -- 32 bit Windows Installer for Python 2.5 and NumPy 1.1
-   [biopython-1.58.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.58.win32-py2.6.exe)
    1,450 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.3
-   [biopython-1.58.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.58.win32-py2.7.exe)
    1,577 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.5

<!-- -->

-   [biopython-1.57.tar.gz](http://biopython.org/DIST/biopython-1.57.tar.gz)
    6,783 Kb -- Source Tarball (2 April 2011)
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

<!-- -->

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

<!-- -->

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

-   [Numeric-24.2.win32-py2.5.exe](https://sourceforge.net/projects/numpy/files/Old%20Numeric/24.2/)
    446 Kb - Windows Installer for Python 2.5

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

