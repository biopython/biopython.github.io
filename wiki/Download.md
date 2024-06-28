---
title: Download
permalink: wiki/Download
layout: page
redirect_from:
 - /download/
---

Current Release - 1.84 - 28 June 2024
========================================

See also [What's
new](https://github.com/biopython/biopython/blob/master/NEWS.rst).

### Files

#### Biopython 1.84

-   [biopython-1.84.tar.gz](http://biopython.org/DIST/biopython-1.84.tar.gz)
    25Mb -- Source Tarball
-   [biopython-1.84.zip](http://biopython.org/DIST/biopython-1.84.zip)
    27Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.84)
-   [Documentation](https://biopython.org/docs/1.84/)

### Installation Instructions

All supported versions of Python include the Python package management
tool ``pip``, which allows an easy installation from the command line on
all platforms. Try:

``` bash
pip install biopython
```

For updating an older version of Biopython try:

``` bash
pip install biopython --upgrade
```

This will remove older versions of Biopython and NumPy before it installs
the recent versions.

Should you wish to uninstall Biopython:

```bash
pip uninstall biopython
```

If pip is not already installed you may need to update your Python, but first try:

``` bash
python -m ensurepip
```

If you need to install under a specific version of Python, try something
like this:

``` bash
python3.9 -m pip install biopython
pypy -m pip install biopython
```

On **Windows**, by default ``python`` and ``pip`` are not on the ``PATH``.
You can re-install Python and tick this option, or give the full path instead.
Try something like this, depending on where your copy of Python is installed:

```
C:\Python39\Scripts\pip install biopython
```

### Other packages

While we generally recommend using ``pip`` to install Biopython using
the wheel packages we provide on PyPI (as above), there are also
[Biopython packages for Conda, Linux, etc](Packages "wikilink").


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
python setup.py build
python setup.py test
python setup.py install
```

or:

``` bash
pip install .
```

If you are still stuck, sign up to the [Biopython mailing
list](Mailing_lists "wikilink") and ask for help there.

### Required Software

-   [Python 3.6, 3.7 or 3.8](http://www.python.org) or PyPy,
    including the Python development header files like ``python.h``
-   [C compiler (if compiling from
    source)](https://docs.python.org/3/using/index.html) You
    need a C compiler supported by ``setuptools``, **gcc** will work fine on
    UNIX-like platforms. This is not needed on Windows if using the
    compiled packages provided. On Mac OS, you should install
    Apple's the compiler tools as described above.
-   [NumPy (Numerical Python)](https://numpy.org/).

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


Old Releases
============

Recent releases of Biopython require NumPy (and not Numeric).
Version 1.76 is the last release to support Python 2.7 and 3.5,
all later releases require Python 3:

-   [biopython-1.83.tar.gz](http://biopython.org/DIST/biopython-1.83.tar.gz)
    19Mb -- Source Tarball
-   [biopython-1.83.zip](http://biopython.org/DIST/biopython-1.83.zip)
    20Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.83)
-   [Tutorial-1.83.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.83.pdf) -- Documentation

<!-- -->

-   [biopython-1.82.tar.gz](http://biopython.org/DIST/biopython-1.82.tar.gz)
    19Mb -- Source Tarball (22 December 2023)
-   [biopython-1.82.zip](http://biopython.org/DIST/biopython-1.82.zip)
    20Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1\
.82)
-   [Tutorial-1.82.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.82.p\
df) -- Documentation

<!-- -->

-   [biopython-1.81.tar.gz](http://biopython.org/DIST/biopython-1.81.tar.gz)
    17Mb -- Source Tarball (12 February 2023)
-   [biopython-1.81.zip](http://biopython.org/DIST/biopython-1.81.zip)
    19Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.81)
-   [Tutorial-1.81.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.81.pdf) -- Documentation

<!-- -->

-   [biopython-1.80.tar.gz](http://biopython.org/DIST/biopython-1.80.tar.gz)
    17Mb -- Source Tarball (18 November 2022)
-   [biopython-1.80.zip](http://biopython.org/DIST/biopython-1.80.zip)
    19Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.80)
-   [Tutorial-1.80.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.80.pdf) -- Documentation

<!-- -->

-   [biopython-1.79.tar.gz](http://biopython.org/DIST/biopython-1.79.tar.gz)
    16Mb -- Source Tarball (1 June 2021)
-   [biopython-1.79.zip](http://biopython.org/DIST/biopython-1.79.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.79)
-   [Tutorial-1.79.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.79.pdf) -- Documentation

<!-- -->

-   [biopython-1.78.tar.gz](http://biopython.org/DIST/biopython-1.78.tar.gz)
    16Mb -- Source Tarball (4 September 2020)
-   [biopython-1.78.zip](http://biopython.org/DIST/biopython-1.78.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.78)
-   [Tutorial-1.78.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.78.pdf) -- Documentation

<!-- -->

-   [biopython-1.77.tar.gz](http://biopython.org/DIST/biopython-1.77.tar.gz)
    16Mb -- Source Tarball (25 May 2020)
-   [biopython-1.77.zip](http://biopython.org/DIST/biopython-1.77.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.77)
-   [Tutorial-1.77.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.77.pdf) -- Documentation

<!-- -->

-   [biopython-1.76.tar.gz](http://biopython.org/DIST/biopython-1.76.tar.gz)
    16Mb -- Source Tarball (20 December 2019)
-   [biopython-1.75.zip](http://biopython.org/DIST/biopython-1.76.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.76)
-   [Tutorial-1.76.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.76.pdf) -- Documentation

<!-- -->

-   [biopython-1.75.tar.gz](http://biopython.org/DIST/biopython-1.75.tar.gz)
    16Mb -- Source Tarball (6 November 2019)
-   [biopython-1.75.zip](http://biopython.org/DIST/biopython-1.75.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.75)
-   [Tutorial-1.75.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.75.pdf) -- Documentation

<!-- -->

-   [biopython-1.74.tar.gz](http://biopython.org/DIST/biopython-1.74.tar.gz)
    16Mb -- Source Tarball (16 July 2019)
-   [biopython-1.74.zip](http://biopython.org/DIST/biopython-1.74.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.74)
-   [Tutorial-1.74.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.74.pdf) -- Documentation

<!-- -->

-   [biopython-1.73.tar.gz](http://biopython.org/DIST/biopython-1.73.tar.gz)
    15Mb -- Source Tarball (18 December 2018)
-   [biopython-1.73.zip](http://biopython.org/DIST/biopython-1.73.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.73)
-   [Tutorial-1.73.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.73.pdf) -- Documentation

<!-- -->

-   [biopython-1.72.tar.gz](http://biopython.org/DIST/biopython-1.72.tar.gz)
    16Mb -- Source Tarball (27 June 2018)
-   [biopython-1.72.zip](http://biopython.org/DIST/biopython-1.72.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.72)
-   [Tutorial-1.72.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.72.pdf) -- Documentation

<!-- -->

-   [biopython-1.71.tar.gz](http://biopython.org/DIST/biopython-1.71.tar.gz)
    16Mb -- Source Tarball (4 April 2018)
-   [biopython-1.71.zip](http://biopython.org/DIST/biopython-1.71.zip)
    17Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.71)
-   [Tutorial-1.71.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.71.pdf) -- Documentation

<!-- -->

-   [biopython-1.70.tar.gz](http://biopython.org/DIST/biopython-1.70.tar.gz)
    15Mb -- Source Tarball (11 July 2017)
-   [biopython-1.70.zip](http://biopython.org/DIST/biopython-1.70.zip)
    16Mb -- Source Zip File
-   [Pre-compiled wheel files on PyPI](https://pypi.python.org/pypi/biopython/1.70)
-   [Tutorial-1.70.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.70.pdf) -- Documentation

<!-- -->

-   [biopython-1.69.tar.gz](http://biopython.org/DIST/biopython-1.69.tar.gz)
    15Mb -- Source Tarball (7 April 2017)
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
-   [Tutorial-1.69.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.69.pdf) -- Documentation

<!-- -->

-   [biopython-1.68.tar.gz](http://biopython.org/DIST/biopython-1.68.tar.gz)
    14Mb -- Source Tarball (26 August 2016)
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
-   [Tutorial-1.68.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.68.pdf) -- Documentation

<!-- -->

-   [biopython-1.67.tar.gz](http://biopython.org/DIST/biopython-1.67.tar.gz)
    14Mb -- Source Tarball (8 June 2016)
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
-   [Tutorial-1.67.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.67.pdf) -- Documentation

<!-- -->

-   [biopython-1.66.tar.gz](http://biopython.org/DIST/biopython-1.66.tar.gz)
    14Mb -- Source Tarball (21 October 2015)
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
-   [Tutorial-1.66.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.66.pdf) -- Documentation

<!-- -->

-   [biopython-1.65.tar.gz](http://biopython.org/DIST/biopython-1.65.tar.gz)
    13Mb -- Source Tarball (17 December 2014)
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
-   [Tutorial-1.65.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.65.pdf) -- Documentation

<!-- -->

-   [biopython-1.64.tar.gz](http://biopython.org/DIST/biopython-1.64.tar.gz)
    12Mb -- Source Tarball (29 May 2014)
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
-   [Tutorial-1.64.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.64.pdf) -- Documentation

<!-- -->

-   [biopython-1.63.tar.gz](http://biopython.org/DIST/biopython-1.63.tar.gz)
    11Mb -- Source Tarball (6 December 2013)
-   [biopython-1.63.zip](http://biopython.org/DIST/biopython-1.63.zip)
    12Mb -- Source Zip File
-   [biopython-1.63.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.63.win32-py2.6.exe)
    2Mb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.7
-   [biopython-1.63.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.63.win32-py2.7.exe)
    2Mb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.7
-   [biopython-1.63.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.63.win32-py3.3.exe)
    2Mb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7
-   [Tutorial-1.63.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.63.pdf) -- Documentation

<!-- -->

-   [biopython-1.63b.tar.gz](http://biopython.org/DIST/biopython-1.63b.tar.gz)
    11,123 Kb -- Source Tarball (*beta release*, 12 November 2013)
-   [biopython-1.63b.zip](http://biopython.org/DIST/biopython-1.63b.zip)
    12,111 Kb -- Source Zip File
-   [biopython-1.63b.win32-py2.6.exe](http://biopython.org/DIST/biopython-1.63b.win32-py2.6.exe)
    1,877 Kb -- 32 bit Windows Installer for Python 2.6 and NumPy 1.7
-   [biopython-1.63b.win32-py2.7.exe](http://biopython.org/DIST/biopython-1.63b.win32-py2.7.exe)
    2,003 Kb -- 32 bit Windows Installer for Python 2.7 and NumPy 1.7
-   [biopython-1.63b.win32-py3.3.exe](http://biopython.org/DIST/biopython-1.63b.win32-py3.3.exe)
    2,005 Kb -- 32 bit Windows Installer for Python 3.3 and NumPy 1.7
-   [Tutorial-1.63b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.63b.pdf) -- Documentation

<!-- -->

-   [biopython-1.62.tar.gz](http://biopython.org/DIST/biopython-1.62.tar.gz)
    11,123 Kb -- Source Tarball (28 August 2013)
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
-   [Tutorial-1.62.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.62.pdf) -- Documentation

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
-   [Tutorial-1.62b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.62b.pdf) -- Documentation

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
-   [Tutorial-1.61.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.61.pdf) -- Documentation

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
-   [Tutorial-1.60.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.60.pdf) -- Documentation

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
-   [Tutorial-1.59.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.59.pdf) -- Documentation

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
-   [Tutorial-1.58.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.58.pdf) -- Documentation

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
-   [Tutorial-1.57.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.57.pdf) -- Documentation

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
-   [Tutorial-1.56.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.56.pdf) -- Documentation

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
-   [Tutorial-1.55.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.55.pdf) -- Documentation

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
-   [Tutorial-1.55b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.55b.pdf) -- Documentation

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
-   [Tutorial-1.54.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.54.pdf) -- Documentation

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
-   [Tutorial-1.54b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.54b.pdf) -- Documentation

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
-   [Tutorial-1.53.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.53.pdf) -- Documentation

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
-   [Tutorial-1.52.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.52.pdf) -- Documentation

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
-   [Tutorial-1.51.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.51.pdf) -- Documentation

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
-   [Tutorial-1.51b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.51b.pdf) -- Documentation

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
-   [Tutorial-1.50.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.50.pdf) -- Documentation

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
-   [Tutorial-1.50b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.50b.pdf) -- Documentation

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
-   [Tutorial-1.49.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.49.pdf) -- Documentation

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
-   [Tutorial-1.49b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.49b.pdf) -- Documentation

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
-   [Tutorial-1.48.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.48.pdf) -- Documentation
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
-   [Tutorial-1.47.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.47.pdf) -- Documentation
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
-   [Tutorial-1.46.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.46.pdf) -- Documentation
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
-   [Tutorial-1.45.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.45.pdf) -- Documentation
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
-   [Tutorial-1.44.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.44.pdf) -- Documentation
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
-   [Tutorial-1.43.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.43.pdf) -- Documentation
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
-   [Tutorial-1.42.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.42.pdf) -- Documentation
-   [biopython-1.41.tar.gz](http://biopython.org/DIST/biopython-1.41.tar.gz)
    3,719 Kb (October 28, 2005)
-   [biopython-1.41.zip](http://biopython.org/DIST/biopython-1.41.zip)
    4,241 Kb
-   [biopython-1.41.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.3.exe)
    1,038 Kb
-   [biopython-1.41.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.41.win32-py2.4.exe)
    1,042 Kb
-   [Tutorial-1.41.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.41.pdf) -- Documentation
-   [biopython-1.40b.tar.gz](http://biopython.org/DIST/biopython-1.40b.tar.gz)
    3,437 Kb (February 18, 2005)
-   [biopython-1.40b.zip](http://biopython.org/DIST/biopython-1.40b.zip)
    3,267 Kb
-   [biopython-1.40b.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.3.exe)
    1,019 Kb
-   [biopython-1.40b.win32-py2.4.exe](http://biopython.org/DIST/biopython-1.40b.win32-py2.4.exe)
    1,023 Kb
-   [Tutorial-1.40b.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.40b.pdf) -- Documentation
-   [biopython-1.30.tar.gz](http://biopython.org/DIST/biopython-1.30.tar.gz)
    3,186 Kb (May 14, 2004)
-   [Tutorial-1.30.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.30.pdf) -- Documentation
-   [biopython-1.24.tar.gz](http://biopython.org/DIST/biopython-1.24.tar.gz)
    3,081 Kb (February 16, 2004)
-   [biopython-1.24.zip](http://biopython.org/DIST/biopython-1.24.zip)
    3,623 Kb
-   [biopython-1.24.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.2.exe)
    892 Kb
-   [biopython-1.24.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.24.win32-py2.3.exe)
    894 Kb
-   [Tutorial-1.24.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.24.pdf) -- Documentation
-   [biopython-1.23.tar.gz](http://biopython.org/DIST/biopython-1.23.tar.gz)
    2,241 Kb (October 18, 2003)
-   [biopython-1.23.zip](http://biopython.org/DIST/biopython-1.23.zip)
    2,719 Kb
-   [biopython-1.23.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.2.exe)
    833 Kb
-   [biopython-1.23.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.23.win32-py2.3.exe)
    842 Kb
-   [Tutorial-1.23.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.23.pdf) -- Documentation
-   [biopython-1.22.tar.gz](http://biopython.org/DIST/biopython-1.22.tar.gz)
    2,214 Kb (October 9, 2003)
-   [biopython-1.22.zip](http://biopython.org/DIST/biopython-1.22.zip)
    2,691 Kb
-   [Tutorial-1.22.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.22.pdf) -- Documentation
-   [biopython-1.21.tar.gz](http://biopython.org/DIST/biopython-1.21.tar.gz)
    2,214 Kb
-   [biopython-1.21.zip](http://biopython.org/DIST/biopython-1.21.zip)
    2,897 Kb
-   [biopython-1.21.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.21.win32-py2.2.exe)
    770 Kb
-   [biopython-1.21.win32-py2.3.exe](http://biopython.org/DIST/biopython-1.21.win32-py2.3.exe)
    832 Kb
-   [Tutorial-1.21.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.21.pdf) -- Documentation
-   [biopython-1.20.tar.gz](http://biopython.org/DIST/biopython-1.20.tar.gz)
    2,101 Kb (July 28, 2003)
-   [biopython-1.20.zip](http://biopython.org/DIST/biopython-1.20.zip)
    2,602 Kb
-   [Tutorial-1.20.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.20.pdf) -- Documentation
-   [biopython-1.10.tar.gz](http://biopython.org/DIST/biopython-1.10.tar.gz)
    1,811 Kb (December 17, 2002)
-   [biopython-1.10.zip](http://biopython.org/DIST/biopython-1.10.zip)
    2,300 Kb
-   [biopython-1.10.win32-py2.2.exe](http://biopython.org/DIST/biopython-1.10.win32-py2.2.exe)
    1,199 Kb
-   [Tutorial-1.10.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.10.pdf) -- Documentation
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
-   [Tutorial-1.00a4.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.00a4.pdf) -- Documentation
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
-   [Tutorial-1.00a3.pdf](http://biopython.org/DIST/docs/tutorial/Tutorial-1.00a3.pdf) -- Documentation
