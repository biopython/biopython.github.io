---
title: Packages
permalink: wiki/Packages
layout: default
---

As per our [Downloads Page](Download "wikilink") page, we generally
recommend using Python's package manager ``pip`` to install Biopython:

``` bash
pip install biopython
```

However, this is not the only option and a separate packaging system may
be more appropriate for your system.

Conda
=====

If your Python is installed using [conda](https://conda.io/docs/), for
example using [miniconda](https://conda.io/miniconda.html) or
[anaconda](https://www.continuum.io/what-is-anaconda), then you should
be able to use Biopython from the conda packages:

``` bash
conda install -c conda-forge biopython
```

or:

``` bash
conda update -c conda-forge biopython
```

We deliberately recommend using [Biopython from the conda-forge
channel](https://anaconda.org/conda-forge/biopython), as this is usually
up to date and covers Windows, Mac OS X and Linux. The default Conda
channel does have Biopython, but is often out of date.

Note Conda is available on Windows, Mac OS X and Linux, and covers far
more than just Python.


Linux Packages
==============

Although we would generally recommend ``pip``, most Linux systems will
have a Biopython package available - it may however by a little out of
date.

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

Biopython is an official Fedora package (since Fedora 5). The package is named
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

Other Unix-like systems
=======================

### FreeBSD and Ports

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

This should automatically fetch and install Biopython (as well as its necessary
dependencies).
