---
title: Development Page for the PopGen Module.
permalink: wiki/PopGen_dev
layout: wiki
---

Introduction
------------

The [`PopGen`](PopGen "wikilink") module contains modules to handle
population genetics data, applications and algorithms.

History and philosophy
----------------------

Most of the existing `Bio.PopGen` features are of non-core population
genetics functionality. This was seen as feature (and not as a bug) in
order to start building a module with functionality where newbie crass
errors would not have dramatic consequences. Currently, with the
experience accumulated is is possible and desirable to concentrate on
core population genetics functionality (i.e., statistics).

Also worth noticing is that we wrap existing functionality whenever
possible. For instance we don't provide our own coalescent simulator,
but we provide wrappers to an existing one which is established and
widely used (SIMCOAL2).

Future Goals
------------

The fundamental goal is to have support for "classic" population
genetics operations (statistics). This should be provided in an
extensible, easy to use and future-proof framework. Code exists (see
below on how to find it), but will probably be refactored. Below there
is also a with list where you can add your desired features.

Code and contributing
---------------------

Your contributions are most welcome. You should follow the general development
guidelines for Biopython.

Wish list
---------

-   support for a binary format - like [HDF5](http://www.pytables.org)
    or this one:
    [snpfile](http://lists.open-bio.org/pipermail/biopython/2008-December/004830.html)
-   support for database: it is frequent to carry analysis on a big
    scale, so it is not unfrequent to use databases to store data
