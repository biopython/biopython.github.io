---
title: GSOC2011 MocapyExt
permalink: wiki/GSOC2011_MocapyExt
layout: wiki
---

Introduction
------------

BioPython is a very popular library in Bioinformatics and Computational
Biology. Mocapy++ is a machine learning toolkit for training and using
Bayesian networks. However, Mocapy++ is implemented in C++ and its
monolithic architecture does not provide any mechanism to plug-in new
node types (probability distributions/densities) without prior source
code conversion to C++ and the recompilation of the library. The goal of
this project is to develop an easy-to-use plug-in system for Mocapy++,
which would allow to load and test probability distributions on the fly.
If a user is working in C++, the user could load/embed arbitrary
probability distributions in Mocapy++ and then proceed to use them in a
familiar manner.

Author & Mentors
----------------

[Justinas V. Daugmaudis](User%3AJustinas_Daugmaudis "wikilink")
vygis.d@gmail.com

**Mentors**

  
Thomas Hamelryck

Eric Talevich

Work Plan
---------

Project Schedule
----------------
