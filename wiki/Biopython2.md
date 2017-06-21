---
title: Biopython2
permalink: wiki/Biopython2
layout: wiki
---

Ideas for Biopython 2
---------------------

This page is mostly a collection of ideas from discussions on the mailling list regarding a possible version 2 of Biopython.

If I missed something, please accept my apologies and change accordingly.

I am trying to keep a neutral tone, but sometime we might want to include major arguments in favor and against some option
so that readers are made aware of the gist of main arguments

# Dependencies

  - Depend on matplotlib, numpy and scipy or
  - Create a list of acceptable dependencies (e.g. requests or pandas) that developers could use

# Documentation

  - There seems to be a consensus that API documentation takes precedence
  - numpydoc and Sphinx proposed
  - Juptyer notebook vs HTML proposed for tutorials.
  - Latex to be discontinued?

# Deprecate Python 2?

Ongoing discussion

# Which Python 3 version?

   - All of them? (no one defended this)
   - 3.5
   - Most recent version when we start developing Bioython 2 (e.g. now would mean 3.6)


# Lower case package and module naming

Lower case for packages and modules seems consensual.

# Top-level policies

 - Name this biopy.*? biopython.*? biop.*?
 - How many modules on the top-level? None or some basic stuff (life exceptions, abstract file management)
 - Automatic import of everything, or just partial? or even nothing?


# Package organization

Which sub-packages? Tightly related to the modular approach (see below)

Anyone wants to write a bit here about options?


# Module system

 - Should there be a module system with core modules allowing third parties to to make extensions, a la [Biogems](
http://biogems.info/)
 - Bow's work might be a [starting point](https://github.com/bow/poc_biopy)


# Removal of outdated sub-packages

Discuss which packages to remove


# Relationship between Biopython 1 and 2

 - Is there any code sharing? If so how? Code sharing should jeopardize desired Bioypthon 2 features and architecture
 - APIs would probably not be compatible

# Feasability

Is there even time to do this?
