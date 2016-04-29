---
title: MotifDev
permalink: wiki/MotifDev
layout: wiki
---

This page describes current state of the development of the Bio.Motif
module

Enhancements currently underway:

-   Expanding the Bio.Motif tutorial on analysis of protein motifs (Dave
    Bridges is workin on this, see
    [http://github.com/davebridges/biopython-biomotif-supplement/tree/master
    his branch on
    github](http://github.com/davebridges/biopython-biomotif-supplement/tree/master_his_branch_on_github "wikilink")
-   writing a simplistic, pure-python, de-novo motif finder
-   writing a wrapper for RSAT tools (http://rsat.ulb.ac.be/rsat/) using
    either local binaries or SOAP

There are also multiple features - suggested mostly by Leighton
Pritchard
[http://lists.open-bio.org/pipermail/biopython-dev/2009-April/005834.html
1](http://lists.open-bio.org/pipermail/biopython-dev/2009-April/005834.html_1 "wikilink")
and
[http://lists.open-bio.org/pipermail/biopython-dev/2009-April/005811.html
2](http://lists.open-bio.org/pipermail/biopython-dev/2009-April/005811.html_2 "wikilink")
- which would extend highly extend usability of the Bio.Motif module:

-   Writing different implementations for different motif types:
    -   consensus sequences (possibly with ambiguity)
    -   alignments (gapped or ungapped-&gt; pwms)
    -   motifs with variable-length gaps
-   writing conversions from alignments to motif instances (possibly
    different for different motifs)
-   integrating with Bio.Seq and Bio.Align, so that the API is more
    natural for the end users
-   making modifications to Seq.startwswith, .endswith .find and .rfind
    to support passing motifs to them. See also
    [this commit](https://github.com/biopython/biopython/commit/49096ecf89d050129ef5b22b66b8bc34fec692ed).

Kind of separate issue (also from Leighton) is to take into account HMM
based motif searching:

-   Supporting HMM-based motifs (with dependencies between positions)
-   particular tools like HMMer

If you want to contribute, please contact bartek@rezolwenta.eu.org or
biopython-dev@biopython.org
