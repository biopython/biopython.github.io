---
title: Analyzing protein sequences with the ProtParam module.
permalink: wiki/ProtParam
layout: wiki
---

Protein sequences can be analysed by several tools, based on the
**ProtParam** tools on the Expasy Proteomics Server. The module is part of
the `SeqUtils` package.

Protein Sequence Format
-----------------------

The `ProteinAnalysis` class takes one argument, the protein sequence as a
string and builds a sequence object using the `Bio.Seq module`. This is
done just to make sure the sequence is a protein sequence and not
anything else.

Example
-------

``` python
>>> from Bio.SeqUtils.ProtParam import ProteinAnalysis
>>> my_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESV
GEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAIL
FLPLPV"
>>> analysed_seq = ProteinAnalysis(my_seq)
>>> analysed_seq.molecular_weight()
17103.1617
>>> analysed_seq.gravy()
-0.597368421052632
>>> analysed_seq.count_amino_acids()
{'A': 6, 'C': 3, 'E': 12, 'D': 5, 'G': 14, 'F': 6, 'I': 5, 'H': 5, 'K': 12, 'M':
 2, 'L': 18, 'N': 7, 'Q': 6, 'P': 8, 'S': 10, 'R': 6, 'T': 13, 'W': 1, 'V': 5,
 'Y': 8}
 >>>
```

Available Tools
---------------

-  `count_amino_acids`: Simply counts the number times an amino acid is repeated
   in the protein sequence. Returns a dictionary {AminoAcid: Number} and also
   stores the dictionary in `self.amino_acids_content`.
-  `get_amino_acids_percent`: The same as `count_amino_acids`, only returns the
   number in percentage of entire sequence. Returns a dictionary and stores the
   dictionary in `self.amino_acids_content_percent`.
-  `molecular_weight`: Calculates the molecular weight of a protein.
-  `aromaticity`: Calculates the aromaticity value of a protein according to Lobry &
   Gautier (1994, [Nucleic Acids Res., 22, 3174-3180](http://dx.doi.org/10.1093/nar/22.15.3174)).
   It is simply the relative frequency of Phe+Trp+Tyr.
-  `instability_index`: Implementation of the method of Guruprasad *et al.*
   (1990, [Protein Engineering, 4, 155-161](http://dx.doi.org/10.1093/protein/4.2.155)).
   This method tests a protein for stability. Any value above 40 means the protein
   is unstable (=has a short half life).
-  `flexibility`: Implementation of the flexibility method of Vihinen *et al.*
  (1994, [Proteins, 19, 141-149](http://dx.doi.org/10.1002/prot.340190207)).
-  `isoelectric_point`: This method uses the module `IsoelectricPoint` to calculate
   the pI of a protein.
-  `secondary_structure_fraction`: This methods returns a list of the fraction
   of amino acids which tend to be in helix, turn or sheet.
     -  Amino acids in helix: V, I, Y, F, W, L.
     -  Amino acids in turn: N, P, G, S.
     -  Amino acids in sheet: E, M, A, L.

   The list contains 3 values: \[Helix, Turn, Sheet\].

Protein Scales
--------------

`protein_scale(Scale, WindowSize, Edge)`:
The method returns a list of values which can be plotted to view the change
along a protein sequence. You can set several parameters that control the
computation of a *scale* profile, such as the *window size* and the
*window edge* relative weight value. Many scales exist. Just add your favorites
to the `ProtParamData` modules.

-   *Scale*: An amino acid scale is defined by a numerical value assigned
    to each type of amino acid. The most frequently used scales are the
    hydrophobicity or hydrophilicity scales and the secondary structure
    conformational parameters scales, but many other scales exist which are
    based on different chemical and physical properties of the amino acids.
-   *WindowSize*: The window size is the length of the interval to use
    for the profile computation. For a window size n, we use the i - (n - 1)/2
    neighboring residues on each side of residue it compute the score for
    residue i. The score for residue is the sum of the scale values for
    these amino acids, optionally weighted according to their position in
    the window.
-   *Edge* : The central amino acid of the window always has a
    weight of 1. By default, the amino acids at the remaining window
    positions have the same weight, but you can make the residue at the
    center of the window have a larger weight than the others by setting the
    edge value for the residues at the beginning and end of the interval to
    a value between 0 and 1. For instance, for Edge=0.4 and a window size of
    5 the weights will be: 0.4, 0.7, 1.0, 0.7, 0.4.
