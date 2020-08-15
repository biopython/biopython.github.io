---
title: Genepop Support in Biopython
permalink: wiki/PopGen_Genepop
layout: wiki
redirect_from:
 - wiki/PopGen_dev_Genepop
---

The `Genepop` module allows to access **Genepop** functionality using a Python
interface. This means that the vast majority of **Genepop**'s methods (exact
tests for Hardy–Weinberg equilibrium, population differentiation,
genotypic disequilibrium, F-statistics, null allele frequencies, allele
size-based statistics for microsatellites and much more) can now be
accessed from Python. **Genepop** needs to be installed for this as the code
is only a wrapper.

A parser to **Genepop** files is also available (and documented in the
Tutorial). **Genepop** formatted files can be processed without **Genepop**.

Two interfaces are supplied: A general, more complex and more efficient
one (`GenePopController`) and a simplified, more easy to use, not complete
and not so efficient version (`EasyController`). `EasyController` might not
be able to handle very large files, by virtue of its interface. On the
other hand it provides utility functions to compute some very simple
statistics like allele counts, which are not directly available in the
general interface.

The more complex interface assumes more proficient Python developers
(e.g., by the use of iterators) and for now it is not documented. But
even for experienced Python developers, `EasyController` can be convenient
as long as the required functionality is exposed in `EasyController` and
its performance is deemed acceptable.

For details on the methods used for calculations, check the **Genepop**
documentation, which provides pointers to all papers from where the
calculations are derived.

EasyController tutorial
=======================

Installation
-----------

In order for the controllers to be used, **Genepop** has to be installed in
the system, it can be downloaded from
[here](http://kimura.univ-montp2.fr/~rousset/Genepop.htm).

Before we start, lets test the installation (for this you need a **Genepop**
formated file):

``` python
from Bio.PopGen.GenePop.EasyController import EasyController

ctrl = EasyController(your_file_here)
print ctrl.get_basic_info()
```

Replace `your_file_here` with the name and path to your file. If you get
a `IOError: Genepop not found` then Biopython cannot find your **Genepop**
executable. If **Genepop** is not on the PATH, you can add it to the
constructor line, i.e.

``` python
ctrl = EasyController(your_file_here, path_to_genepop_here)
```

If everything is working, now we can go on and use **Genepop**. For the
examples below, we will use the **Genepop** file
[`big.gen`](https://raw.githubusercontent.com/biopython/biopython/master/Tests/PopGen/big.gen)
made available with the unit tests. We will also assume that there is a
`ctrl` object initialized with the relevant file chosen.

We start by getting some basic info

``` python
pop_names, loci_names = ctrl.get_basic_info()
```

Returns the list of population names and loci names available on the
file.

**Caveat:** Most existing **Genepop** files provide erroneous data regarding
population names. In many cases that information might not be trusted.
Assessing population information is, most of the times, done by the
relative position of the population in the file, not the name. So the
first population is the file is index 0, the second index 1, and so
on...

Statistics
----------

### Heterozygosity

Lets get heterozygosity info for a certain population and a certain
allele:

``` python
(exp_homo, obs_homo, exp_hetero, obs_hetero) = ctrl.get_heterozygosity_info(0, "Locus2")
```

Will get expected and observed homozygosity and heterozygosity for
population 0 and Locus2 (of the file big.gen, if you are using another
file, adjust the population position and locus name accordingly).

### Existing alleles

It is possible to get the list of all alleles of a certain locus in a
certain population:

``` python
allele_list = ctrl.get_alleles(0, "Locus2")
```

*allele\_list* will be \[3, 20\] (i.e., alleles 3 and 20 are on the
population).

The number of alleles is simply getting with `len(allele_list)`.

It is also possible to get the list of all alleles of a certain locus
for all populations:

``` python
all_allele_list = ctrl.get_alleles_all_pops("Locus2")
```

*all\_allele\_list* will be \[3, 20\].

### Allele and genotype frequencies

It is possible to get the frequency of alleles in a certain population:

``` python
allele_data = ctrl.get_allele_frequency(0, "Locus2")
```

*allele\_data* will be (62, {3: 0.88700000000000001, 20: 0.113}). That is
there are 62 genes. 88.7% are 3 and 11.3% are 20.

We can get similar information for genotypes (diploid data). Expected
frequencies will also be reported:

``` python
genotype_list = ctrl.get_genotype_frequency(0, "Locus2")
```

*genotype\_list* will be: \[(3, 3, 24, 24.3443), (20, 3, 7,
6.3114999999999997), (20, 20, 0, 0.34429999999999999)\]

Lets interpret the first element: There are 24 individuals which have a
genotype of (3, 3), whereas the expected number of individuals with that
genotype is 24.2443.

### F statistics

Lets start with general multilocus F statistics:

``` python
Fis, Fst, Fit = ctrl.get_multilocus_f_stats()
```

This gets multilocus Fis, Fst and Fit.

Lets get that (and a bit more) per locus:

``` python
Fis, Fst, Fit, Qintra, Qinter = ctrl.get_f_stats("Locus2")
```

This gets single locus Fis, Fst and Fit, Qintra and Qinter.

There are specific sections below for Fst and Fis (where pairwise and
population specific variants are introduced). On the Fis section Qintra
and Qinter are explained.

### Fst

Lets get the pairwise Fst for a certain locus:

``` python
pair_fst = ctrl.get_avg_fst_pair_locus("Locus4")
```

Will return a map where the key is the pair composed of *population1*,
*population2* (the population Id). *population2* is always LOWER than
*population1*. Example: the pairwise Fst for Locus4 between population 0
and population 3 is given by `pair_fst[(3,0)]`.

You can also get the multilocus pairwise Fst estimate:

``` python
multilocus_fst = ctrl.get_avg_fst_pair()
```

This will return the same data structure as above but with a multilocus
pairwise Fst.

### Fis

We will now get the Fis of a certain locus/population plus a few other
statistics:

``` python
allele_dict, summary_fis = ctrl.get_fis(0, "Locus2")
```

Lets have a detailed look at the output of `get_fis`:

``` python
summary_fis = (62, -0.1111, -0.11269999999999999)

allele_dict = {3: (55, 0.8871, -0.1111), 20: (7, 0.1129, -0.1111)}
```

*summary\_fis* holds a triple with: total number of alleles, Cockerham and
Weir Fis, Robertson and Hill Fis.

*allele\_dict* holds for each allele (being each allele the key), number
of repetitions of the allele, frequency and Cockerham and Weir Fis.

So, from the above results the following can be read: there are 62 genes
with 2 different alleles (55 are of type 3, and 7 of type 20). Type 3 has
frequency 0.89 and type 20 of 0.11. All CW Fis are -0.111 and the RH Fis is
-0.112.

Lets now get multilocus Fis:

``` python
pop_list = ctrl.get_avg_fis()
```

*pop\_list* will return an element per population. Each element is a
quadruple containing:

1.  population name (again, population names are not to be trusted)
2.  1 - QIntra: Gene diversity between individuals
3.  1 - QInter: Gene diversity among individuals within populations
4.  Fis

### Migration

We can get an estimation of the number of migrants:

``` python
samp_size, priv_allele_freq, mig10, mig25, mig50, migcorr = ctrl.estimate_nm()
```

*samp\_size* is mean sample size, *priv\_allele\_freq* is the mean frequency
of private alleles, *mig10* is the number of migrants for Ne=10, *mig25* for
Ne=25, *mig50* for Ne=50 and *migcorr* is the number of migrants after
correcting for expected size.

Tests
-----

Tests are normally computationally intensive as they are normally based
on a Markov Chain algorithm. In some cases full enumeration approaches
are available but those can only be applied for locus with a very low
number of alleles. **This means that most tests will take quite some
time to complete**.

For more details about Markov Chain parameters below (dememorization,
batched and iterations) please consult the **Genepop** manual. Also consult
the manual to understand when full enumeration is applicable.

### Hardy-Weinberg equilibrium

Lets start by testing Hardy-Weinberg equilibrium for each loci in each
population:

``` python
loci_map = ctrl.test_hw_pop(1, "excess")
```

The second parameter can be *probability*, *excess* or *deficiency*.
*probability* is the standard Haldane HW test. Use *deficiency* when you
are interested in heterozygote deficiency or *excess* if you are
interested in excess.

The output is a map where the key is the locus name. The content is a
tuple containing P-value, Standard Error, Fis (Weir and Cockerham), Fis
(Robertson and Hill) and steps.

``` python
pop_test, loc_test, all_test = ctrl.test_hw_global("deficiency")
```

Use *deficiency* when you are interested in heterozygote deficiency or
*excess* if you are interested in excess. *probability* does not apply
here like in test\_hw\_pop.

The output is a triple:

1.  *pop\_test* is a list with an element per population including
    P-value, Standard error and switches.
2.  *loc\_test* is the same list, but with one element per locus including
    locus name, P-value, Standard error and switches.
3.  *all\_test* are the overall results consisting of a triple P-value,
    standard error and switches.

### Linkage Disequilibrium

We can test if 2 loci are in linkage disequilibrium using the log
likelihood ratio statistic (G-test).

``` python
chi2, df, pval = ctrl.test_ld_all_pair(
    "Locus1", "Locus2", dememorization=1000, batches=10, iterations=100
)
```

Returns the Chi square value, degrees of freedom and P value for the G
statistic.

Isolation By Distance (IBD)
---------------------------

Isolation By Distance (IBD) analysis **requires** a special form of
**Genepop** files:

1.  One individual per population
2.  The name of the individual has to be its coordinates

Example:

```
...
Pop
0 15, 0201  0303 0102 0302 1011
Pop
0 30, 0202  0301 0102 0303 1111
Pop
0 45, 0102  0401 0202 0102 1010
Pop
0 60, 0103  0202 0101 0202 1011
Pop
0 75, 0203  0204 0101 0102 1010
POP
15 15, 0102 0202 0201 0405 0807
...
```

Note that the example file that we are using, **cannot be used for this
case**.

There is a single call for IBD analysis:

``` python
estimate, distance, (a, b), (bb, bblow, bbhigh) = ctrl.calc_ibd(
    self, is_diplo=True, stat="a", scale="Log", min_dist=0.00001
)
```

-  `is_diplo` specifies if data is diploid (`True`) or haploid (`False`).
-  `stat` is either `'a'` or `'e'` (see the **Genepop** manual for details).
-  `scale` is either `'Log'` or `'Linear'`. `'Log'` is used for 2D coordinates and
   `'Linear'` for 1D.

Only pairwise comparisons above `min_dist` are used to compute regression
coefficients.

The method returns:

-  *estimate*, a triangular matrix containing genetic distances among samples
   according to the chosen statistic.
-  *distance*, a triangular matrix containing distances (log or linear) among
   smamples.
-  *a* and *b *are the parameter fits for the regression. *bblow* and *bbhigh* are
   the bootstrap confidence intervals for the *b* parameter (*bb* should be
   very close to *b*).

Interpretation of the triangular matrices should be done like this:
Pythonwise, a matrix is implemented with a list of lists of numbers,
like this:

``` python
[[0.1], [0.2, 0.3], [0.4, 0.5, 0.6]]
```

The above data structure corresponds to the following triangular matrix:

```
      1    2    3
2   0.1
3   0.2  0.3
4   0.4  0.5  0.6
```
