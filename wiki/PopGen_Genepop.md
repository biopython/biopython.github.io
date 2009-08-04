---
title: PopGen Genepop
permalink: wiki/PopGen_Genepop
layout: wiki
---

Two interfaces are supplied: A general, more complex and more efficient
one (GenePopController) and a simplified, more easy to use, not complete
and not so efficient version (EasyController). EasyController might not
be able to handle very large files, by virtue of its interface, on the
other hand it provides utility functions to compute some very simple
statistics like allele counts, which are not directly available in the
general interface.

The more complex interface assumes more proficient Python developers
(e.g., by the use of iterators) and for now it is not documented. But
even for experienced Python developers, EasyController can be convenient
as long as the required functionality is exposed in EasyController and
its performance is deemed acceptable.

In order for the controllers to be used, Genepop has to be installed in
the system, it can be dowloaded from
[here](http://kimura.univ-montp2.fr/~rousset/Genepop.htm).

EasyController tutorial
-----------------------

Before we start, lets test the installation (for this you need a genepop
formated file):

``` python
from Bio.PopGen.GenePop.EasyController import EasyController

ctrl = EasyController(your_file_here)
print ctrl.get_basic_info()
```

Replace your\_file\_here with the name and path to your file. If you get
a **IOError: Genepop not found** then Biopython cannot find your Genepop
executable. If Genepop is not on the PATH, you can add it to the
constructor line, i.e.

``` python
ctrl = EasyController(your_file_here, path_to_genepop_here)
```

If everything is working, now we can go on and use Genepop. For the
examples below, we will use the genepop file
[big.gen](http://biopython.open-bio.org/SRC/biopython/Tests/PopGen/big.gen)
made available with the unit tests. We will also assume that there is a
ctrl object initialized with the relevant file chosen.

We start by getting some basic info

``` python
pop_names, loci_names = ctrl.get_basic_info()
```

Returns the list of population names and loci names available on the
file.

**Caveat:** Most existing Genepop files provide erroneous data regarding
population names. In many cases that information might not be trusted.
Assessing population information is, most of the times, done by the
relative position of the population in the file, not the name. So the
first population is the file is index 0, the second index 1, and so
on...

Lets get heterozygosity info for a certain population and a certain
allele:

``` python
(exp_homo, obs_homo, exp_hetero, obs_hetero) = ctrl.get_heterozygosity_info(0,"Locus2")
```

Will get expected and observed homozygosity and heterozygosity for
population 0 and Locus2 (of the file big.gen, if you are using another
file, adjust the population position and locus name accordingly).

It is possible to get the list of all alleles of a certain locus in a
certain population:

``` python
allele_list = ctrl.get_alleles(0,"Locus2")
```

The statistic allele count is simply getting len(allele\_list).

It is also possible to get the list of all alleles of a certain locus
for all populations:

``` python
all_allele_list = ctrl.get_alleles_all_pops("Locus2")
```

We will now get the Fis of a certain locus/population plus a few other
statistics:

``` python
allele_dict, summary_fis =  ctrl.get_fis(0,"Locus2")
```

Lets have a detailed look the output of get\_fis:

``` python
summary_fis = (62, -0.1111, -0.11269999999999999)

allele_dict = {
    3: (55, 0.8871, -0.1111),
    20: (7, 0.1129, -0.1111)
}
```

summary\_fis holds a triple with: total number of alleles, Cockerham and
Weir Fis, Robertson and Hill Fis.

allele\_dict holds for each allele (being each allele the key), number
of repetitions of the allele, frequency and Cockerham and Weir Fis.

So, from the above results the following can be read: there are 62 genes
with 2 different allees (55 are of type 3, and 7 of type 20). 3 has
frequency 0.89 and 20 0.11. All CW Fis are -0.111 and the RH Fis is
-0.112.

If the objective is just to get allele frequencies, then there is an
easier way:

``` python
total_alleles, count = ctrl.get_allele_frequency(0,"Locus2")
```

total\_alleles returns the number of total alleles and count is a
dictionary whose key is the allele id and the value is the frequency.

``` python
print ctrl.estimate_nm()
```

``` python
print ctrl.get_avg_fst_pair_locus("Locus4")
```

``` python
print ctrl.get_avg_fst_pair()
```

``` python
print ctrl.get_avg_fis()
```

``` python
print ctrl.test_hw(1, "excess")
```

``` python
print ctrl.get_multilocus_f_stats()
```

``` python
print ctrl.get_f_stats("Locus2")
```

### Tests

``` python
print ctrl.test_ld_all_pair("Locus1", "Locus2",
    dememorization=1000, batches=10, iterations=100)
```

### Isolation By Distance (IBD)
