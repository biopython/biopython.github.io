---
title: PopGen Genepop
permalink: wiki/PopGen_Genepop
layout: wiki
---

Two interfaces are supplied: A general, more complex and more efficient
one (GenePopController) and a simplified, more easy to use, not complete
and not so efficient version (EasyController). EasyController might not
be able to handle very large files, by virtue of its interface.

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
