---
title: Interfacing with "Phylogenetic Analysis by Maximum Likelihood" (PAML) package.
permalink: wiki/PAML
layout: wiki
tags:
 - Wiki Documentation
---

This module provides an interface to the
[PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (Phylogenetic
Analysis by Maximum Likelihood) package of programs. Currently,
interfaces to the programs *codeml*, *baseml* and *yn00* as well as a
Python re-implementation of *chi2* have been included.

Availability
------------

This module is available in Biopython 1.58 and later.

To use this module, you must have PAML version 4.1 or greater installed
on your system.

codeml
------

The *codeml* interface is provided in the following module:

``` python
from Bio.Phylo.PAML import codeml
```

### The Codeml object

#### Initialization

The interface is implemented as an object which maintains program
options. In order to run *codeml*, typically one provides a control file
which indicates the locations of an alignment file, a tree file and an
output file, along with a series of options dictating how the software
is to be run. In the `Codeml` object, the three file locations are stored
as string attributes. The three file locations, as well as the location
of the desired working directory, may be specified in the `Codeml`
constructor as follows:

``` python
cml = codeml.Codeml(alignment = "align.phylip", tree = "species.tree",
                    out_file = "results.out", working_dir = "./scratch")
```

They may also be set individually:

``` python
cml = codeml.Codeml()
cml.alignment = "align.phylip"
cml.tree = "species.tree"
cml.out_file = "results.out"
cml.working_dir = "./scratch"
```

Note that all file locations are converted to locations relative to the
working directory. PAML programs have a limit on the lengths of file
location strings; converting to relative locations will shorten these
strings, allowing you generally to sidestep this limitation.

#### Options

The *codeml* runtime options are stored in a dictionary object that is
keyed by the option names. For more information on the usage of these
options, please refer to the PAML user manual.

The complete set of options and their current values may be printed:

``` python
>>> cml.print_options()
verbose = None
CodonFreq = None
cleandata = None
fix_blength = None
NSsites = None
fix_omega = None
clock = None
ncatG = None
runmode = None
fix_kappa = None
fix_alpha = None
Small_Diff = None
method = None
Malpha = None
aaDist = None
RateAncestor = None
aaRatefile = None
icode = None
alpha = None
seqtype = None
omega = None
getSE = None
noisy = None
Mgene = None
kappa = None
model = None
ndata = None
```

Setting an option to a value of `None` will cause it to be ignored by
*codeml* (it will be omitted from the final control file). Options may
be set by the `set_option()` function and their values may be retrieved
by the `get_option()` function:

``` python
>>> cml.set_options(clock=1)
>>> cml.set_options(NSsites=[0, 1, 2])
>>> cml.set_options(aaRatefile="wag.dat")
>>> cml.get_option("NSsites")
[0, 1, 2]
```

The `NSsites` option deserves special attention: in a *codeml* control
file, NSsites models are entered as a space-delimited list of numbers,
such as 0 1 2, whereas in the `Codeml` object it is stored as a Python
list.

Finally, options may be read in from an existing *codeml* control file
or they may be written to a new file. Writing to a file isn't necessary,
as this is done automatically when executing the `run()` method (see
below). The control file to read is provided as an argument to the
`read_ctl_file()` method, while the `write_ctl_file()` method writes to
the `Codeml` object's `ctl_file` attribute

``` python
>>> cml.read_ctl_file("codeml.ctl")
>>> cml.print_options()
verbose = 1
CodonFreq = 2
cleandata = 1
fix_blength = None
NSsites = 0
fix_omega = 0
clock = 0
ncatG = 8
runmode = 0
fix_kappa = 0
fix_alpha = 1
Small_Diff = 5e-07
method = 0
Malpha = 0
aaDist = 0
RateAncestor = 1
aaRatefile = dat/jones.dat
icode = 0
alpha = 0.0
seqtype = 2
omega = 0.4
getSE = 0
noisy = 9
Mgene = 0
kappa = 2
model = 2
ndata = None

>>> cml.ctl_file = "control2.ctl"
>>> cml.write_ctl_file()
```

#### Running the program

Executing the object's `run()` method will run *codeml* with the current
options and will return the parsed results in a dictionary object (see
below). You can also pass a number of optional arguments to the `run()`
method:

-   `verbose` (boolean): *codeml*'s screen output is suppresed by default;
    set this argument to `True` to see all of the output as it
    is generated. This is useful if *codeml* is failing and you need to
    see its error messages.
-   `parse` (boolean): set to `False` to skip parsing the results. `run()`
    will instead return `None`
-   `ctl_file` (string): provide a path to an existing control file
    to execute. The file is not parsed and read into `Codeml`'s
    options dictionary. If set to `None` (default), the options dictionary
    is written to a control file, which is then used by *codeml*.
-   `command` (string): provide a path to the *codeml* executable. This is
    set to "codeml" by default, so if the program is in your system
    path, that should suffice. If it's not in your system path or, for
    example, if you use multiple versions of PAML, you may instead
    provide the full path to the executable.

If the *codeml* process exits with an error, a `PamlError` exception will
be raised.

### read() and Results

As previously stated, the `Codeml.run()` method returns a dictionary
containing results parsed from *codeml*'s main output file.
Alternatively, an existing output file may be parsed using the `read()`
function:

``` python
results = codeml.read()
```

The results dictionary is organized hierarchically and is (in most
cases) keyed in accordance with the terminology used in the output file.
Numeric values are automatically converted to numeric Python types. As
the output of *codeml* varies widely depending on the parameters used,
the contents of the results dictionary will vary accordingly. Similarly,
in the case of a runtime error, the results dictionary may be empty or
missing contents. Thus, it is advisable to use Python's `dict.get(key)`
method rather than `dict[key]` in order to handle missing elements
gracefully.

All possible keys of the results dictionary are organized as follows:

-   **version** : the program version number that generated the
    results
-   **model** : the model description
-   **codon model** : the codon frequency model used
-   **site-class model** : the site-class model used, in the case of a
    single model
-   **lnL max** : the maximum log-likelihood
-   **NSsites** : NSsites results
    -   **model\_number** (integer) : the NSsites model number
        -   **description** : description of the NSsites model
        -   **lnL** : log-likelihood
        -   **parameter list** : string list of model parameters
            (found beneath the "lnL" line in the output file)
        -   **SEs** : string list of standard error values
        -   **tree length** : tree length
        -   **dS tree** : tree with dS estimates as branch labels
        -   **dN tree** : tree with dN estimates as branch labels
        -   **omega tree** : tree with dN/dS estimates as branch
            labels
        -   **parameters** : model parameters
            -   **rates** : rate estimates
            -   **kappa** : kappa estimate
            -   **omega** : dN/dS estimate(s)
            -   **genes** : results for genes in files containing
                multiple alignments
                -   **gene number** (integer)
                    -   **kappa** : gene kappa estimate
                    -   **omega** : gene omega estimate
            -   **dN** : dN estimate
            -   **dS** : dS estimate
            -   **site classes** : results for multiple site classes
                -   **site class number** (integer)
                    -   **proportion** : proportion of sites belonging
                        to this class
                    -   **omega** : omega estimate for this class
                    -   **branch types** : results for branch types
                        under Clade Model C or Branch Site A
                        -   **branch type number** (integer) : omega
                            estimate for this branch type (Clade
                            Model C)
                        -   **foreground** : omega estimate for the
                            foreground branch (Branch Site A)
                        -   **background** : omega estimate for the
                            background branch (Branch Site A)
            -   **branches** : branch-specific results
                -   **branch** : branch description (ie "6..7")
                    -   **t** : t estimate for this branch
                    -   **N** : N estimate for this branch
                    -   **S** : S estimate for this branch
                    -   **omega** : omega estimate for this branch
                    -   **dN** : dN estimate for this branch
                    -   **dS** : dS estimate for this branch
                    -   **N\*dN** : N\*dN for this branch
                    -   **S\*dS** : S\*dS for this branch
-   **pairwise** : results for pairwise comparisons
    -   **sequence1** : the first sequence's name
        -   **sequence2** : the second sequence's name (note that order
            doesn't matter; seq1 vs seq2 and seq2 vs seq1 are
            both stored)
            -   **t** : t estimate for this comparison
            -   **S** : S estimate for this comparison
            -   **N** : N estimate for this comparison
            -   **omega** : dN/dS estimate for this comparison
            -   **dN** : dN estimate for this comparison
            -   **dS** : dS estimate for this comparison
-   **distances** : amino acid distances
    -   **raw** : raw amino acid distances
        -   **sequence1** : the first sequence's name
        -   **sequence2** : the distance to the second sequence's amino
            acid sequence
    -   **ml** : maximum likelihood amino acid distances
        -   **sequence1** : the first sequence's name
        -   **sequence2** : the distance to the second sequence's amino
            acid sequence

baseml & yn00
-------------

Interfaces to *baseml* and *yn00* are provided in the following modules:

``` python
from Bio.Phylo.PAML import baseml, yn00

bml = baseml.Baseml()
yn = yn00.Yn00()
```

`Baseml` and `Yn00` share the same methods and attributes as `Codeml`, and are
thus used in the same manner. It should be noted, however, that `Yn00`
does not have a tree attribute, as *yn00* does not require a tree file.

### Baseml

The parameters available in the Baseml options are:

``` python
>>> bml.print_options()
verbose = None
cleandata = None
fix_blength = None
nparK = None
model_options = None
clock = None
ncatG = None
runmode = None
fix_kappa = None
fix_alpha = None
noisy = None
method = None
fix_rho = None
Malpha = None
nhomo = None
RateAncestor = None
icode = None
rho = None
alpha = None
getSE = None
Small_Diff = None
Mgene = None
kappa = None
model = None
ndata = None
```

All of the possible contents of the `Baseml` results file are:

-   **version** : the version number of Baseml used to generate the
    results
-   **lnL max** : the maximum log-likelihood
-   **lnL** : the log-likelihood
-   **tree length** : the length of the estimated tree
-   **tree** : the estimated tree
    -   **parameters** : model parameters
        -   **parameter list** : string list of model parameters
            (found beneath the "lnL" line in the output file)
        -   **SEs** : string list of standard error values
        -   **kappa** : kappa estimate(s)
        -   **branches** : branch-specific estimates
            -   **branch number** (integer)
                -   **t** : t estimate for this branch
                -   **kappa** : kappa estimate for this branch
                -   **TS** : TS estimate for this branch
                -   **TV** : TV estimate for this branch
        -   **rate parameters** : list of rate parameters
        -   **rates** : list of rates
        -   **Q matrix** : Q matrix
            -   **matrix** : the actual matrix stored as a list of
                lists (4x4)
            -   **average Ts/Tv** : average Ts/Tv for the matrix
        -   **alpha** : alpha estimate
        -   **base frequencies** : dictionary of base frequencies
            -   **base** ("A", "T", "C" or "G") : frequency of base
        -   **rate frequencies** : list of rate frequencies
        -   **nodes** : branch-specific frequency parameters
            -   **node number** (integer)
                -   **root** : boolean indicating if the node is the
                    root
                -   **frequency parameters** : list of frequency
                    parameters for this node
                -   **base frequencies** : dictionary of base
                    frequencies
                    -   **base** ("A", "T", "C" or "G") : frequency of
                        base

### Yn00

The parameters available in the `Yn00` options are:

``` python
>>> yn.print_options()
commonf3x4 = None
weighting = None
icode = None
ndata = None
verbose = None
```

All of the possible contents of the `Yn00` results file are:

-   **sequence1** : first sequence name
    -   **sequence2** : second sequence name (order does not matter;
        both seq1 vs seq2 and seq2 vs seq1 are stored)
        -   **NG86** : results for the Nei & Gojobori (1986) method
            -   **omega** : dN/dS estimate
            -   **dN** : dN estimate
            -   **dS** : dS estimate
        -   **YN00** : results for the Yang & Nielsen (2000) method
            -   **S** : S estimate
            -   **N** : N estimate
            -   **t** : t estimate
            -   **kappa** : kappa estimate
            -   **omega** : dN/dS estimate
            -   **dN** : dN estimate
            -   **dN SE**: dN standard error
            -   **dS** : dS estimate
            -   **dS SE** : dS standard error
        -   **LWL85** : results for the LWL85 method
            -   **dS** : dS estimate
            -   **dN** : dN estimate
            -   **omega** : dN/dS estimate
            -   **S** : S estimate
            -   **N** : N estimate
        -   **LWL85m** : results for the LWL85m method
            -   **dS** : dS estimate
            -   **dN** : dN estimate
            -   **omega** : dN/dS estimate
            -   **S** : S estimate
            -   **N** : N estimate
            -   **rho** : rho estimate
        -   **LPB93** : results for the LPB93 method
            -   **dS** : dS estimate
            -   **dN** : dN estimate
            -   **omega** : dN/dS estimate

chi2
----

The `chi2` module offers an easy method to retrieve *p*-values from a
Chi-squared cumulative distribution function for likelihood ratio tests,
which are performed frequently when using PAML programs. As of the
current version of PAML, the *chi2* program does not allow passing both
a test statistic and the degrees of freedom as command-line arguments.
Thus, the `chi2` module currently consists of a re-implementation of
Ziheng Yang's original code in pure Python. For most cases, this should
not affect you, however using very large numbers of degrees of freedom,
such as in the FMutSel vs FMutSel0 model test in *codeml* (41 degrees of
freedom), this Python version runs significantly slower than the
original.

To retrieve a *p*-value, simply import the module

``` python
 from Bio.Phylo.PAML.chi2 import cdf_chi2
```

And use the `cdf_chi2()` function:

``` python
>>> df = 2
>>> statistic = 7.21
>>> cdf_chi2(df, statistic)
0.027187444813595696
```
