---
title: GSOC2011 Mocapy
permalink: wiki/GSOC2011_Mocapy
layout: wiki
---

Mocapy++ is a machine learning toolkit for training and using Bayesian
networks. It has been used to develop probabilistic models of
biomolecular structures. The goal of this project is to develop a Python
interface to Mocapy++ and integrate it with Biopython. This will allow
the training of a probabilistic model using data extracted from a
database. The integration of Mocapy++ with Biopython will provide a
strong support for the field of protein structure prediction, design and
simulation.

Introduction
------------

Discovering the structure of biomolecules is one of the biggest problems
in biology. Given an amino acid or base sequence, what is the three
dimensional structure? One approach to biomolecular structure prediction
is the construction of probabilistic models. A Bayesian network is a
probabilistic model composed of a set of variables and their joint
probability distribution, represented as a directed acyclic graph. A
dynamic Bayesian network is a Bayesian network that represents sequences
of variables. These sequences can be time-series or sequences of
symbols, such as protein sequences. Directional statistics is concerned
mainly with observations which are unit vectors in the plane or in
three-dimensional space. The sample space is typically a circle or a
sphere. There must be special directional methods which take into
account the structure of the sample spaces. The union of graphical
models and directional statistics allows the development of
probabilistic models of biomolecular structures. Through the use of
dynamic Bayesian networks with directional output it becomes possible to
construct a joint probability distribution over sequence and structure.
Biomolecular structures can be represented in a geometrically natural,
continuous space. Mocapy++ is an open source toolkit for inference and
learning using dynamic Bayesian networks that provides support for
directional statistics. Mocapy++ is excellent for constructing
probabilistic models of biomolecular structures; it has been used to
develop models of protein and RNA structure in atomic detail. Mocapy++
is used in several high-impact publications, and will form the core of
the molecular modeling package Phaistos, which will be released soon.
The goal of this project is to develop a highly useful Python interface
to Mocapy++, and to integrate that interface with the Biopython project.
Through the Bio.PDB module, Biopython provides excellent functionality
for data mining biomolecular structure databases. Integrating Mocapy++
and Biopython will allow training a probabilistic model using data
extracted from a database. Integrating Mocapy++ with Biopython will
create a powerful toolkit for researchers to quickly implement and test
new ideas, try a variety of approaches and refine their methods. It will
provide strong support for the field of biomolecular structure
prediction, design, and simulation.

Author & Mentors
----------------

[Michele Silva](https://github.com/mchelem) (michele.silva@gmail.com)

**Mentors**


Thomas Hamelryck

Eric Talevich

Project Schedule
----------------

### Work Plan

'''Gain understanding of SEM and directional statistics '''

-   Review the theory behind machine learning for bioinformatics, Markov
    chain Monte Carlo and dynamic Bayesian networks.

<!-- -->

-   Build the theoretical background on the algorithms used in Mocapy++,
    such as parameter learning of Bayesian networks using Stochastic
    Expectation Maximization (SEM).

'''Study Mocapy++'s use cases '''

-   Read several papers and attempt to replicate part of the experiments
    described using Mocapy++.

<!-- -->

-   Get a better understanding of biological sequence analysis done
    through probabilistic models of proteins and nucleic acids.

'''Work with Mocapy++ '''

-   Understand Mocapy++'s internal architecture and algorithms by
    exploring its source code and running its test cases.

<!-- -->

-   Research other applications of Mocapy++ in Bioinformatics.

'''Design Mocapy++'s Python interface '''

-   Explore the source code of Biopython to understand its design
    and implementation. The Mocapy++ interface to be included in
    Biopython must be made compatible with the methods of solving
    problems in Biopython.

<!-- -->

-   Design a Python interface for Mocapy++, based on its data structures
    and algorithms. Examine Mocapy++'s use cases and existing test cases
    to provide guidance for the interface design.

'''Implement Python bindings '''

-   Implement test cases in Python using the new interface to Mocapy++.

<!-- -->

-   Implement python bindings for the defined interface.

'''Explore Mocapy++'s applications '''

-   Develop example applications that involve data mining of
    biomolecular structure databases using Biopython.

<!-- -->

-   Formulate probabilistic models using Python-Mocapy++. Apply the
    models to solve biological problems. Examples of problems that can
    be solved using dynamic Bayesian networks include deciding if a pair
    of sequences is evolutionarily related, finding sequences which are
    homologous to a known evolutionary family and predicting RNA
    secondary structure.

Project Code
------------

Hosted at [the gSoC11 Mocapy
branch](http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/)

Project Progress
----------------

### Options to create Python bindings to C++ code

#### Swig

There is already an effort to provide bindings for Mocapy++ using Swig.
However, Swig is not the best option if performance is to be required.
The Sage project aims at providing an open source alternative to
Mathematica or Maple. Cython was developed in conjunction with Sage (it
is an independent project, though), thus it is based on Sage's
requirements. They tried Swig, but declined it for performance issues.
According to the [Sage programming
guide](http://sage.math.washington.edu/tmp/sage-2.8.12.alpha0/doc/prog/node36.html)
"The idea was to write code in C++ for SAGE that needed to be fast, then
wrap it in SWIG. This ground to a halt, because the result was not
sufficiently fast. First, there is overhead when writing code in C++ in
the first place. Second, SWIG generates several layers of code between
Python and the code that does the actual work". This was written back in
2004, but it seems things didn't evolve much. The only reason I would
consider Swig is for future including Mocapy++ bindings on BioJava and
BioRuby projects.

#### Boost Python

Boost Python is comprehensive and well accepted by the Python community.
I would go for it for its extensive use and testing. I would decline it
for being hard to debug and having a complicated building system. I
don't think it would be worth including a boost dependency just for the
sake of creating the Python bindings, but since Mocapy++ already depends
on Boost, using it becomes a more attractive option. In my personal
experience, Boost Python is very mature and there are no limitations on
what one can do with it. When it comes to performance, Cython still
overcomes it. Have a look at the [Cython C++ wrapping
benchmarks](http://blog.behnel.de/index.php?p=38) and check the timings
of [Cython against Boost Python](http://www.behnel.de/cycppbench/).
There are also previous [benchmarks comparing Swig and Boost
Python](http://telecom.inescporto.pt/~gjc/pybindgen-benchmarks/).

#### Cython

It is incredibly faster than other options to create python bindings to
C++ code, according to several benchmarks available on the web. Check
the Simple benchmark between Cython and Boost.Python. It is also very
clean and simple, yet powerful. Python's doc on porting extension
modules mentions cython: "If you are writing a new extension module, you
might consider Cython." Cython has now support for efficient interaction
with numpy arrays. it is a young, but developing language and I would
definitely give it a try for its leanness and speed.

Since Boost is well supported and Mocapy++ already relies on it, we
decided to use Boost.Python for the bindings.

For further information see [Mocapy++Biopython - Box of
ideas](https://docs.google.com/document/d/1E72Qysp3pMd69hSYfIXJKgBLdeSMbzSol9RYD2rKHlI/edit?hl=pt_BR&authkey=CPmFxK0H).

### Bindings Prototype

The source code for the prototype is on the gSoC11 branch:
<http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/bindings_prototype/>

Bindings for a few Mocapy++ features and a couple of examples to find
possible implementation and performance issues.

**Procedure**

-   Implemented the examples hmm\_discrete and
    discrete\_hmm\_with\_prior in Python, assuming the interface
    Mocapy++ already provides.

<!-- -->

-   Implemented the bindings to provide a minimum subset of
    functionality, in order to run the implemented examples.

<!-- -->

-   Compared the performance of C++ and Python versions.

Mocapy++’s interface remained unchanged, so the tests look similar to
the ones in Mocapy/examples.

In the prototype the bindings were all implemented in a single module.
For the actual implementation, we could mirror the src packages
structure, having separated bindings for each package such as discrete,
inference, etc.

It was possible to implement all the functionality required to run the
examples. It was not possible to use the
[vector\_indexing\_suite](http://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/indexing.html)
when creating bindings for vectors of MDArrays. A few operators (in the
MDArray) must be implemented in order to export indexable C++ containers
to Python.

Two Mocapy++ examples that use discrete nodes were implemented in
Python. There was no problem in exposing Mocapy’s data structures and
algorithms. The performance of the Python version is very close to the
original Mocapy++.

For additional details have a look at the [Mocapy++ Bindings
Prototype](https://docs.google.com/document/d/1JPkCbvJ9Gk3b6LmQ68__UUt4Yz2P1-c286p6FVEqyXs/edit?hl=pt_BR&authkey=CJTCpZgL)
report.

### Bindings Implementation

#### Bindings for the core functions and data structures

''' Data structures '''

Mocapy uses an internal data structure to represent arrays: MDArray. In
order to make it easier for the user to interact with Mocapy's API, it
was decided to provide an interface that accepts numpy arrays.
Therefore, it was necessary to implement a translation between a numpy
array and an MDArray.

The translation from MDArray to python was done through the use of
Boost.Python
[to\_python\_converter](http://www.boost.org/doc/libs/1_43_0/libs/python/doc/v2/to_python_converter.html).
We've implemented a template method convert\_MDArray\_to\_numpy\_array,
which converts an MDArray of any basic type to a corresponding numpy
array. In order to perform the translation the original array's shape
and internal data are copied into a new numpy array.

The numpy array was created using the [Numpy Array
API](http://docs.scipy.org/doc/numpy/reference/c-api.array.html). The
creation of a new
[PyArrayObject](http://docs.scipy.org/doc/numpy/reference/c-api.types-and-structures.html#PyArrayObject)
using existing data (PyArray\_SimpleNewFromData) doesn't copy the array
data, it just stores a pointer to it. Thus, one can only free the data
when there is no reference to the object. This was done through the use
of a [Capsule](http://docs.python.org/c-api/capsule.html). Besides
encapsulating the data, the capsule also stores a destructor to be used
when the array is destroyed. The PyArrayObject has a field named "base"
which points to the capsule.

The translation from Python to C++, i.e. creating an MDArray from a
numpy array is slightly more complex. Boost.Python will provide a chunk
of memory into which the new C++ object must constructed in-place. See
the [How to write boost.python
converters](http://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/)
article for more details.

A translation between std::vector of basic types (double, int...) and
Python list was also implemented. For std::vector of custom types, such
as Node, the translation to a Python list was not performed. If done the
same ways as for basic types, a type error is raised: "TypeError: No
to\_python (by-value) converter found for C++ type". When using
[vector\_indexing\_suite](http://www.boost.org/doc/libs/1_46_1/libs/python/doc/v2/indexing.html#vector_indexing_suite_class)
this problem was already solved. See [Wrapping
std::vector<AbstractClass*>](http://mail.python.org/pipermail/cplusplus-sig/2005-July/008865.html).
The only inconvenience of using the vector\_indexing\_suite is creating
new types such as vector\_Node, instead of using a standard Python list.

The code for the translations is in the
[mocapy\_data\_structures](http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/bindings/mocapy_data_structures.cpp?revision=340&view=markup)
module.

''' Core functions '''

The mocapy Python packages follow Mocapy's current source tree. For each
package, a shared library with the bindings was created. This makes
compilation faster and debug easier. Also, if a single library was
created it wouldn't be possible to define packages.

Each of the libraries is called libmocapy\_<nameofthepackage>. For
example, libmocapy\_gaussian provides bindings for the gaussian nodes
and probability distributions. The libmocapy\_data\_structures is used
by other libraries and, therefore, must be imported first. This is done
on the Python side. Each of the libmocapy\_\* libraries is imported in
the corresponding package. See [Creating
Packages](http://www.boost.org/doc/libs/1_46_1/libs/python/doc/tutorial/doc/html/python/techniques.html#python.creating_packages).

The bindings code can be found in the [Bindings
directory](http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/bindings/).

Currently, tests to the just created interface are being developed.
There are a few tests already implemented under the framework package:
[mocapy/framework/tests](http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/python/mocapy/framework/tests/)

#### Bindings for the remaining Mocapy++ functionality

''' Data structures '''

While implementing the bindings for the remaining Mocapy++ functionality
there were problems with methods that take pointers and references to an
mdarray:

-   It is not possible to call a method which takes a pointer if the
    object is created on the python side. See [how to call a function
    that expects a
    pointer?](http://stackoverflow.com/questions/3881457/boostpython-howto-call-a-function-that-expects-a-pointer).

<!-- -->

-   It is not possible to automatically translate a non const reference.
    The custom rvalue converters only match functions with the following
    signatures:

<cpp> void foo(std::vector<double> const& array); // pass by
const-reference

void foo(std::vector<double> array); // pass by value </cpp>

For further details see [How can I wrap functions which take C++
containers as
arguments?](http://www.boost.org/doc/libs/1_46_1/libs/python/doc/v2/faq.html#question2)

The mdarray is created in python using a numpy.array that is translated
to c++ using [custom
converters](http://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/faq.html#custom_string).
The custom converters are registered in the global Boost.Python registry
near the top of the module initialization function. Once flow control
has passed through the registration code the automatic conversions from
and to Python.

Because of this automatic conversions, it was necessary to create
wrappers for functions which take pointers as arguments and change the
functions which take references, to get const references. Because
Mocapy++ is not [const
correct](http://en.wikipedia.org/wiki/Const-correctness), changes are
needed to use the const references properly. While the changes are being
done, some const\_cast have been used. When using const\_cast one must
be aware [it is not always
safe](http://stackoverflow.com/questions/357600/is-const-cast-safe).

The [call
policies](http://www.boost.org/doc/libs/1_46_1/libs/python/doc/tutorial/doc/html/python/functions.html#python.call_policies)
were also reviewed. When using an incorrect return value policy, you
won't get a compile error, but your code will crash at runtime.

''' Examples '''

Mocapy++'s examples were implemented in Python, using the exposed API
and data type conversions.
<http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/python/examples/>

#### Testing and Improving Mocapy Bindings

Before integrating to Biopython, some unit testing was required, to
detect possible errors and make sure future changes that break
functionality won't go unnoticed.

For every Python package, it was created a "tests" directory which
contains the unit tests created for each module. Here is one example of
the tests created for the framework package:
<http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/python/mocapy/framework/tests/>

While testing the code, a few issues were detected:

-   **Object ownership:**

When passing an object created on the C++ side to a method that takes a
pointer as an argument, one should be careful about the life time of
that object.

For example, the set\_random\_gen method takes a pointer to a RandomGen
object. The following code works just fine.

``` python
random_gen = RandomGen()
node.set_random_gen(random_gen=random_gen)
```

But if instead of doing that, we do the following:

``` python
node.set_random_gen(random_gen=RandomGen())
```

The reference count has not been incremented and therefore the object
can be destroyed.

The way to solve the problem is to make sure the C++ object is held by
auto\_ptr: <cpp> class\_&lt;RandomGen, std::auto\_ptr<RandomGen>
&gt;("RandomGen") </cpp>

Then make a thin wrapper function which takes an auto\_ptr parameter:
<cpp> void node\_set\_random\_gen(Node& node, std::auto\_ptr<RandomGen>
random\_gen) {

`   node.set_random_gen(random_gen.get());`
`   node.release();`

} </cpp> For further details, see [How can I wrap a function which needs
to take ownership of a raw
pointer?](http://www.boost.org/doc/libs/1_46_0/libs/python/doc/v2/faq.html#ownership)

Pointers returned via
[manage\_new\_object](http://wiki.python.org/moin/boost.python/CallPolicy#manage_new_object)
will also be held by auto\_ptr, so the transfer-of-ownership works
correctly. When using this call policy the caller is responsible for
deleting the C++ object from the heap.

-   *' Translation from numpy.array to a float mdarray*'

If the numpy array is an integer array, the translation creates an
mdarray<int> and this is passed to a method which expects an mdarray of
floats. This generates incorrect results.

The way to deal with that from the user perspective is either using
floating pointer numbers to create the array or setting the ndtype
parameter when creating the array:

``` python
x = numpy.array([[1, 2, 3, 4, 5, 6]], dtype=numpy.float64)
```

#### Building and Distributing Mocapy as a Package

[Distutils](http://docs.python.org/distutils/index.html) was used to
distribute Mocapy's Python modules.

Besides distributing the python code, it was also necessary to build the
extension modules. [Building Extensions with
boost.python](http://wiki.python.org/moin/boost.python/BuildingExtensions)
describes ways to build extensions using distutils.

Mocapy's setup.py can be found at
<http://mocapy.svn.sourceforge.net/viewvc/mocapy/branches/gSoC11/python/setup.py?revision=418&view=markup>

Using the setup script, mocapy installation is done in a few steps:

-   Build the mocapy library using cmake (usual procedure described in
    mocapy docs);
-   Issue "python setup.py build", to build the extension modules;
-   Issue "python setup.py install", to install the package (normally,
    the install procedure does the step above in case you didn't).

### Integration with Biopython

#### API Design

In order to use Mocapy in conjunction with Biopython, a new module for
PDB-specific features was added to Bio.PDB. This is where the API is
being designed.

Mocapy is added as an optional dependency in Biopython. Inside the
function or module that requires Mocapy, "import mocapy" is wrapped in a
try/except block. A MissingPythonDependencyError is issued if the import
fails.

Things that are being studied to be included in the module:

-   extract the backbone dihedral angles from a given set of structures;
-   use this data to train a TorusDBN-like model;
-   automatically decide on the best model using the BIC criterion.

#### Barnacle

Frellsen J, Moltke I, Thiim M, Mardia KV, Ferkinghoff-Borg J, et al.
2009 [A Probabilistic Model of RNA Conformational
Space](http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000406).
PLoS Comput Biol 5(6): e1000406. <doi:10.1371/journal.pcbi.1000406>.

RNA 3-D structure prediction methods require an accurate energy function
and a conformational sampling procedure. Barnacle focuses on the problem
of conformational sampling.

The aim of BARNACLE (BAyesian network model of RNA using Circular
distributions and maximum Likelihood Estimation) is to capture both the
marginal distributions of each of the angles and the local dependencies
between them. Barnacle describes RNA structure in a natural continuous
space. It can be used purely as a proposal distribution, but also as an
energy term enforcing realistic local conformations. The model combines
a dynamic Bayesian network (DBN) with directional statistics.

<img src="Journal.pcbi.1000406.g002.png" title="Barnacle DBN (doi:10.1371/journal.pcbi.1000406.g002)" alt="Barnacle DBN (doi:10.1371/journal.pcbi.1000406.g002)" width="600" />

The DBN represents nine consecutive dihedral angles, where the seven
central angles originate from a single nucleotide. Each slice j (a
column of three variables) corresponds to one dihedral angle in an RNA
fragment. The variables in each slice are: an angle identifier, Dj, a
hidden variable, Hj, and an angular variable, Aj. The angle identifier
keeps track of which dihedral angle is represented by a slice, while the
angular node models the actual dihedral angle value. The hidden nodes
induce dependencies between all angles along the sequence (and not just
between angles in consecutive slices).

The original source code for Barnacle, which contains an embedded
version of Mocapy written in Python, can be found at
<http://sourceforge.net/projects/barnacle-rna>.

The modified version of Barnacle, changed to work with the Mocapy
bindings can be found at
<https://github.com/mchelem/biopython/tree/master/Bio/PDB/Barnacle>.

Here is an example of use:

``` python
model = Barnacle("ACCU")
model.sample()
print("log likelihood = ", model.get_log_likelihood())
model.save_structure("structure01.pdb")
```

#### TorusDBN

Wouter Boomsma, Kanti V. Mardia, Charles C. Taylor, Jesper
Ferkinghoff-Borg, Anders Krogh, and Thomas Hamelryck. A generative,
probabilistic model of local protein structure. Proc Natl Acad Sci U S
A. 2008 July 1; 105(26): 8932–8937.
<http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2440424/>

TorusDBN aims at predicting the 3D structure of a biomolecule given its
amino-acid sequence. It is a continuous probabilistic model of the local
sequence–structure preferences of proteins in atomic detail. The
backbone of a protein can be represented by a sequence of dihedral angle
pairs, φ and ψ that are well known from the [Ramachandran
plot](http://en.wikipedia.org/wiki/Ramachandran_plot). Two angles, both
with values ranging from −180° to 180°, define a point on the torus.
Hence, the backbone structure of a protein can be fully parameterized as
a sequence of such points.

<img src="Torus_dbn.png" title="TorusDBN (doi: 10.1073/pnas.0801715105)" alt="TorusDBN (doi: 10.1073/pnas.0801715105)" width="600" />

The circular nodes represent stochastic variables. The rectangular boxes
along the arrows illustrate the nature of the conditional probability
distribution between them. A hidden node emits angle pairs, amino acid
information, secondary structure labels and cis/trans information.

The TorusDBN model is originally implemented as part of the backboneDBN
package, which is freely available at
<http://sourceforge.net/projects/phaistos/>.

A new version of the TorusDBN model was implemented in the context of
this project and can be found at
<https://github.com/mchelem/biopython/tree/master/Bio/PDB/TorusDBN>.

The TorusDBNTrainer can be used to train a model with a given training
set:

``` python
trainer = TorusDBNTrainer()
trainer.train(training_set)  # training_set is a list of files
model = trainer.get_model()
```

Then the model can be used to sample new sequences:

``` python
model.set_aa("ACDEFGHIK")
model.sample()
print(model.get_angles())  # The sampled angles.
```

When creating a model, it is possible to create a new DBN specifying the
size of the hidden node or loading the DBN from a file.

``` python
model = TorusDBNModel()
model.create_dbn(hidden_node_size=10)
model.save_dbn("test.dbn")
```

``` python
model = TorusDBNModel()
model.load_dbn("test.dbn")
model.set_aa("ACDEFGHIK")
model.sample()
print(model.get_angles())  # The sampled angles.
```

It is also possible to choose the best size for the hidden node using
the find\_optimal\_model method:

``` python
trainer = TorusDBNTrainer()
hidden_node_size, IC = trainer.find_optimal_model(training_set)
model = trainer.get_model()
```

IC is either the [Bayesian Information
Criterion](http://en.wikipedia.org/wiki/Bayesian_information_criterion)
(BIC) or the [Akaike Information
Criterion](http://en.wikipedia.org/wiki/Akaike_information_criterion)
(AIC) (Defaults to BIC. AIC can be specified by setting the use\_aic
flag).

For more details on the model API, see the test files:
<https://github.com/mchelem/biopython/blob/master/Tests/test_TorusDBNTrainer.py>
and
<https://github.com/mchelem/biopython/blob/master/Tests/test_TorusDBNModel.py>.

### Performance

A few performance measurements were made comparing test cases
implemented both in C++ and in Python. The tests were run in a computer
with the following specification: Core 2 Duo T7250 2.00GHz, Memory Dual
Channel 4.0GB (2x2048) 667 MHz DDR2 SDRAM, Hard Drive 200GB 7200RPM.

There were no significant performance differences. For both
implementations the methods responsible for consuming most cpu time were
the same:

<img src="Hmm_discrete.png" title="fig:DBN with discrete nodes, C++ implementation " alt="DBN with discrete nodes, C++ implementation " width="400" />
<img src="Hmm_discrete_py.png" title="fig:DBN with discrete nodes, Python implementation " alt="DBN with discrete nodes, Python implementation " width="400" />

The profiling tests were made using
[Callgrind](http://valgrind.org/info/tools.html#callgrind) and
visualized using [Kcachegrind](http://kcachegrind.sourceforge.net/).

Here are the average running time of the examples available with Mocapy
(10 runs):

| Test name                  | C++ (s) | Python (s) |
|----------------------------|---------|------------|
| hmm\_simple                | 0.52    | 0.58       |
| hmm\_discrete              | 48.12   | 43.45      |
| discrete\_hmm\_with\_prior | 55.95   | 50.09      |
| hmm\_dirichlet             | 340.72  | 353.98     |
| hmm\_factorial             | 0.01    | 0.12       |
| hmm\_gauss\_1d             | 53.97   | 63.39      |
| hmm\_gauss                 | 16.02   | 16.96      |
| hmm\_multinomial           | 134.64  | 125.83     |
| hmm\_poisson               | 11.00   | 10.60      |
| hmm\_vonmises              | 7.22    | 7.36       |
| hmm\_torus                 | 53.79   | 53.65      |
| hmm\_kent                  | 61.35   | 61.06      |
| hmm\_bippo                 | 40.66   | 41.81      |
| infenginehmm               | 0.01    | 0.12       |
| infenginemm                | 0.01    | 0.15       |

#### TorusDBN

Even though the PDB files are read, parsed and transformed in a format
mocapy can understand, the most time consuming methods are the ones
performing mathematical operations during the sampling process
(Chebyshev and exp, for example).

<img src="TorusDBN.png" title="Training of the TorusDBN model " alt="Training of the TorusDBN model " width="400" />

The model has been trained with a training set consisting of about 950
chains with maximum 20% homology, resolution below 1.6 Å and R-factor
below 25%. It took about 67 minutes to read and train the whole dataset.

The resulting DBN is available at
<https://github.com/mchelem/biopython/blob/master/Tests/TorusDBN/pisces_dataset.dbn>
and can be loaded directly into the model as explained in the TorusDBN
section above.

### Future work

The summer is over, but the work continues... There are still a lot of
things I intend to work on:

-   Test the trained models to check their effectiveness in protein
    structure prediction.

<!-- -->

-   Try to reduce dynamic allocation as it is responsible for a lot of
    running time.

<!-- -->

-   Guarantee there are no memory leaks in the bindings.

