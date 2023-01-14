---
title: Deprecation policy
permalink: wiki/Deprecation_policy
layout: wiki
---

As bioinformatics and computational biology are developing quickly, some
of previously developed Biopython modules may no longer be relevant in
today's world. To keep the code base clean, we aim to deprecate and
remove such code from Biopython, while avoiding any nasty surprises for
users who may be relying on older code.

We use the same process to deprecate obsolete modules, methods, functions
and classes.

We keep a plain text file in the Biopython source code to record these
changes, named
[`DEPRECATED.rst`](http://github.com/biopython/biopython/blob/master/DEPRECATED.rst).

This is the current policy for deprecating and removing code from
Biopython:

-   First, ask on the [Biopython mailing list](Mailing_lists "wikilink")
    whether a given piece of code has any users. Please keep in mind that
    not all users are following the biopython mailing list.
-   Consider declaring the module as "obsolete" for a release
    *before* deprecation. No code changes, just:
    -   Note this in the DEPRECATED file,
    -   Add "(OBSOLETE)" to the first line of the module docstring,
    -   Use the module docstring to explain why the module is obsolete
        and what should be used instead.
    -   As use of obsolete code should be avoided in the Biopython code
        base, start replacing obsolete code by the newly recommended
        code. This also serves to verify whether the code can actually
        be deprecated and removed. To track changes more easily, multiple
        pull requests may be used to replace obsolete code.
    -   Start reorganizing or rewriting the Biopython documentation to
        point users to the newly recommended code, and away from the
        deprecated code. However, during the deprecation process a
        description of the obsolete code should remain in the
        documentation until the code is actually removed.
-   If there are no apparent users, then actually deprecate it:
    -   Note this in the DEPRECATED file,
    -   Add "(DEPRECATED)" to the first line of the module docstring
    -   Use the module docstring to explain any migration needed,
        ideally with examples or a reference to the tutorial.
    -   Most importantly, add a ``BiopythonDeprecationWarning`` to the
        code (Python's ``DeprecationWarning`` is silent by default):

``` python
import warnings
from Bio import BiopythonDeprecationWarning

warnings.warn(
    "Bio.SomeModule has been deprecated, and we intend to remove it"
    " in a future release of Biopython. Please use the SomeOtherModule"
    " instead, as described in the Tutorial. If you would like to"
    " continue using Bio.SomeModule, please contact the Biopython"
    " developers via the mailing list or GitHub.",
    BiopythonDeprecationWarning,
)
```

-   Silence the new warning in the tests, for example if your tests use
    ``from Bio import SomeModule`` replace that with:

``` python
from Bio import BiopythonDeprecationWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import SomeModule
```

-   In principle, we require that two Biopython releases carrying the
    deprecation warning are made before the code can be
    actually removed.
-   In addition, at least one year should pass between the first
    Biopython release carrying the deprecation warning, and the first
    Biopython release in which the code has been actually removed.

See here for the discussion on the mailing list:

<http://lists.open-bio.org/pipermail/biopython-dev/2008-November/004920.html>
