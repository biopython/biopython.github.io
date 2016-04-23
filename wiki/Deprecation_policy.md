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

We keep a plain text file in the Biopython source code to record these
changes, available on 
[GitHub](http://github.com/biopython/biopython/blob/master/DEPRECATED) or
[here](http://biopython.open-bio.org/SRC/biopython/DEPRECATED).

This is the current policy for deprecating and removing code from
Biopython:

-   First, ask on the biopython and biopython-dev [mailing
    lists](Mailing_lists "wikilink") whether a given piece of code has
    any users. Please keep in mind that not all users are following the
    biopython-dev mailing list.
-   Consider declaring the module as "obsolete" for a release
    *before* deprecation. No code changes, just:
    -   Note this in the DEPRECATED file,
    -   Add "(OBSOLETE)" to the first line of the module docstring,
    -   Use the module docstring to explain why the module is obsolete
        and what should be used instead.
-   If there are no apparent users, then actually deprecate it:
    -   Note this in the DEPRECATED file,
    -   Add "(DEPRECATED)" to the first line of the module docstring
    -   Use the module docstring to explain any migration needed,
        ideally with examples or a reference to the tutorial.
    -   Most importantly, add a DeprecationWarning to the code:

```
import warnings
warnings.warn("Bio.SomeModule has been deprecated, and we intend to remove it"
              " in a future release of Biopython. Please use the SomeOtherModule" 
              " instead, as described in the Tutorial. If you would like to"
              " continue using Bio.SomeModule, please contact the Biopython"
              " developers via the mailing list.",
              DeprecationWarning)
```

-   In principle, we require that two Biopython releases carrying the
    deprecation warning are made before the code can be
    actually removed.
-   In addition, at least one year should pass between the first
    Biopython release carrying the deprecation warning, and the first
    Biopython release in which the code has been actually removed.

See here for the discussion on the mailing list:

<http://lists.open-bio.org/pipermail/biopython-dev/2008-November/004920.html>
