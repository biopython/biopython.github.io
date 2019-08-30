---
title: Contributing to Biopython.
permalink: wiki/Contributing
layout: wiki
---

### A Guide to Contributing to Biopython

So you want to contribute to Biopython, huh? Great! New contributions
are the lifeblood of the project. However, if done incorrectly, they can
quickly suck up valuable developer time. (We have day jobs too!) This is
a short guide to the recommended way to contribute code to Biopython.

See also the chapter about contributing in the
[Tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
([PDF](http://biopython.org/DIST/docs/tutorial/Tutorial.pdf)).

If you get as far as making a pull request with changes, you should
read [CONTRIBUTING.rst](https://github.com/biopython/biopython/blob/master/CONTRIBUTING.rst)
which describes how the automated testing including style checks are
done.

### Non-code contributions

Even if you don't feel ready or able to contribute code, you can still
help out. There always things that can be improved on the
[documentation](Documentation "wikilink") (even just proof reading, or
telling us if a section isn't clear enough). We also need people on the
[mailing lists](Mailing_lists "wikilink") to help out with beginner's
questions, or to participate in debates about new features. Maybe you
can propose general examples for the [wiki
cookbook](Category%3ACookbook "wikilink")?

### Finding a Project

The best contributions are the code that you have already been using in
your daily research. This should be code that you think might be useful
for other people and is already free of bugs. If you are thinking of
sending this in, go on to Step 2, Submitting Code!

Otherwise, there are still many ways to contribute to Biopython, both
involving coding and not. Some things that you can help on include:

-   **More unit tests:** Some of our modules still only have partial
    unit test coverage.
-   **Support for More Programs:** There are many different
    bioinformatics programs being developed. Identify one that does not
    currently have support in Biopython and add support for it.
-   **Support for More File Formats:** We can read and write lots of
    different file formats, but there are always more. For sequences and
    alignments look at the [SeqIO](SeqIO "wikilink") and
    [AlignIO](AlignIO "wikilink") pages first. Note that HTML parsers
    for specific websites are discouraged as these require long
    term maintenance.
-   **Support for Databases:** Identify a biological database that does
    not currently have support in Biopython and add support for it.
-   **Add New Data Type:** You can add code that works with a new type
    of data. This is a tough area, though. Creating a new robust and
    useful data type is difficult, and we may be hesitant to add
    something unless it's already tried and tested.
-   **Add New Algorithm:** You can add a new well-known algorithm that
    might be useful for other biologists.
-   **Parser Verification:** As Biopython supports more and more
    databases, the difficulty in maintaining the format
    parsers increases. These formats are changing very quickly. Thus, we
    need to periodically verify that the parsers are still working. For
    example, the GenBank parser needs to be checked to make sure it
    handles each new dump of GenBank.
-   **Regression Tests:** Biopython uses a regression testing framework
    to make sure code is working. Although most of the functionality in
    Biopython is tested, there are still some holes.
-   **Documentation:** The tutorial is not complete and can use
    some work. New users can be especially helpful here, as you learn
    new packages. Our API documentation could also do with some work,
    see coding conventions below.
-   **News Postings:** If you keep up with the mailing lists (which are
    fairly low volume), we need someone to help summarize important
    posts and events as [news items](News "wikilink").

### Submitting Code

In general, we will consider any code that is applicable to biological
or chemical data. Please do not submit code whose functionality largely
overlaps with code already in Biopython, unless there is an obvious
improvement and you have a plan for integrating the code.

Before you submit it, please check that:

-   It is generalized and likely to be useful for other things.
-   The dependencies are reasonable. Adding support for a third party
    command line tool is fine, but requiring additional python libraries
    needs discussion. Dependencies on commercial closed-source software
    probably won't be accepted to Biopython.
-   The code will be licensed with the Biopython license.
-   You must have the legal right to contribute it and license it under
    the Biopython license.
-   You are enthusiastic about maintaining it and responding to
    bug reports.
-   You included docstrings in the code, and are willing to
    write documentation.
-   You have written, or are willing to write, a unit test for the
    new code.

If all these terms seem acceptable, please send a description of your
code to [biopython@biopython.org](mailto:biopython@biopython.org) (see [Mailing
lists](Mailing_lists "wikilink")). Be sure to subscribe to biopython mailing list
before sending a message, otherwise your message will be discarded by
the mail server (this was done to avoid spam on the mailing list). Don't
send the code directly to the biopython mailing list. Instead,
please use [our GitHub page](https://github.com/biopython/biopython/) by
creating an issue and either attaching the file(s) or linking to your
GitHub branch.

### Coding conventions

Biopython tries to follow the coding conventions laid out in
[PEP 8](http://www.python.org/dev/peps/pep-0008/) and
[PEP 257](http://www.python.org/dev/peps/pep-0257/). The important
highlights are:

-   Classes should be in AllFirstLetterUppercase style.
-   Functions should be in lowercase\_separated\_by\_underscores style.
-   Variables are either in lowercase\_separated\_by\_underscores or
    lowercasemungedtogether style, depending on your preferences and the
    length of the variable.
-   \_single\_leading\_underscores to indicate internal functions or
    classes that shouldn't be called directly be a user.
-   Tabs are bad. Most people in the Python community now dislike tabs
    and instead prefer using 4 spaces for indentation. Most editors can
    help you take care of this (Emacs python-mode uses the 4 space rule,
    for instance). Tools/scripts/reindent.py in the Python distribution
    will help get rid of tabs in files.

The one notable exception is module names, where we tend to use title
case. With hindsight this is unfortunate, but we are constrained by
backwards compatibility.

[Epydoc](http://epydoc.sourceforge.net/) and/or
[Sphinx apidoc](https://www.sphinx-doc.org/en/master/man/sphinx-apidoc.html)
is being used to generate
automatic documentation of the source code so it definitely is useful to
put helpful comments in your code so that they will be reflected in [the
API documentation](http://biopython.org/DIST/docs/api) (in addition to
all the normal reasons to document code).

We generally don't do anything fancy to try and format the comments in
the code -- but they are interpreted as reStructuredText markup which
allows things like bullet points and italics. This isn't fancy, but
it's effective and easier then trying to deal with the myriad of
different ways to try and structure source code comments.

However, there are a few tricks to make your documentation look its
best. The main ones are:

-   Modules, classes and function documentation should start with a one
    line description. This must end with a period.

Here's an example of a module documented so that epydoc will be happy:

``` python
"""This is a one line description of the module followed with a period.

More information about the module and its goals and usage.
"""

class MyClass:
    """One line description of the class followed by a period.

    More information about the class -- its purpose, usage, and
    implementation.
    """

    def my_function(self, spam):
        """A terse description of my function followed with a period.

        A longer description with all kinds of additional goodies. This may
        include information about what the function does, along with
        what parameters it will be passed and what it returns. You know,
        information so people know how to use the function.
        """
        #the code ...
```
