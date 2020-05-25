---
title: How to build a release for Biopython.
permalink: wiki/Building_a_release
layout: wiki
---

Build Biopython with many small steps!!

Setup required for a new release manager
----------------------------------------

The instructions below require that you have access to a few servers and
the code repository. When you start, be sure to have write access to:

1. [GitHub Biopython source code repository](https://github.com/biopython/biopython)
2. [GitHub biopython.org website repository](https://github.com/biopython/biopython.github.io)
3. [GitHub biopython.org/DIST/ repository](https://github.com/biopython/DIST)
4. [GitHub biopython wheel repository](https://github.com/biopython/biopython-wheels)
5. [OBF WordPress Blog](https://www.open-bio.org)
6. [Biopython on PyPI](https://pypi.python.org/pypi/biopython)

If you don't have any of the above, please ask.

We assume you have cloned these under ``~/repositories/`` and that the
``git origin`` is the official Biopython copy of the repository.

We assume you are running Linux (but macOS should be fine too),
and have the following tools installed (plus as many of the
Biopython optional dependencies as possible for local testing):

1. Python 3
2. git
3. [twine](https://github.com/pypa/twine/), installed with ``pip install twine``
4. LaTeX, including assorted packages like comments and preprint.
5. [hevea](http://hevea.inria.fr/), I am currently using version 1.10+9 of 2008-12-17

Final commit(s)
---------------

1. Using git make sure I have the latest code:

   ``` bash
   $ cd ~/repositories/biopython
   $ git checkout master
   $ git pull origin master
   ```

2. make sure the `README.rst` file is still up to date

3. add any important info to `NEWS.rst` or `DEPRECATED.rst` - you can get a
   log of recent git changes like this (adjust the date accordingly):

   ``` bash
   $ git log --since="2016/01/01" --reverse --pretty="medium"
   ```

4. make sure `CONTRIB.rst` still current

5. make sure `setup.py` and `MANIFEST.in` are still up to date

   - Are there any new modules/files which should get installed?

6. bump version numbers and set the release data:

   - Biopython version - edit `Bio/__init__.py`
   - Biopython Tutorial - update the date/version line in the
     `Doc/Tutorial.tex` file
   - Biopython `NEWS.rst` - fill in the release date
   - Make sure to commit the modified files to github, e.g.

   ``` bash
   $ git commit Bio/__init__.py Doc/Tutorial.tex NEWS.rst -m "Call this Biopython 1.68"
   ```

7. do a final check to make sure things are checked in:

   ``` bash
   $ rm -r build
   $ rm Tests/*.pyc
   $ make clean -C Doc
   $ git status
   ```

8. build Biopython and run a a final regression test:

   ``` bash
   drevil:~biopython> python setup.py build
   drevil:~biopython> python setup.py test
   ```
   Running the tests simultaneously is risky as two threads may both try to read/write to the same temp files.

9. Push this to gitub, all being well this commit will be tagged as the release
   (barring no problems uncovered while building the documentation, or
   with the manifest while testing the tar-ball):

   ``` bash
   $ git push origin master
   ```

Making and testing the tar-ball
-------------------------------

10. check out a clean version somewhere else:

    ``` bash
    drevil:~tmp1/> git clone https://github.com/biopython/biopython.git
    drevil:~tmp1/> cd biopython
    ```

11. make the documentation PDF, text and HTML files in Doc:

    ``` bash
    drevil:~tmp1/biopython/> make -C Doc
    drevil:~tmp1/biopython/> make clean -C Doc
    ```

12. make the source distribution

    ``` bash
    drevil:~tmp1/biopython> python setup.py sdist --formats=gztar,zip
    ```

13. untar the file somewhere else

    ``` bash
    drevil:~tmp1/biopython/> cd ..
    drevil:~tmp1/> tar -xzvf biopython/dist/biopython-1.71.tar.gz
    drevil:~tmp1/> cd biopython-1.71
    ```

    Check to make sure it includes the HTML and PDF files under Doc

14. make sure I can build and test it

    ``` bash
    drevil:~tmp1/biopython-1.71/> python setup.py build
    drevil:~tmp1/biopython-1.71/> python setup.py test
    drevil:~tmp1/biopython-1.71/> python setup.py install --prefix /tmp/test-install
    ```

    A typical source of failure here (on the tests) is the lack of example
    files being added to the source distribution: add them to `MANIFEST.in`

Checking the compiled documentation
-----------------------------------

15. Since Biopython 1.74, Sphinx has handled the API documentation via continuous
    integration, but you still have to update the Tutorial on the website manually.

  - Update the tutorial:

    ``` bash
    $ cd ~/repositories/DIST/docs/tutorial/
    $ cp .../tmp1/biopython/Doc/Tutorial.html .
    $ cp .../tmp1/biopython/Doc/Tutorial.pdf .
    $ git commit Tutorial.html Tutorial.pdf -m "Tutorial for Biopython 1.71"
    ```

    - Push this to GitHub Pages to update the website:

    ``` bash
    $ git push origin gh-pages
    ```

    - Check this is live at <http://biopython.org/DIST/docs/api/Bio-module.html>,
      <http://biopython.org/DIST/docs/tutorial/Tutorial.html>, and
      <http://biopython.org/DIST/docs/tutorial/Tutorial.pdf>

Making wheels
-------------


16. Now we use https://github.com/biopython/biopython-wheels to build wheels,
    by updating the ``BUILD_COMMIT`` line in ``.travis.yml`` and ``appveyor.yaml``
    to the new release's commit hash (which all being well will get a git tag).

    ``` bash
    $ cd ~/repositories
    $ git clone git@github.com:biopython/biopython-wheels.git
    $ cd biopython-wheels/
    $ git submodule update --init
    $ emacs .travis.yml  # update BUILD_COMMIT=... line
    $ emacs appveyor.yml  # update BUILD_COMMIT: ... line
    $ git commit .travis.yml appveyor.yml -m "Build Biopython 1.xx"
    $ git push origin master
    ```

    Check the wheels build
    [on TravisCI for Linux and Mac](https://travis-ci.org/biopython/biopython-wheels/builds) and
    [on AppVeyor for Windows](https://ci.appveyor.com/project/biopython/biopython-wheels/history).

    You don't seem to need to update the ``biopython`` git submodule, but if you
    need to update ``multibuild`` this seems to work.

    ``` bash
    $ git submodule foreach git pull origin master
    $ git commit -a -m "Update submodules"
    $ git push origin master
    ```

17. Successful wheels will be on
    [Rackspace](https://a365fff413fe338398b6-1c8a9b3114517dc5fe17b7c3f8c63a43.ssl.cf2.rackcdn.com/),
    download them from there to your ``~/repository/biopython/DIST/`` folder.
    We will upload these to PyPI later using Twine.

18. If you have a Windows machine, remove any prior Biopython installations,
    and confirm the Windows wheel file(s) work.

Tagging the release, and uploading
----------------------------------

19. Back in the main repository, tag the release:

    ``` bash
    $ cd  .../tmp1/biopython/
    $ git tag biopython-171
    $ git push origin master --tags
    ```

20. Upload the new release tar-ball and zip to the website via GitHub Pages `DIST` repository.

    ``` bash
    $ cp dist/biopython-1.68.* ~/repositories/DIST/
    $ cd ~/repositories/DIST/
    $ git add biopython-1.71.*
    $ git commit biopython-1.71.* -m "Downloads for Biopython 1.71"
    $ shasum -a 256 biopython-1.71.*
    $ md5sum biopython-1.71.*
    $ git commit --amend # paste checksums into comment
    $ git push origin gh-pages
    ```

21. Upload to the python package index (except for beta/alpha level
releases):

    ``` bash
    $ cd  ~/repositories/biopython/
    $ pip install twine
    $ twine upload dist/biopython-1.71.tar.gz
    $ twine upload dist/biopython-1.71-*.whl
    ```

    - You need to have a login on pypi and be registered with Biopython to be
      able to upload the new version

    - Check this is live at <https://pypi.python.org/pypi/biopython/>

22. Update the website:

    - If you haven't already, clone the ``biopython.github.io`` repository,
      (otherwise make sure your copy is up to date):

    ``` bash
    $ cd ~/repositories
    $ git clone git@github.com:biopython/biopython.github.io.git
    ```

    - Update the website:

    ``` bash
    $ cd ~/repositories/biopython.github.io
    $ emacs _config.yml
    $ emacs wiki/Biopython.md
    $ emacs wiki/Download.md
    $ emacs wiki/Documentation.md # update API docs link
    $ git commit _config.yml wiki/Biopython.md wiki/Download.md -m "Biopython 1.68 released"
    ```

    - before you announce the release, be sure to send your announcement
      text to the [Biopython mailing list](mailto:biopython@biopython.org) for
      proof-reading/final corrections.
    - add to [main page](Main_Page "wikilink") and [downloads
      page](Download "wikilink") (through the wiki), make sure the links
      work.

23. Announcement:

    - post the announcement on the [www.open-bio.org](https://www.open-bio.org)
      blog (making sure to use the Biopython category which will update the
      [news page](News "wikilink") and [twitter](http://twitter.com/Biopython)
      via the news feed)
    - send an email to biopython-announce@biopython.org
      (see [mailing lists](Mailing_lists "wikilink"))
    - forward the email to Linux packagers e.g.
      debian-med@lists.debian.org

23. Conda-Forge should automatically open a pull request to update the
    package once it appears on PyPI. Check for a new pull request on
    [github.com/conda-forge/biopython-feedstock](https://github.com/conda-forge/biopython-feedstock)
    which once merged will upload the new release to [anaconda.org/conda-forge/biopython](https://anaconda.org/conda-forge/biopython)

Post release version bump
-------------------------

24. Bump version numbers again

    - Update `Bio/__init__.py` version
    - Biopython Tutorial - update the date/version line in the
      `Doc/Tutorial.tex` file
    - Make sure to commit the modified files to github.

    Include the suffix ``.dev0`` to indicate this is a development version
    e.g. if you have `__version__ = "1.68"`, make it `1.69.dev0`
