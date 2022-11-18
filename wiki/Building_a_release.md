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

1. [OBF WordPress Blog](https://www.open-bio.org)
2. [Biopython on PyPI](https://pypi.python.org/pypi/biopython)
3. These repositories (e.g. via membership of the *Releases* team, or have someone ready to merge pull requests):

  - [Source code repository](https://github.com/biopython/biopython)
  - [Main website repository](https://github.com/biopython/biopython.github.io)
  - [Website DIST repository](https://github.com/biopython/DIST)
  - [Wheel building repository](https://github.com/biopython/biopython-wheels)

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
5. [hevea](http://hevea.inria.fr/), I am currently using version 2.32 of 2012-07-04

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
   $ git log --since="2020/05/25" --reverse --pretty="medium"
   ```

4. make sure `CONTRIB.rst` still current

5. make sure `setup.py` and `MANIFEST.in` are still up to date

   - Are there any new modules/files which should get installed?

6. bump version numbers and set the release data:

   - Biopython version - edit `Bio/__init__.py`
   - Biopython `NEWS.rst` - fill in the release date
   - Make sure to commit the modified files to github, e.g.

   ``` bash
   $ git commit Bio/__init__.py NEWS.rst -m "Call this Biopython 1.78"
   ```

7. do a final check to make sure things are checked in:

   ``` bash
   $ rm -r build
   $ rm Tests/*.pyc
   $ make clean -C Doc
   $ git status
   ```

8. build Biopython and run a final regression test:

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
    drevil:~tmp1/> tar -xzvf biopython/dist/biopython-1.78.tar.gz
    drevil:~tmp1/> cd biopython-1.78
    ```

    Check to make sure it includes the HTML and PDF files under Doc

14. make sure I can build and test it

    ``` bash
    drevil:~tmp1/biopython-1.78/> python setup.py build
    drevil:~tmp1/biopython-1.78/> python setup.py install --prefix /tmp/test-install
    drevil:~tmp1/biopython-1.78/> cd Tests && python run_tests.py
    ```

    A typical source of failure here (on the tests) is the lack of example
    files being added to the source distribution: add them to `MANIFEST.in`

Checking the compiled documentation
-----------------------------------

15. Since Biopython 1.74, Sphinx has handled the API documentation via continuous
    integration, but you still have to update the "latest" symlink.

    ``` bash
    $ cd ~/repositories/docs/
    $ git fetch origin
    $ git checkout gh-pages  # should only be this one branch
    $ git checkout rebase origin/gh-pages  # get any changes
    $ rm latest
    $ ln -s 1.78 latest
    $ git commit latest -m "Update 'latest' symlink to point at 1.78"
    $ git push origin gh-pages
    ```

16. Update Tutorial and PDB FAQ on the website manually.

      - Update the Tutorial and PDB FAQ:

      ``` bash
      $ cd ~/repositories/DIST/docs/tutorial/
      $ cp ../../../biopython/Doc/biopdb_faq.pdf .
      $ cp ../../../biopython/Doc/Tutorial.html Tutorial-1.80.html
      $ cp ../../../biopython/Doc/Tutorial.pdf Tutorial-1.80.pdf
      $ rm Tutorial.html Tutorial.pdf
      $ ln -s Tutorial-1.80.html Tutorial.html
      $ ln -s Tutorial-1.80.pdf Tutorial.pdf
      $ git add Tutorial-1.80.html Tutorial-1.80.pdf
      $ git commit Tutorial-1.80.html Tutorial-1.80.pdf Tutorial.html Tutorial.pdf biopdb_faq.pdf -m "Tutorial and FAQ for Biopython 1.80"
      ```

      - Push this to GitHub Pages to update the website:

      ``` bash
      $ git push origin gh-pages
      ```

      - Check this is live at <http://biopython.org/DIST/docs/api/Bio-module.html>,
        <http://biopython.org/DIST/docs/tutorial/Tutorial.html>,
        <http://biopython.org/DIST/docs/tutorial/Tutorial.pdf>, and
        <http://biopython.org/DIST/docs/tutorial/biopdb_faq.pdf>

Making wheels
-------------

17. Now we use https://github.com/biopython/biopython-wheels to build wheels,
    by updating the ``git checkout`` line in ``.github/workflows/cibuildwheel.yml``
    to the new release's commit hash (which all being well will get a git tag).

    ``` bash
    $ cd ~/repositories
    $ git clone git@github.com:biopython/biopython-wheels.git
    $ cd biopython-wheels/
    $ git submodule update --init
    $ emacs .github/workflows/cibuildwheel.yml  # update git checkout line
    $ git commit .github/workflows/cibuildwheel.yml -m "Build Biopython 1.xx"
    $ git push origin master
    ```

    Check the wheels build on the [GitHub Actions runs](https://github.com/biopython/biopython-wheels/actions).

    You don't seem to need to update the ``biopython`` git submodule, but if you
    need to this seems to work.

    ``` bash
    $ git submodule foreach git pull origin master
    $ git commit -a -m "Update submodules"
    $ git push origin master
    ```

18. Successful wheels will in an ``artifact.zip`` file available in the footer of the
    run via [GitHub Actions runs](https://github.com/biopython/biopython-wheels/actions).
    Download this and unzip to your ``~/repository/biopython/DIST/`` folder.
    We will upload these to PyPI later using Twine.

19. If you have a Windows machine, remove any prior Biopython installations,
    and confirm the Windows wheel file(s) work.

Tagging the release, and uploading
----------------------------------

20. Back in the main repository, tag the release:

    ``` bash
    $ cd  .../tmp1/biopython/
    $ git tag biopython-178
    $ git push origin master --tags
    ```

21. Upload the new release tar-ball and zip to the website via GitHub Pages `DIST` repository.

    ``` bash
    $ cp dist/biopython-1.78.* ~/repositories/DIST/
    $ cd ~/repositories/DIST/
    $ git add biopython-1.78.*
    $ git commit biopython-1.78.* -m "Downloads for Biopython 1.78"
    $ shasum -a 256 biopython-1.78.*
    $ md5sum biopython-1.78.*
    $ git commit --amend # paste checksums into comment
    $ git push origin gh-pages
    ```

22. Upload to the python package index (except for beta/alpha level releases):

    ``` bash
    $ cd  ~/repositories/biopython/
    $ pip install twine
    $ twine upload dist/biopython-1.78.tar.gz
    $ twine upload dist/biopython-1.78-*.whl
    ```

    - You need to have a login on pypi and be registered with Biopython to be
      able to upload the new version

    - Check this is live at <https://pypi.python.org/pypi/biopython/>

23. Update the website:

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
    $ git commit _config.yml wiki/Biopython.md wiki/Download.md -m "Biopython 1.78 released"
    ```

    - before you announce the release, be sure to send your announcement
      text to the [Biopython mailing list](mailto:biopython@biopython.org) for
      proof-reading/final corrections.
    - add to [main page](Main_Page "wikilink") and [downloads
      page](Download "wikilink") (through the wiki), make sure the links
      work.

24. Announcement:

    - post the announcement on the [www.open-bio.org](https://www.open-bio.org)
      blog (making sure to use the Biopython category which will update the
      [news page](News "wikilink") and [twitter](http://twitter.com/Biopython)
      via the news feed)
    - send an email to biopython-announce@biopython.org
      (see [mailing lists](Mailing_lists "wikilink"))
    - forward the email to Linux packagers e.g.
      debian-med@lists.debian.org

25. Conda-Forge should automatically open a pull request to update the
    package once it appears on PyPI. Check for a new pull request on
    [github.com/conda-forge/biopython-feedstock](https://github.com/conda-forge/biopython-feedstock)
    which once merged will upload the new release to [anaconda.org/conda-forge/biopython](https://anaconda.org/conda-forge/biopython)

Post release version bump
-------------------------

26. Bump version numbers again

    - Update `Bio/__init__.py` version
    - Start entry in `NEWS.rst` for next version
    - Make sure to commit the modified files to github.

    Include the suffix ``.dev0`` to indicate this is a development version
    e.g. if you had `__version__ = "1.78"`, make it `1.79.dev0`
