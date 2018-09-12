---
title: How to build a release for Biopython.
permalink: wiki/Building_a_release
layout: wiki
---

Build Biopython in 22 easy steps!!

Setup required for a new release manager
----------------------------------------

The instructions below require that you have access to a few servers and
the code repository. When you start, be sure to have write access to:

1. [GitHub Biopython source code repository](https://github.com/biopython/biopython)
2. [GitHub biopython.org website repository](https://github.com/biopython/biopython.github.io)
3. [GitHub biopython.org/DIST/ repository](https://github.com/biopython/DIST)
4. [GitHub biopython wheel repository](https://github.com/biopython/biopython-wheels)
5. [OBF WordPress Blog](https://news.open-bio.org)
6. [Biopython on PyPI](https://pypi.python.org/pypi/biopython)

If you don't have any of the above, please ask.

The instructions proper
-----------------------

These instructions are for a Unix machine, with a Windows machine also
needed to test and prepare the Windows installers.

1. make sure I have the latest code:

   ``` bash
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

7. do last check to make sure things are checked in:

   ``` bash
   $ rm -r build
   $ rm Tests/*.pyc
   $ make clean -C Doc
   $ git status
   ```

8. build Biopython and do last regression test:

   ``` bash
   drevil:~biopython> python setup.py build
   drevil:~biopython> python setup.py test
   ```

   Ideally do this with a clean checkout on your Windows machine too.
   Assuming you have setup your compilers etc appropriately just do
   this:

   ```
   C:\python27\python setup.py build
   C:\python27\python setup.py test

   C:\python33\python setup.py build
   C:\python33\python setup.py test

   C:\python34\python setup.py build
   C:\python34\python setup.py test
   ```

   If you are using MinGW, do not forget to add `--compiler=mingw32`
   (or make it the default on distutils, see the step on building on
   Windows machines). Also If you are using a modern MinGW compiler,
   then distutils of Python will use an option (`-mno-cywgin`) that
   is deprecated and will break gcc. A possible solution is to
   [edit distutils](http://bugs.python.org/issue12641), but on
   recent Python (3.3 as tested) this seems to have been addressed.

   Running the tests simultaneously is risky, as two threads may both
   try to read/write to the same temp files.

9. check out clean version somewhere else:

   ``` bash
   drevil:~tmp1/> git clone git://github.com/biopython/biopython.git
   drevil:~tmp1/> cd biopython
   ```

10. make documentation PDF, text and HTML files in Doc:

    ``` bash
    drevil:~tmp1/biopython/> make -C Doc
    drevil:~tmp1/biopython/> make clean -C Doc
    ```

11. make the source distribution

    ``` bash
    drevil:~tmp1/biopython> python setup.py sdist --formats=gztar,zip
    ```

12. untar the file somewhere else

    ``` bash
    drevil:~tmp1/biopython/> cd ..
    drevil:~tmp1/> tar -xzvf biopython/dist/biopython-1.71.tar.gz
    drevil:~tmp1/> cd biopython-1.71
    ```

    - Check to make sure it includes the HTML and PDF files under Doc

13. make sure I can build and test it

    ``` bash
    drevil:~tmp1/biopython-1.71/> python setup.py build
    drevil:~tmp1/biopython-1.71/> python setup.py test
    drevil:~tmp1/biopython-1.71/> python setup.py install --prefix /tmp/test-install
    ```

    A typical source of failure here (on the tests) is the lack of example
    files being added to the source distribution: add them to `MANIFEST.in`

14. Update API documentation using Epydoc (this can often report otherwise overlooked
    errors).

    - If you haven't already, clone the ``DIST`` repository - otherwise first
      fetch any upstream changes.

    ``` bash
    $ cd ~/repositories
    $ git clone git@github.com:biopython/DIST.git
    $ cd ~/repositories/DIST
    ```

    - Remove version of the API documentation which we're going to replace:

    ``` bash
    $ git rm docs/api/*
    ```

    - Go to the `/tmp/test-install/lib/python2.7/site-packages` directory. This is the
      directory created under your source installation after the install step above.
      Running epydoc in your git tree works, but can miss some packages due to import
      errors.

    ``` bash
    $ cd /tmp/test-install/lib/python2.7/site-packages
    $ grep "__version__" Bio/__init__.py
    __version__ = "1.71"
    $ epydoc -v -o ~/repositories/DIST/docs/api/ -u http://biopython.org -n Biopython --docformat restructuredtext Bio BioSQL
    $ grep "__version__" ~/repositories/DIST/docs/api/Bio-pysrc.html
    <a name="L13"></a><tt class="py-lineno"> 13</tt>  <tt class="py-line"><tt class="py-name">__version__</tt> <tt class="py-op">=</tt> <tt class="py-string">"1.71"</tt> </tt>
    ```

    - Assuming no mismatches in version number (which would suggest epydoc is not looking
      at the new test installation), commit these new HTML files to the `gh-pages` branch:

    ``` bash
    $ cd ~/repositories/DIST
    $ git add docs/api/*
    $ git commit -m "epydoc for Biopython 1.68"
    ```

    - Update the tutorial too:

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

15. Back in the main repository, tag the release:

    ``` bash
    $ cd  .../tmp1/biopython/
    $ git tag biopython-171
    $ git push origin master --tags
    ```

16. Now we use https://github.com/biopython/biopython-wheels to build wheels,
    by updating the ``BUILD_COMMIT`` line in ``.travis.yml`` and ``appveyor.yaml``
    to the new release's git tag.

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

    Successful wheels will be on
    [on Rackspace](https://a365fff413fe338398b6-1c8a9b3114517dc5fe17b7c3f8c63a43.ssl.cf2.rackcdn.com/),
    download them from there to your ``~/repository/biopython/dist`` folder.
    We will upload these to PyPI later using Twine.

17. Remove any prior Biopython installations on your windows machine,
and confirm the Windows wheel file(s) work. Then copy them to your Linux
machine.

18. Upload the new release to the website via GitHub Pages `DIST` repository.

    ``` bash
    $ cp dist/biopython-1.68.* ~/repositories/DIST/
    # Also copy in the Windows files
    $ cd ~/repositories/DIST/
    $ git add biopython-1.71.*
    $ git commit biopython-1.71.* -m "Downloads for Biopython 1.71"
    $ shasum -a 256 biopython-1.71.*
    $ md5sum biopython-1.71.*
    $ git commit --amend # paste checksums into comment
    $ git push origin gh-pages
    ```

19. Upload to the python package index (except for beta/alpha level
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

20. Update the website and announce the release:

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
    $ git commit _config.yml wiki/Biopython.md wiki/Download.md -m "Biopython 1.68 released"
    ```

    - before you announce the release, be sure to send your announcement
      text to the [Biopython mailing list](mailto:biopython@biopython.org) for
      proof-reading/final corrections.
    - add to [main page](Main_Page "wikilink") and [downloads
      page](Download "wikilink") (through the wiki), make sure the links
      work
    - post the announcement on
      [news.open-bio.org](http://news.open-bio.org) (which will update the
      [news page](News "wikilink") and
      [twitter](http://twitter.com/Biopython) via the news feed)
    - add the new version to
      [RedMine](https://redmine.open-bio.org/projects/biopython)
    - send email to biopython-announce@biopython.org
      (see [mailing lists](Mailing_lists "wikilink"))
    - forward the email to Linux packagers e.g.
      debian-med@lists.debian.org

21. Conda-Forge should automatically open a pull request to update the
    package once it appears on PyPI. Check for a new pull request on
    [github.com/conda-forge/biopython-feedstock](https://github.com/conda-forge/biopython-feedstock)
    which once merged will upload the new release to [anaconda.org/conda-forge/biopython](https://anaconda.org/conda-forge/biopython)

22. Bump version numbers again

    - Update `Bio/__init__.py` version
    - Biopython Tutorial - update the date/version line in the
      `Doc/Tutorial.tex` file
    - Make sure to commit the modified files to github.

    Include the suffix ``.dev0`` to indicate this is a development version
    e.g. if you have `__version__ = "1.68"`, make it `1.69.dev0`
