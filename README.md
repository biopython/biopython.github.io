---
title: README
---

This is an attempt to automatically convert a Mediawiki XML export
from http://biopython.org into markdown (using pandoc) as a git
repository to be hosted using GitHub Pages (rendered using Jekyll).
https://help.github.com/articles/using-jekyll-with-pages/

The conversion is via a Python script (calling pandoc and git), see:
https://github.com/peterjc/mediawiki_to_git_md

I'm using a manually compiled table mapping MediaWiki usernames
to GitHub accounts - if I have mis-identified you, please email
me during this testing period and I'll remove the false mapping.
I failed to match about 20 accounts on the wiki (mostly single
contributions).

This site ought to be viewable via https://peterjc.github.io/
and this page as *html* at https://peterjc.github.io/README.html
and https://github.com/peterjc/peterjc.github.io/blob/master/README.md
in the original source *markdown* view on GitHub.

Branches: GitHub pages will automatically show the ``master`` branch
on https://peterjc.github.io/ which I am therefore using for live
testing of the automated imports. This means I will regularly
re-write git history with replacement ``master`` branches.
Each time I will return to the ``pre_auto_import`` branch.
