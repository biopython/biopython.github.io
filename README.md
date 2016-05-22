---
title: Biopython.org website README
layout: default
---

The git repository at <https://github.com/biopython/biopython.github.io/>
produces the Biopython website at <https://biopython.org>
using GitHub Pages (rendered using Jekyll, see
<https://help.github.com/articles/using-jekyll-with-pages/> for details).

Almost all of the content is under the URL prefix ``/wiki/`` because
this was based on an automated conversion of the old MediaWiki website,
using <https://github.com/peterjc/mediawiki_to_git_md> to turn all
the changes in the MediaWiki XML export file into into markdown (using
pandoc) as a git repository.

The old MediaWiki usernames were manually mapped to GitHub accounts.
About 20 accounts on the wiki (mostly single contributions) could not
be identified, but the old username is still logged in the git commits.

Note the website content under ``/DIST/`` is hosted in a separate
GitHub Pages project repository <https://github.com/biopython/DIST>
covering the Biopython releases and assorted documentation files.
