---
title: Biopython.org - README
layout: default
---

# README

### Few things to note:
- The git repository at <https://github.com/biopython/biopython.github.io/>
produces the Biopython website at <https://biopython.org> using GitHub Pages (rendered using Jekyll, see <https://help.github.com/articles/using-jekyll-with-pages/> for details).

- Almost all of the content is under the URL prefix ``/wiki/`` because this was based on an automated conversion of the old MediaWiki website, using <https://github.com/peterjc/mediawiki_to_git_md> to turn all the changes in the MediaWiki XML export file into into markdown (using pandoc) as a git repository.

- The old MediaWiki usernames were manually mapped to GitHub accounts. About 20 accounts on the wiki (mostly single contributions) could not be identified, but the old username is still logged in the git commits.

- Note the website content under ``/DIST/`` is hosted in a separate GitHub Pages project repository <https://github.com/biopython/DIST> covering the Biopython releases and assorted documentation files.

### To contribute to this site:
- If you haven't already, begin by creating an account for yourself [here](http://lists.open-bio.org/mailman/listinfo/biopython-dev/) at the Biopython Developers List.
- Introduce yourself to the Biopython community by sending an email to <biopython-dev@biopython.org>.
- If you haven't already, create a [Github](https://github.com/) account [here](https://github.com/join?source=header-home).
- Then fork a copy of <https://github.com/biopython/biopython.github.io> repository which holds the source for this site.
- The above step should create a copy of the repository under your Github user account.
- Clone this copy of the repository to your local machine.
- Follow the instructions (specific to your OS) to setup a local *Jekyll 3* environment (see <https://jekyllrb.com/docs/installation/>) so you can test your changes locally.
- Then create a descriptive local branch and do your magic - fix a bug/outstanding issue and/or improve/implement a feature.
- Test your changes locally and push to your local repository.
- Finally submit a pull request!