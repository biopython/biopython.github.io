---
title: GitUsage
permalink: wiki/GitUsage
layout: wiki
---

These are (draft) general guidelines for Biopython development using
git. We're still working on the finer details etc.

This document is meant as an outline of the way Biopython is developed.
It should include all essential technical information as well as typical
procedures and usage scenarios. It should be helpful for core
developers, potential code contributors, testers and everybody
interested in Biopython code.

<b>This version is an unofficial draft and is subject to change.</b>

Relevance
=========

If you just want to grab the latest (not yet officially released)
Biopython from our repository, see our [source code
page](SourceCode "wikilink"). This page is about actually using git for
tracking changes.

If you have found a problem with Biopython, and think you know how to
fix it, then we suggest following the simple route of filing a
bug and describe your
fix. Ideally, you would upload a patch file showing the differences
between the latest version of Biopython (from our repository) and your
modified version. Working with the command line tools *diff* and *patch*
is a very useful skill to have, and is almost a precursor to working
with a version control system.

You shouldn't go to the trouble of creating your own git fork unless you
are intending to make more than a simple one off contribution.

Technicalities
==============

This section describes technical introduction into git usage including
required software and integration with Github. If you want to start
contributing to Biopython, you definitely need to install git and learn
how to obtain a branch of Biopython. If you want to share your changes
easily with others, you should also sign up for a Github account and
read the corresponding section of the manual. Finally, if you are
engaged in one of the collaborations on experimental Biopython modules,
you should look also into code review and branch merging.

Installing Git
--------------

You will need to install Git on your computer. [Git](http://git-scm.com/)
is available for all major operating systems. Please use the appropriate
installation method as described below.

### Linux

Git is now packaged in all major Linux distributions, you should find it
in your package manager.

#### Ubuntu/Debian

You can install Git from the `git-core` package. e.g.,

``` bash
sudo apt-get install git-core
```

You'll probably also want to install the following packages: `gitk`,
`git-gui`, and `git-doc`

#### Redhat/Fedora/Mandriva

git is also packaged in rpm-based linux distributions.

``` bash
yum install gitk
```

should do the trick for you in any recent fedora/mandriva or
derivatives

### Mac OS X

Download the `.dmg` disk image from
<http://code.google.com/p/git-osx-installer/>

### Windows

Download the official installers from
[Windows installers](https://git-scm.com/download/win)

Testing your git installation
-----------------------------

If your installation succeeded, you should be able to run

``` bash
git --help
```

in a console window to obtain information on git usage. If this fails,
you should refer to git
[documentation](https://git-scm.com/doc) for troubleshooting.

Creating a GitHub account (Optional)
------------------------------------

Once you have Git installed on your machine, you can obtain the code and
start developing. Since the code is hosted at GitHub, however, you may
wish to take advantage of the site's offered features by signing up for
a GitHub account. While a GitHub account is completely optional and not
required for obtaining the Biopython code or participating in
development, a GitHub account will enable all other Biopython developers
to track (and review) your changes to the code base, and will help you
track other developers' contributions. This fosters a social,
collaborative environment for the Biopython community.

If you don't already have a GitHub account, you can create one
[here](https://github.com/join).
Once you have created your account, upload an SSH public key by clicking
on '[SSH and GPG keys](https://github.com/settings/keys)' after logging in. For more
information on generating and uploading an SSH public key, see [this
GitHub guide](https://help.github.com/en/articles/connecting-to-github-with-ssh).

Working with the source code
============================

In order to start working with the Biopython source code, you need to
obtain a local clone of our git repository. In git, this means you will
in fact obtain a complete clone of our git repository along with the
full version history. Thanks to compression, this is not much bigger
than a single copy of the tree, but you need to accept a small overhead
in terms of disk space.

There are, roughly speaking, two ways of getting the source code tree
onto your machine: by simply "cloning" the repository, or by "forking"
the repository on GitHub. They're not that different, in fact both will
result in a directory on your machine containing a full copy of the
repository. However, if you have a GitHub account, you can make your
repository a public branch of the project. If you do so, other people
will be able to easily review your code, make their own branches from it
or merge it back to the trunk.

Using branches on Github is the preferred way to work on new features
for Biopython, so it's useful to learn it and use it even if you think
your changes are not for immediate inclusion into the main trunk of
Biopython. But even if you decide not to use github, you can always
change this later (using the .git/config file in your branch.) For
simplicity, we describe these two possibilities separately.

Cloning Biopython directly
--------------------------

Getting a copy of the repository (called "cloning" in Git terminology)
without GitHub account is very simple:

``` bash
git clone https://github.com/biopython/biopython.git
```

This command creates a local copy of the entire Biopython repository on
your machine (your own personal copy of the official repository with its
complete history). You can now make local changes and commit them to
this local copy (although we advise you to use named branches for this,
and keep the master branch in sync with the official Biopython code).

If you want other people to see your changes, however, you must publish
your repository to a public server yourself (e.g. on GitHub).

Forking Biopython with your GitHub account
------------------------------------------

If you are logged in to GitHub, you can go to the Biopython repository
page:

[https://github.com/biopython/biopython/tree/master](https://github.com/biopython/biopython/tree/master)

and click on a button named 'Fork'. This will create a fork (basically a
copy) of the official Biopython repository, publicly viewable on GitHub,
but listed under your personal account. It should be visible under a URL
that looks like this:

https://github.com/yourusername/biopython

Since your new Biopython repository is publicly visible, it's considered
good practice to change the description and homepage fields to something
meaningful (i.e. different from the ones copied from the official
repository).

If you haven't done so already, setup an SSH key and [upload it to
github](https://github.com/settings/keys) for
authentication.

Now, assuming that you have git installed on your computer, execute the
following commands locally on your machine. This "url" is given on the
GitHub page for your repository (if you are logged in):

``` bash
git clone git@github.com:yourusername/biopython.git
```

Where `yourusername`, not surprisingly, stands for your GitHub username.
You have just created a local copy of the Biopython repository on your
machine.

You may want to also link your branch with the official distribution
(see below on how to keep your copy in sync):

``` bash
git remote add upstream https://github.com/biopython/biopython.git
```

To add additional contributors to your repository on GitHub (i.e. people
you want to be able to commit to it), select 'edit' and then add them to
the 'Repository Collaborators' section. You will need to know their
username on GitHub.

If you haven't already done so, tell git your name and the email address
you are using on GitHub (so that your commits get matched up to your
GitHub account). For example,

``` bash
git config --global user.name "David Jones"
git config --global user.email "d.jones@example.com"
```

Setting up a coding-style checker
---------------------------------

Biopython tries to follow the coding conventions laid out in PEP8 and PEP257.

Before starting to work on the code, we ask you to install some tools for automated
checks. This includes a git pre-commit hook so that each of your commits (see below)
will automatically be checked for violations of Biopython's agreed coding style.
Commits with violations will be blocked. Thus you ensure that a later submission to
Biopython (a pull request, see below) will not be stopped by our automatic online
style-checks.

See the [CONTRIBUTING.rst](https://github.com/biopython/biopython/blob/master/CONTRIBUTING.rst) file for more.

Making changes locally
----------------------

Now you can make changes to your local repository - you can do this
offline, and you can commit your changes as often as you like. In fact,
you should commit as often as possible, because smaller commits are much
better to manage and document.

First of all, create a new branch to make some changes in, and switch to
it:

``` bash
git branch demo-branch
git checkout demo-branch
```

To check which branch you are on, use:

``` bash
git branch
```

Let us assume you've made changes to the file Bio/x.py. Try this:

``` bash
git status
```

So commit this change you first need to explicitly add this file to your
change-set:

``` bash
git add Bio/x.py
```

and now you commit:

``` bash
git commit -m "added feature Y in Bio.x"
```

Your commits in Git are local, i.e. they affect only your working branch
on your computer, and not the whole Biopython tree or even your fork on
GitHub. You don't need an internet connection to commit, so you can do
it very often.

Pushing changes to Github
-------------------------

If you are using Github, and you are working on a clone of your own
branch, you can very easily make your changes available for others.

Once you think your changes are stable and should be reviewed by others,
you can push your changes back to the GitHub server:

``` bash
git push origin demo-branch
```

*This will not work if you have cloned directly from the official
Biopython branch, since only the core developers will have write access
to the main repository.*

Merging upstream changes
------------------------

We recommend that you don't actually make any changes to the **master**
branch in your local repository (or your fork on github). Instead, use
named branches to do any of your own work. The advantage of this
approach it is the trivial to pull the upstream **master** (i.e. the
official Biopython branch) to your repository.

Assuming you have issued this command (you only need to do this once):

``` bash
git remote add upstream https://github.com/biopython/biopython.git
```

Then all you need to do is:

``` bash
git checkout master
git pull upstream master
```

Provided you never commit any change to your local **master** branch,
this should always be a simple *fast forward* merge without any
conflicts. You can then deal with merging the upstream changes from your
local master branch into your local branches (and you can do that
offline).

If you have your repository hosted online (e.g. at github), then push
the updated master branch there:

``` bash
git push origin master
```

Submitting changes for inclusion in Biopython
---------------------------------------------

If you think you changes are worth including in the main Biopython
distribution, then file an (enhancement) bug on our bug
tracker, and include a
link to your updated branch (i.e. your branch on GitHub, or another
public Git server). You could also attach a patch to the bug. If the
changes are accepted, one of the Biopython developers will have to check
this code into our main repository.

On GitHub itself, you can inform keepers of the main branch of your
changes by sending a 'pull request' from the main page of your branch.
Once the file has been committed to the main branch, you may want to
delete your now redundant bug fix branch on GitHub. Branches can be
deleted by selecting 'edit' and then 'delete repository' from the bottom
of the edit page.

If other things have happened since you began your work, it may require
merging when applied to the official repository's master branch. In this
case we might ask you to help by rebasing your work:

``` bash
git fetch upstream
git checkout demo-branch
git rebase upstream/master
```

Hopefully the only changes between your branch and the official repository's
master branch are trivial and git will handle everything automatically.
If not, you would have to deal with the clashes manually. If this works,
you can update the pull request by replacing the existing (pre-rebase)
branch:

``` bash
git push origin demo-branch --force
```

If however the rebase does not go smoothly, give up with the following command
(and hopefully the Biopython developers can sort out the rebase or merge for you):

``` bash
git rebase --abort
```

Evaluating changes
------------------

Since git is a fully distributed version control system, anyone can
integrate changes from other people, assuming that they are using
branches derived from a common root. This is especially useful for
people working on new features who want to accept contributions from
other people.

This section is going to be of particular interest for the Biopython
core developers, or anyone accepting changes on a branch.

For example, suppose Eric has some interesting changes on his public
repository:

https://github.com/etal/biopython.git

You must tell git about this by creating a reference to this remote
repository:

``` bash
git remote add eric https://github.com/etal/biopython.git
```

Now we can fetch *all* of Eric's public repository with one line:

``` bash
git fetch eric
remote: Counting objects: 138, done.
remote: Compressing objects: 100% (105/105), done.
remote: Total 105 (delta 77), reused 0 (delta 0)
Receiving objects: 100% (105/105), 27.53 KiB, done.
Resolving deltas: 100% (77/77), completed with 24 local objects.
From https://github.com/etal/biopython
 * [new branch]      bug2754    -> eric/bug2754
 * [new branch]      master     -> eric/master
 * [new branch]      pdbtidy    -> eric/pdbtidy
 * [new branch]      phyloxml   -> eric/phyloxml
```

Now we can run a diff between any of our own branches and any of Eric's
branches. You can list your own branches with:

``` bash
git branch
* master
  ...
```

Remember the asterisk shows which branch is currently checked out.

To list the remote branches you have setup:

``` bash
git branch -r
 eric/bug2754
 eric/master
 eric/pdbtidy
 eric/phyloxml
 upstream/master
 origin/HEAD
 origin/master
 ...
```

For example, to show the difference between your **master** branch and
Eric's **master** branch:

``` bash
git diff master eric/master
...
```

If you are both keeping your **master** branch in sync with the upstream
Biopython repository, then his **master** branch won't be very
interesting. Instead, try:

``` bash
git diff master eric/pdbtidy
...
```

You might now want to merge in (some) of Eric's changes to a new branch
on your local repository. To make a copy of the branch (e.g. pdbtidy)
in your local repository, type:

``` bash
git checkout --track eric/pdbtidy
```

If Eric is adding more commits to his remote branch and you want to update
your local copy, just do:

``` bash
git checkout pdbtidy  # if you are not already in branch pdbtidy
git pull
```

If you later want to remove the reference to this particular branch:

``` bash
git branch -r -d eric/pdbtidy
Deleted remote branch eric/pdbtidy (79b5974)
```

Or, to delete the references to all of Eric's branches:

``` bash
git remote rm eric
git branch -r
  upstream/master
  origin/HEAD
  origin/master
  ...
```

Alternatively, from within GitHub you can use the fork-queue to cherry
pick commits from other people's forked branches. See [this GitHub blog
post](https://github.com/blog/270-the-fork-queue) for details. While this
defaults to applying the changes to your current branch, you would
typically do this using a new integration branch, then fetch it to your
local machine to test everything, before merging it to your main branch.

Committing changes to main branch
=================================

This section is intended for Biopython developers, who are allowed to
commit changes to the Biopython main "official" branch. It describes the
typical activities, such as merging contributed code changes both from
git branches and patch files.

Prerequisites
-------------

Currently, the main Biopython branch is hosted on github. In order to
make changes to the main branch you need a GitHub account and you need
to be added as a collaborator to the Biopython account. This needs to be
done only once. If you have a GitHub account, but you are not yet a
collaborator and you think you should be (for example, you had a cvs
account on open-bio server): ask Peter to be added (this is meant for
regular contributors, so in case you have only a single change to make,
please consider submitting your changes through one of developers).

Once you are a collaborator, you can pull Biopython official branch
using the private url. If you want to make a new repository (linked to
the main branch), you can just clone it:

``` bash
git clone git@github.com:biopython/biopython.git
```

It creates a new directory "biopython" with a local copy of the official
branch. It also sets the "origin" to the GitHub copy This is the
recommended way (at least for the beginning) as it minimizes the risk of
accidentally pushing changes to the official GitHub branch.

Alternatively, if you already have a working git repo (containing your
branch and your own changes), you can add a link to the official branch
with the git "remote command"... but we'll not cover that here.

In the following sections, we assume you have followed the recommended
scenario and you have the following entries in your .git/config file:

```
[remote "origin"]
       url = git@github.com:biopython/biopython.git

[branch "master"]
       remote = origin
```

Committing a patch
------------------

If you are committing from a patch, it's also quite easy. First make
sure you are up to date with official branch:

``` bash
git checkout master
git pull origin
```

Then do your changes, i.e. apply the patch:

``` bash
patch -r someones_cool_feature.diff
```

If you see that there were some files added to the tree, please add them
to git:

``` bash
git add Bio/Tests/some_new_file
```

Then make a commit (after adding files):

``` bash
git commit -a -m "committed a patch from a kind contributor adding feature X"
```

After your changes are committed, you can push to github:

``` bash
git push origin
```

Committing from someone's git branch
------------------------------------

Assume you want to merge changes someone has committed to a git
repository which was at some point cloned from the official Biopython
branch. He needs to make his repository available to you (read-only) by
giving you a URL. Typically this will be on GitHub (but it may be any
public git url). Let us assume that the url is (which happens to be my
clone of Biopython):

https://github.com/barwil/biopython.git

First, you need to get the code from this repository:

``` bash
git remote add Bartek https://github.com/barwil/biopython.git
git fetch Bartek
```

Then you can see what branches are there:

``` bash
git branch -r
 Bartek/master
 Bartek/motif_docs
 Bartek/test-branch
```

Let's say you want to merge changes from test-branch. You need to make
sure you are up to date with the official branch:

``` bash
git checkout master
git pull origin
```

And then you can do the actual merge:

``` bash
git pull Bartek test-branch
```

And (assuming you are OK with the results of git diff and git status),
you can push to the public repository on GitHub (please don't try that
with this exemplary data):

``` bash
git push origin
```

After you're done, you can remove the reference to the remote repo:

``` bash
git remote rm Bartek
```

Tagging the official branch
---------------------------

If you want to put tag on the current Biopython official branch (this is
usually done to mark a new release), you need to follow these steps:

First make sure you are up to date with official branch:

``` bash
git checkout master
git pull origin
```

Then add the actual tag:

``` bash
git tag new_release
```

And push it to github:

``` bash
git push --tags origin master
```

Additional Resources
====================

There are a lot of different nice guides to using Git on the web:

-   [Understanding Git
    Conceptually](https://www.sbf5.com/~cduan/technical/git/)
-   [git ready: git tips](http://gitready.com/)
-   <http://http://cheat.errtheblog.com/s/git>
-   <https://docs.scipy.org/doc/numpy-1.15.1/dev/gitwash/development_workflow.html> Numpy is also
    evaluating git
-   <https://github.github.com/training-kit/downloads/github-git-cheat-sheet>
-   <https://lab.github.com/courses>
-   [Pro Git](https://git-scm.com/book/en/v2)

