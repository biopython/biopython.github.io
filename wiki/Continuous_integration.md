---
title: Continuous integration
permalink: wiki/Continuous_integration
layout: wiki
---

This is a temporary page documenting efforts to provide a continuous
integration platform for Biopython. This effort is ad-hoc for now and
this page will probably be deleted in the future.

Currently using [Buildbot](http://buildbot.net). The main caveat is the
requirement of a public accessible server running buildbot (twister
based)

Workflow
--------

[Buildbot's
architecture](http://buildbot.net/buildbot/docs/current/full.html) might
be important to understand the next paragraphs.

The starting point for any continuous integration will be any change
made to the main repository. In Biopython's case, this is based on
github.

github can inform buildbot that a change was committed using [Post
receive service hooks](http://help.github.com/post-receive-hooks/)
(Admin&gt;Service Hooks on the Biopython github interface). This means
that whenever there is a change, github can POST that to a list of
specified URLs. One of those URLs will be a script that will call the
buildbot master, informing that there are changes to the source.
Buildbot can now act.

[github/buildbot
integration](http://www.apparatusproject.org/blog/2009/06/github-and-buildbot-continuous-integration/)
