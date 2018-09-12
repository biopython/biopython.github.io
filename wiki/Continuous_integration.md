---
title: Continuous integration
permalink: wiki/Continuous_integration
layout: wiki
---

This is a temporary page documenting efforts to provide a continuous
integration platform for Biopython. This effort is ad-hoc for now and
this page will probably be deleted in the future.

Are you interested in helping testing Biopython for your preferred
platform? Then read on.

We have a [Buildbot](http://buildbot.net) system up to do integration
testing. The information provided here is split in two main sections:

1.  Documentation for people maintaining volunteer testing machines. Of
    interest to users donating resources for testing.

<!-- -->

1.  Information regarding server configuration. Of interest to both
    developers interested in the integration testing status and
    administrators of the buildbot central server.

Checking the integration status
===============================

Currently you have 3 options

1.  Visit our [buildbot web site](http://testing.open-bio.org/biopython)
2.  Subscribe to the integration [RSS
    feed](http://testing.open-bio.org/biopython/rss)
3.  Subscribe to the [Biopython mailing list](mailto:biopython@biopython.org)
    where important errors are automatically forward (under construction)

Volunteer integration testing
=============================

General instructions
--------------------

You will need to install buildbot on your system (*ONLY* the slave
part). Please see below more specific instructions for your operating
system.

You will also need git installed.

After that you should contact the Biopython development mailing list
saying what you are willing to test (platform, python version,
applications) to get a username and a password. With your username and
password you can now configure your buildslave. Something like this:

```bash
buildslave create-slave your_slave_directory testing.open-bio.org:9989
username password
```

Then you should start your slave with

``` bash
buildslave start your_slave_directory
```

Remember to have the external applications that you wish to test (e.g.
blast) correctly installed, if you want to test all of Biopython you
might want to have a look at the [List of applications executed via
Biopython](List_of_applications_executed_via_Biopython "wikilink") and
install them all. You should also have all dependencies needed by
Biopython installed (e.g. NumPy).

**Please read the security notes at the end.**

Linux
-----

Installing buildbot-slave is easy. There is just a pitfall: in some
distributions buildbot and buildbot-slave are the same program (just
called buildbot). In this case your create-slave command is handled by
buildbot (not buildslave)

On linux you should make sure that the version of python you want to
test is available by the name pythonM.m. E.g, if you offer to test for
python 2.6, your python binary should be available by the name python2.6
This is the case in most distributions.

Windows
-------

After you have your version of Python installed, you will need to
install buildbot. We recommend using Python 2.7. `easy_install`
(setuptools) will allow you to install all dependencies except `pywin32`
and maybe `twister` (just use the binary installer for the packages).

See: <http://buildbot.net/trac/wiki/RunningBuildbotOnWindows>

If you install buildbot slave as a service (recommended), you might want
to run the python setup script on the link above (buildbot_service.py)
on a shell (cmd.exe) with Administrator privileges (right-click on
cmd.exe and chose "Run as Administrator").

Mac
---

If when running `buildslave start` the slave doesn't seem to start, try
adding python to the list of applications allowed to accept incoming
connections in the Mac OS X advanced firewall settings.

Jython
------

Server configuration
====================

Buildbot has a server/client architecture where a central server
schedules builds that are actually build and tested by the clients. The
clients are typically volunteer machines. The main requirement of a
public accessible server running Buildbot (`twister` based).

One note for people trying to install a Builbot server: It seems that
Buildbot configuration files change a bit from version to version. Be
sure to use the correct documentation for your version. Our examples are
for 0.8.1p1.

Workflow
--------

Understanding [Buildbot's
architecture](http://buildbot.net/buildbot/docs/current/full.html) might
be important to understand the next paragraphs.

The starting point for any continuous integration will be any change
made to the main repository. In Biopython's case, this is based on
github.

GitHub can inform buildbot that a change was committed using [Post
receive service hooks](http://help.github.com/post-receive-hooks/)
(Admin&gt;Service Hooks on the Biopython GitHub interface). This means
that whenever there is a change, GitHub can POST that to a list of
specified URLs. One of those URLs will be a script that will call the
buildbot master, informing that there are changes to the source.
Buildbot can now act.

Buildbot is split in two parts: The master and (potentially several)
slaves.

The master:

1.  Is informed of repository changes
2.  Makes decisions on what should be tested and when (schedules work
    for the slaves)
3.  Informs the outside world of results (via email, web, irc, ...)
4.  Has to be visible as a server with a public address

The slaves:

1.  Do the actual integration testing
2.  Can be donated by volunteers that want to help test for a specific
    platform
3.  Can be run anywhere
4.  We need at least 3 (one per major OS)

Before Buildbot
---------------

Before buildbot we have git. We need to inform buildbot of git changes.

[//]: # (BROKEN Recommended read is: [github/buildbot integration](http://www.apparatusproject.org/blog/2009/06/github-and-buildbot-continuous-integration/).)

While, from an architectural point of view it might seem
convoluted/complex. The implementations is actually quite simple: a
simple entry in github, a small cgi script and configuring buildbot to
accept github's source.

The only problem is that Buildbot has no native support for GitHub (at
least to be informed of changes -- buildbot can download GitHub code),
so the solution is the aforementioned cgi script. For reference the
correct source in buildbot is the generic PBChangeSource (see below)

Builbot: general considerations
-------------------------------

Typical usage. Command line (configure, start, reconfig). Different
versions, different configurations

Buildbot master
---------------

Buildbot master is arguably the most complex part for configuration.
Most decisions are made at this level. The suggestions below are just
that, suggestions. A possible starting point.

The configuration is divided in:

### Slaves

The list of slaves. Only listed slaves will be distributed work. Each
slave has a password.

### Sources

The build sources. In our case Biopython's github. As there is no native
support, we use the standard PBChangeSource.

### Schedulers

Scheduling builders is where many decisions will fall upon. There are
many options like a fast scheduling to test very recent changes with
simple updates as soon as possible, or a daily build, doing a complete
download.

Schedulers also seem to change from version to version of buildbot.

Currently there is a single periodic scheduler (once a day).

### Builders

A builder is responsible for testing the source in a certain
environment. The environment can be defined by things like the
Python/Jython version, the OS, but also if it is a fast or slow build.
The number of builders can explode quite easily: lets say 5 Python
versions (2.5, 2.6, 2.7, 3.1 plus Jython 2.5.2) times 3 OSes times 2
build types (fast and slow) and we are at 30 builds.

The two most important parameters of a build are: the slave that
implements the build and the list of steps the slave has to do. There
are typically:

1.  Get the source (either a light update or a full download), for this
    buildbot has GitHub support.
2.  Compile the source
3.  Run the tests

### Status targets

Targets are systems that relay information on the test outcomes. Typical
targets are:

1.  Buildbot webserver. A web address where developers and users can
    check the results of previous builds
2.  Email addresses. To be warned of outcomes of builders
3.  IRC channels
4.  ...

There is a main buildbot webserver. There are two types of users: the
general audience that can check the results. Admins which can stop the
server, restart builds, etc...

### Identity

### Debugging

Security issues
===============

Slaves are identified by the means of a password. If a password is
known, there is not much problem: The only thing a rogue slave can do is
provide back erroneous integration information.

A much more serious issue is a rogue (cracked) server. Servers can
instruct slaves to run arbitrary code, thus a compromised server can
also compromise the slaves. As such, volunteers running slaves should
probably establish separate accounts for the slaves and restrict the
access of such accounts to Biopython requirements only (preferably no
write access outside the slave scratch area).
