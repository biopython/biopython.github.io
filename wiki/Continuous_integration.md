---
title: Continuous integration
permalink: wiki/Continuous_integration
layout: wiki
---

Currently all our continuous integration is done with free services
linked into the GitHub ecosystem. There are linked badges on our
[main repository's README page](https://github.com/biopython/biopython/blob/master/README.rst):

* [![Build Status](https://travis-ci.org/biopython/biopython.svg?branch=master)](https://travis-ci.org/biopython/biopython/branches)
  Linux testing with TravisCI, see <https://github.com/biopython/biopython/blob/master/.travis.yml>
* [![Build Status](https://img.shields.io/appveyor/ci/biopython/biopython/master.svg)](https://ci.appveyor.com/project/biopython/biopython/history)
  Windows testing with AppVeyor, see <https://github.com/biopython/biopython/blob/master/.appveyor.yml>

These automated tests in turn collect test coverage which is reported
on another separate free service linked into the GitHub ecosystem:

* [![Test Coverage](https://img.shields.io/codecov/c/github/biopython/biopython/master.svg)](https://codecov.io/github/biopython/biopython/)
  Test coverage with CodeCov.io

Historically Biopython used to run a [Buildbot](http://buildbot.net) server
hosted on an [Open Bioinformatics Foundation](https://www.open-bio.org)
server which displayed the results of nightly test jobs submitted to a
range of volunteer worker machines, including Linux, Windows and Mac OS
machines.


