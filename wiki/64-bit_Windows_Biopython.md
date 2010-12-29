---
title: 64-bit Windows Biopython
permalink: wiki/64-bit_Windows_Biopython
layout: wiki
---

Compiling and using a 64 bit version of Biopython on Windows
============================================================

64 Bit Windows machines normally use a 32 bit version of
Biopython/Python/NumPy. This page is a scratchpad documenting how to
generate a 64 bit version.

Firstly, such a version is already available, courtesy of Christoph
Gohlke at [Unofficial Windows Binaries for Python Extension
Packages](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

To use a 64-bit version the following is needed: 64-bit Python and 64
bit NumPy. As of today, there is no official version 64 bit of NumPy.
Chris makes one available at the address above.

If you intend to compile Biopython from scratch then a **64 bit C
compiler is also needed** e.g.
[mingw64](http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/sezero_20101003/).

If you have both 32 bit and 64 bit versions of the above installed than
some care is needed to avoid mixups. Using the wrong compiler seems to
happen quite a bit (if you have both installed).

Using the mingw64 suite
-----------------------

There is a minor hurdle to sort before being able to use the mingw64
suite: Python requires msvcr90, but mingw 64 does not provide a stub
(libmsvcr90.a) for it. This problem doesn't exist in the 32 bit version
as the stub is provided. To get the stub to the following.

1.  Be sure to have a mingw64 version with gendef.exe. Some builds do
    not have it. Get one that does
2.  do copy
    C:\\Windows\\winsxs\\amd64\_microsoft.vc90.crt\_1fc8b3b9a1e18e3b\_9.0.21022.8\_none\_750b37ff97f4f68b\\msvcr90.dll .
3.  gendef.exe msvcr90.dll
4.  dlltool.exe -D msvcr90.dll -l libmsvcr90.a -d msvcr90.def
5.  copy libmsvcr90.a MINGW64\_PATH\\lib

TBC...
