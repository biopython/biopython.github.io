---
title: PDBParser
permalink: wiki/PDBParser
layout: wiki
---

This is a draft page for the PDBParser class.

Benchmark
=========

A performance benchmark of the parser was carried out to evaluate wether
the development of new features degraded the overall parsing speed.

Datasets
--------

[CATH Domain
Collection](http://release.cathdb.info/v3.4.0/CathDomainList) - 11330
Structures containing only coordinate information (no Element assigned).

[Protein Data Bank
Collection](ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/) -
72836 Structures containing both headers and coordinate information.

Versions Tested
---------------

[Biopython 1.49](http://biopython.org/DIST/biopython-1.49.zip) (Nov.
2008)

[Biopython 1.57+](https://github.com/biopython/biopython) (May 2011) |
Element column auto-assignment.

[Biopython PDB
branch](https://github.com/JoaoRodrigues/biopython/tree/pdb_enhancements)
(May 2011 @ Github) | Warnings module replaced
\_handle\_pdb\_exception() && Other minor changes

Benchmarking Script
-------------------

The following script was used to benchmark the parser. The garbage
collection module - <b>gc</b> - was necessary to avoid dead objects
still in memory causing the machine to start swapping.

``` python
#!/usr/bin/env python
""" Script to Benchmark Bio.PDB PDBParser """

import sys, os, warnings

# Parsing Function
def parse_structure(path):
    """ Parses a PDB file """

    s = P.get_structure("test", path)

    return 0


def fancy_output(tps):
    """ Outputs the results in a nicer way """

    print("# Bio.PDB PDBParser Benchmark")
    print()
    print("Structure \tLenght \tTime Spent (ms)")
    for i, s in enumerate(tps):
        print(
            " %s\t (%s) \t%3.3f" % (os.path.basename(pdb_library[i]), pdb_length[i], s)
        )
    print()
    print("Total time spent: %5.3fs" % (sum(tps) / 1000))
    print("Average time per structure: %5.3fms" % (sum(tps) / len(tps)))


if __name__ == "__main__":

    import time, gc
    from Bio.PDB import PDBParser

    P = PDBParser(
        PERMISSIVE=1
    )  # For the pdb_enhancements branch benchmarking, PERMISSIVE was set to 2 (silence warnings).

    library_path = sys.argv[1]

    pdb_library = [os.path.join(library_path, f) for f in os.listdir(library_path)]
    pdb_length = [
        len(set([l[17:26] for l in open(f) if l.startswith("ATOM")]))
        for f in pdb_library
    ]  # Unique counting of residues
    sys.stderr.write(
        "Loaded %s structures (Average Length: %4.3f residues)\n"
        % (len(pdb_length), (sum(pdb_length) / float(len(pdb_length))))
    )

    tps = []
    # Run the Test
    for i, pdb_file in enumerate(pdb_library):
        sys.stderr.write(
            "[%s] %i Structure(s) Parsed \n" % (os.path.basename(pdb_file), i + 1)
        )
        a = time.time()
        parse_structure(pdb_file)
        b = time.time() - a
        tps.append(b * 1000)
        gc.collect()
    # Output Results
    fancy_output(tps)
```

Results
-------

### CATH Dataset

Average Structure Length: 146 residues

#### Biopython 1.49

```
Total Time Spent: 530.686s
Average Time per Structure: 46.84ms/structure
Average Structures per Second: 21.38 structures/s
Failed to parse 0 structures due to errors.

Length                N. Structures   Average Time Spent  ms
< 100                 3660            25.09
100 =< x < 200        5296            44.67
200 =< x < 500        2330            83.40
500 =< x < 1000       43              177.10
>= 1000               1               320.10

TOTAL                 11330           46.84
```

[Link to full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_149.time) |
[Plot of the full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_149.png)

#### Biopython 1.57+

```
Total Time Spent: 686.176s
Average Time per Structure: 60.56ms/structure
Average Structures per Second: 16.51 structures/s
Failed to parse 0 structures due to errors.

Length                N. Structures   Average Time Spent  ms
< 100                 3660            32.55
100 =< x < 200        5296            57.75
200 =< x < 500        2330            107.54
500 =< x < 1000       43              236.62
>= 1000               1               486.602

TOTAL                 11330           60.56
```

[Link to full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_current.time) |
[Plot of the full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_current.png)

#### Biopython PDB Branch

```
Total Time Spent: 695.405s
Average Time per Structure: 61.37ms/structure
Average Structures per Second: 16.29 structures/s
Failed to parse 0 structures due to errors.

Length                N. Structures   Average Time Spent  ms
< 100                 3660            33.21
100 =< x < 200        5296            58.45
200 =< x < 500        2330            108.76
500 =< x < 1000       43              234.37
>= 1000               1               797.583

TOTAL                 11330           61.38
```

[Link to full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_pdb_enhancements.time) |
[Plot of the full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_CATH-biopython_pdb_enhancements.png)

### PDB Dataset

Average Structure Length: 589 residues

Failed to parse 2 structures due to errors:

1. 3NH3 (negative occupancy) - fix: [http://bit.ly/ks6PDN](http://bit.ly/ks6PDN) (still being discussed)
2. 2WMW (invalid ANISOU field) - fix: [http://bit.ly/ld9BWs](http://bit.ly/ld9BWs)

#### Biopython 1.49

```
Total Time Spent: 27801.934s (7.72h)
Average Time per Structure: 381.706ms/structure
Average Structures per Second: 2.62 structures/s

Length                N. Structures   Average Time Spent  (ms)
< 100                 8402            410.270
100 =< x < 200        12226           461.744
200 =< x < 500        26493           182.027
500 =< x < 1000       15878           290.693
>= 1000               9837            942.513

TOTAL                 72836           381.706
```

[Link to full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_PDB-biopython_1.49.time)

#### Biopython 1.57+

```
Total Time Spent: 29516.480s  (~8.20h)
Average Time per Structure: 405.246 ms/structure
Average Structures per Second: 2.47 structures/s

Length                N. Structures   Average Time Spent  (ms)
< 100                 8402            451.819
100 =< x < 200        12226           505.933
200 =< x < 500        26493           190.991
500 =< x < 1000       15878           304.453
>= 1000               9837            980.047

TOTAL                 72836           405.246
```

[Link to full
results](http://nmr.chem.uu.nl/~joao/f/benchmark_PDB-biopython_1.57.time)

#### Biopython PDB Branch (In progress)

In progress
