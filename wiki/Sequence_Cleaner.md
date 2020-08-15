---
title: Sequence Cleaner
permalink: wiki/Sequence_Cleaner
layout: wiki
tags:
 - Cookbook
redirect_from:
 - wiki/Create_an_Article_to_this_category
---

Description
-----------

I want to share my script using Biopython to clean sequences up. You
should know that analyzing poor data takes CPU time and interpreting the
results from poor data takes people time, so it's always important to
make a preprocessing.

Let me call my script as “Sequence\_cleaner” and the big idea is to
remove duplicate sequences, remove too short sequences (the user
defines the minimum length) and remove sequences which have too many
unknown nucleotides (N) (the user defines the % of N it allows ) and in
the end the user can choose if he/she wants to have a file as output or
print the result.

Script
------

``` python
import sys
from Bio import SeqIO


def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    # Create our hash table to add the sequences
    sequences = {}

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (
            len(sequence) >= min_length
            and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
        ):
            # If the sequence passed in the test "is it clean?" and it isn't in the
            # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
            # If it is already in the hash table, we're just gonna concatenate the ID
            # of the current sequence to another one that is already in the hash table
            else:
                sequences[sequence] += "_" + seq_record.id

    # Write the clean sequences

    # Create a file in the same directory where you ran this script
    with open("clear_" + fasta_file, "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + fasta_file)


userParameters = sys.argv[1:]

try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], float(userParameters[1]))
    elif len(userParameters) == 3:
        sequence_cleaner(
            userParameters[0], float(userParameters[1]), float(userParameters[2])
        )
    else:
        print("There is a problem!")
except:
    print("There is a problem!")
```

Using the command line, you should run:

``` bash
python sequence_cleaner.py input[1] min_length[2] min[3]
```

There are 3 basic parameters:

-   \[1\]: your fasta file
-   \[2\]: the user defines the minimum length. Default value 0, it means you
    don't have to care about the minimum length
-   \[3\]: the user defines the % of N is allowed. Default value 100, all
    sequences with 'N' will be in your ouput, set value to 0 if you want no
    sequences with "N" in your output

For example:

``` bash
python sequence_cleaner.py Aip_coral.fasta 10 10
```

FYI: if you don't care about the 2nd and the 3rd parameters you are just
gonna remove the duplicate sequences

Questions, Suggestions or Improvement
-------------------------------------

Send an email to <genivaldo.gueiros@gmail.com>
