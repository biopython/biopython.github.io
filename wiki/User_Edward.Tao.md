---
title: User:Edward.Tao
permalink: wiki/User%3AEdward.Tao
layout: wiki
---

Extract intergenic regions from genome file (genbank).

modified from script of 2009 Iddo Friedberg & Ian MC Fleming
(http://biopython.org/wiki/Intergenic\_regions)

The function:

``` python
def get_interregions(infiles_gbff, intergene_length=1):
   import sys
   from Bio import SeqIO
   from Bio.SeqRecord import SeqRecord
   import os
   databuff = []
   for line in file(infiles_gbff):
       if line.strip() != '//':
           databuff.append(line)
       else:
           genbank_path = open('temp.gb', 'w')
           for dataline in databuff:
               genbank_path.write(dataline)
           genbank_path.write('//')
           genbank_path.close()
           seq_record = SeqIO.parse(open('temp.gb'), "genbank").next()
           cds_list_plus = []
           cds_list_minus = []
           intergenic_records = []
           for feature in seq_record.features:
               print(feature.location._start.position)
               if feature.type == 'CDS':
                   mylocal_tag = feature.qualifiers['locus_tag'][0]
                   myprotein_id = feature.qualifiers['protein_id'][0]
                   mystart = feature.location._start.position
                   myend = feature.location._end.position
                   if feature.strand == -1:
                       cds_list_minus.append((mylocal_tag, mystart, myend, myprotein_id, -1))
                   elif feature.strand == 1:
                       cds_list_plus.append((mylocal_tag, mystart, myend, myprotein_id, 1))
                   else:
                       sys.stderr.write("No strand indicated %d-%d. Assuming +\n" %
                                        (mystart, myend))
                       cds_list_plus.append( ( mystart, myend, 1))

           for i, pospair in enumerate( cds_list_plus[1:]):
               last_end = cds_list_plus[i][2]
               this_start = pospair[1]
               strand = pospair[4]
               locus_tag = pospair[0]
               protein_id = pospair[3]
               if this_start - last_end >= intergene_length:
                   intergene_seq = seq_record.seq[last_end:this_start]
                   strand_string = "+"
                   intergenic_records.append(SeqRecord(intergene_seq, id="%s-ign-%d" % (locus_tag, i), description="%s %d-%d %s" % (protein_id, last_end + 1, this_start, strand_string)))
           for i, pospair in enumerate(cds_list_minus[1:]):
               last_end = cds_list_minus[i][2]
               this_start = pospair[1]
               strand = pospair[4]
               locus_tag = pospair[0]
               protein_id = pospair[3]
               if this_start - last_end >= intergene_length:
                   intergene_seq = seq_record.seq[last_end:this_start]
                   strand_string = "-"
                   intergenic_records.append(SeqRecord(intergene_seq, id="%s %d" % (locus_tag, i), description="%s %d-%d %s" % (protein_id, last_end + 1, this_start, strand_string)))
           outpath = os.path.splitext(os.path.basename(infiles_gbff))[0] + "_ign.fasta"
           SeqIO.write( intergenic_records, open( outpath, "aw" ), "fasta")

           databuff = []
   return
```
