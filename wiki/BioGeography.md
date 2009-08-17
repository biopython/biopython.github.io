---
title: BioGeography
permalink: wiki/BioGeography
layout: wiki
---

Introduction
------------

BioGeography is a module under development by [Nick
Matzke](User%3AMatzke "wikilink") for a [Google Summer of Code
2009](http://socghop.appspot.com/program/home/google/gsoc2009) project.
It is run through NESCENT's [Phyloinformatics Summer of Code
2009](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009).
See the project proposal at: [Biogeographical Phylogenetics for
BioPython](http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798250).
The mentors are [Stephen Smith](http://blackrim.org/) (primary), [Brad
Chapman](http://bcbio.wordpress.com/), and [David
Kidd](http://evoviz.nescent.org/). The source code is in the
Bio/Geography directory of the [Geography fork of the nmatzke branch on
GitHub](http://github.com/nmatzke/biopython/tree/Geography), and you can
see a timeline and other info about ongoing development of the module
[here](http://biopython.org/wiki/BioGeography). The new module is being
documented on [the BioPython
wiki](http://www.biopython.org/wiki/Main_Page) as
[BioGeography](http://biopython.org/wiki/BioGeography).

**Abstract:** Create a BioPython module that will enable users to
automatically access and parse species locality records from online
biodiversity databases; link these to user-specified phylogenies;
calculate basic alpha- and beta-phylodiversity summary statistics,
produce input files for input into the various inference algorithms
available for inferring historical biogeography; convert output from
these programs into files suitable for mapping, e.g. in Google Earth
(KML files).

Summary of functions
--------------------

Introduction
------------

BioGeography is a module under development by [Nick
Matzke](User%3AMatzke "wikilink") for a [Google Summer of Code
2009](http://socghop.appspot.com/program/home/google/gsoc2009) project.
It is run through NESCENT's [Phyloinformatics Summer of Code
2009](https://www.nescent.org/wg_phyloinformatics/Phyloinformatics_Summer_of_Code_2009).
See the project proposal at: [Biogeographical Phylogenetics for
BioPython](http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798250).
The mentors are [Stephen Smith](http://blackrim.org/) (primary), [Brad
Chapman](http://bcbio.wordpress.com/), and [David
Kidd](http://evoviz.nescent.org/). The source code is in the
Bio/Geography directory of the [Geography fork of the nmatzke branch on
GitHub](http://github.com/nmatzke/biopython/tree/Geography), and you can
see a timeline and other info about ongoing development of the module
[here](http://biopython.org/wiki/BioGeography). The new module is being
documented on [the BioPython
wiki](http://www.biopython.org/wiki/Main_Page) as
[BioGeography](http://biopython.org/wiki/BioGeography).

**Abstract:** Create a BioPython module that will enable users to
automatically access and parse species locality records from online
biodiversity databases; link these to user-specified phylogenies;
calculate basic alpha- and beta-phylodiversity summary statistics,
produce input files for input into the various inference algorithms
available for inferring historical biogeography; convert output from
these programs into files suitable for mapping, e.g. in Google Earth
(KML files).

Summary of functions
--------------------

Tutorial
--------

Bio.Geography is a module for gathering and processing biogeographical
data. The major motivation for the module is to assist analyses of
evolutionary biogeography. A variety of inference algorithms are
available for such analyses, such as
[DIVA](http://www.ebc.uu.se/systzoo/research/diva/manual/dmanual.html)
and [lagrange](http://code.google.com/p/lagrange/). The inputs to such
programs are typically (a) a phylogeny and (b) the areas inhabited by
the species at the tips of the phylogeny. A researcher who has gathered
data on a particular group will likely have direct access to species
location data, but many large-scale analyses may require gathering large
amounts of occurrence data. Automated gathering/processing of occurrence
data has a variety of other applications as well, including species
mapping, niche modeling, error-checking of museum records, and
monitoring range changes.

Occurrence data is derived mainly from museum collections. The major
source of such data is the [Global Biodiversity Information
Facility](http://www.gbif.org/) (GBIF). GBIF serves occurrence data
recorded by hundreds of museums worldwide. GBIF occurrence data can be
[searched manually](http://data.gbif.org/occurrences/), and results
downloaded (see examples on GBIF website) in various formats:
spreadsheet, Google Earth KML, or the XML DarwinCore format.

GBIF can also be accessed via an API. Bio.Geography can process manually
downloaded DarwinCore results, or access GBIF directly.

### Parsing a local (manually downloaded) GBIF DarwinCore XML file

For one-off uses of GBIF, you may find it easiest to just download
occurrence data in spreadsheet format (for analysis) or KML (for
mapping). But for analyses of many groups, or for repeatedly updating an
analysis as new data is added to GBIF, automation is desirable.

A manual search conducted on the GBIF website can return results in the
form of an XML file adhering to the
[DarwinCore](http://en.wikipedia.org/wiki/Darwin_Core) data standard. An
example file can be found in biopython's Tests/Geography directory, with
the name *utric\_search\_v2.xml*. This file contains over 1000
occurrence records for *Utricularia*, a genus of carnivorous plant.

Save the utric\_search\_v2.xml file in your working directory (or
download a similar file from GBIF). Here are suggested steps to parse
the file with Bio.Geography's GbifXml module:

``` python
from Bio.Geography.GbifXml import GbifXmlTree, GbifSearchResults

from Bio.Geography.GenUtils import fix_ASCII_file

xml_fn = 'utric_search_v2.xml'
```

First, in order to display results to screen in python, we need to
convert the file to plain ASCII (GBIF results contain all many of
unusual characters from different languages, and no standardization of
slanted quotes and the like; this can cause crashes when attempting to
print to screen in python).

xml\_fn\_new = fix\_ASCII\_file(xml\_fn)

This creates a new file with the string "\_fixed.xml" added to the
filename.

Next, we will parse the XML file into an ElementTree (a python object
which contains the data from the XML file as a nested series of lists
and dictionaries).

from xml.etree import ElementTree as ET xmltree = ET.parse(xml\_fn\_new)

We can then store the element tree as nn object of Class GbifXmlTree:
gbif\_recs\_xmltree = GbifXmlTree(xmltree)

Then, with the xmltree stored, we parse it into individual records
(stored in individual objects of class GbifObservationRecord), which are
then stored as a group in an object of class GbifSearchResults.

recs = GbifSearchResults(gbif\_recs\_xmltree) recs.latlongs\_to\_obj()

The list of individual observation records can be accessed at
recs.obs\_recs\_list:

print recs.obs\_recs\_list\[0:4\], '...'

To get the data for the first individual record:

rec = recs.obs\_recs\_list\[0\]

dir(rec)

rec.lat will return the latitude, rec.long the longitude, etc. Certain
data attributes are not found in all GBIF records; if they are missing,
the field in question will contain "None".

To print all of the records in a tab-delimited table format:

recs.print\_records()

### Checking how many matching records are hosted by GBIF

Before we go through the trouble of downloading thousands of records, we
may wish to know how many there are in GBIF first. The user must set up
a dictionary containing the fields and search terms as keys and items,
respectively. I.e.,

from GbifXml import GbifXmlTree, GbifSearchResults params = {'format':
'darwin', 'scientificname': 'Utricularia'}

"'format': 'darwin'" specifies that GBIF should return the results in
DarwinCore format. 'scientificname' specifies the genus name to search
on. The full list of search terms can be found on GBIF's [Occurrence
record data service](http://data.gbif.org/tutorial/services), which is
linked from the [Using data from the GBIF
portal](http://data.gbif.org/tutorial/services).

Once you have specified your search parameters, initiate a new
GbifSearchResults object and run get\_numhits to get the number of hits:

recs = GbifSearchResults() numhits = recs.get\_numhits(params)

As of August 2009, 1141 matching records existed in GBIF matching
"Utricularia."

### Downloading an indvidual record

Individual records can be downloaded by key. To download an individual
record:

key = 175067484 recs3 = GbifSearchResults() xmlrec =
recs3.get\_record(key) print xmlrec

### Summary statistics for phylogenetic trees with TreeSum

Biogeographical regions are often characterized by alpha and
beta-diversity statistics: basically, these are indices of the number of
species found within or between regions. Given a phylogeny for organisms
in a region, phylogenetic alpha- and beta-diversity statistics can be
calculated. This has been implemented in a thorough way in the phylocom
package by Webb et al., but for some purposes it is useful to calculate
the statistics directly in python.

Here, we need to start with a Newick tree string:

trstr2 = "(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436,
t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171,
t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662,
t1:0.130208)F:0.0318288)D:0.0273876);"

to2 = Tree(trstr2)

Then, we create a tree summary object:

ts = TreeSum(to2)

The function test\_Tree will run the metrics (MPD = Mean Phylogenetic
Distance, NRI = Net Relatedness Index, MNPD = Mean Nearest Neighbor
Phylogenetic Distance, NTI = Nearest Taxon Index, PD = total
Phylogenetic distance) and output to screen:

ts.test\_Tree()

By subsetting a tree to taxa only existing within a region, statistics
can be calculated by region.
