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

Work Plan
---------

Note: all major functions are being placed in the file geogUtils.py for
the moment. Also, the immediate goal is to just get everything basically
working, so details of where to put various functions, what to call
them, etc. are being left for later.

Code usage: For a few things, an entire necessary function already
exists (e.g. for reading a shapefile), and re-inventing the wheel seems
pointless. In most cases the material used appears to be open source
(e.g. previous Google Summer of Code). For a few short code snippets
found online in various places I am less sure. In all cases I am noting
the source and when finalizing this project I will go back and determine
if the stuff is considered copyright, and if so email the authors for
permission to use.

### May, week 1: Functions to read locality data and place points in geographic regions (Tasks 1-2)

#### readshpfile

Parses polygon, point, and multipoint shapefiles into python objects
(storing latitude/longitude coordinates and feature names, e.g. the
region name associated with each polygon)

#### extract\_latlong

Parse a manually downloaded GBIF record, extracting latitude/longitude
and taxon names

#### shapefile\_points\_in\_poly, tablefile\_points\_in\_poly

Input geographic points, determine which region (polygon) each range
falls in (via point-in-polygon algorithm); also output points that are
unclassified, e.g. some GBIF locations were mis-typed in the source
database, so a record will fall in the middle of the ocean.

#### Code

-   [Code fulfilling these tasks is uploaded
    here](http://github.com/nmatzke/biopython/commit/4d963a65ce48b9d50327f191dedcc76abbb149be),
    along with an example script and data files to run.

### June, week 1: Functions to search GBIF and download occurrence records

Note: creating functions for all possible interactions with GBIF is not
possible in the time available, I will just focus on searching and
downloading basic record occurrence record data.

#### access\_gbif

utility function invoked by other functions, user inputs parameters and
the GBIF response in XML/DarwinCore format is returned. The relevant
GBIF web service, and the search commands etc., are here:
<http://data.gbif.org/ws/rest/occurrence>

#### get\_hits

Get the actual hits that are be returned by a given search, returns
filename were they are saved

#### get\_xml\_hits

Like get\_hits, but returns a parsed XML tree

#### fix\_ASCII

files downloaded from GBIF contain HTML character entities & unicode
characters (e.g. umlauts mostly) which mess up printing results to
prompt in Python, this fixes that

#### paramsdict\_to\_string

converts user's search parameters (in python dictionary format; see here
for params <http://data.gbif.org/ws/rest/occurrence> ) to a string for
submission via access\_gbif

#### xmlstring\_to\_xmltree(xmlstring)

Take the text string returned by GBIF and parse to an XML tree using
ElementTree. Requires the intermediate step of saving to a temporary
file (required to make ElementTree.parse work, apparently).

#### element\_items\_to\_dictionary

If the XML tree element has items encoded in the tag, e.g. key/value or
whatever, this function puts them in a python dictionary and returns
them.

#### extract\_numhits

Search an element of a parsed XML string and find the number of hits, if
it exists. Recursively searches, if there are subelements.

#### print\_xmltree

Prints all the elements & subelements of the xmltree to screen (may
require fix\_ASCII to input file to succeed)

#### Deleted (turns out this was unnecessary): gettaxonconceptkey

user inputs a taxon name and gets the GBIF key back (useful for
searching GBIF records and finding e.g. synonyms and daughter taxa). The
GBIF taxon concepts are accessed via the taxon web service:
<http://data.gbif.org/ws/rest/taxon>

#### Code

-   [Code fulfilling these tasks is uploaded
    here](http://github.com/nmatzke/biopython/commits/Geography), along
    with an example script and data files to run.

### June, week 2: Functions to get GBIF records

Added functions download & parse large numbers of records, get
TaxonOccurrence gbifKeys, and search with those keys.

#### get\_record

Retrieves a single specified record in DarwinCore XML format, and
returns an xmltree for it.

#### extract\_occurrence\_elements

Returns a list of the elements, picking elements by TaxonOccurrence;
this should return a list of elements equal to the number of hits.

#### extract\_taxonconceptkeys\_tolist

Searches an element in an XML tree for TaxonOccurrence gbifKeys, and the
complete name. Searches recursively, if there are subelements. Returns
list.

#### extract\_taxonconceptkeys\_tofile

Searches an element in an XML tree for TaxonOccurrence gbifKeys, and the
complete name. Searches recursively, if there are subelements. Returns
file at outfh.

#### get\_all\_records\_by\_increment

Download all of the records in stages, store in list of elements.
Increments of e.g. 100 to not overload server. Currently stores results
in a list of tempfiles which is returned (could return a list of handles
I guess).

#### Code

Updated functions have been pushed to Github
[here](http://github.com/nmatzke/biopython/commit/5df9025ea5cd3458915db982c69422345e1da8d7)

### June, week 3: Functions to read user-specified Newick files (with ages and internal node labels) and generate basic summary information.

(note: I have scripts doing all of these functions already, so the work
is integrating them into a Biopython module, testing them, etc.)

#### read\_ultrametric\_Newick(newickstr)

Read a Newick file into a tree object (a series of node objects links to
parent and daughter nodes), also reading node ages and node labels if
any.

#### list\_leaves(phylo\_obj)

Print out all of the leaves in above a node object

#### treelength(node)

Gets the total branchlength above a given node by recursively adding
through tree.

#### phylodistance(node1, node2)

Get the phylogenetic distance (branch length) between two nodes.

#### get\_distance\_matrix(phylo\_obj)

Get a matrix of all of the pairwise distances between the tips of a
tree.

#### get\_mrca\_array(phylo\_obj)

Get a square list of lists (array) listing the mrca of each pair of
leaves (half-diagonal matrix)

#### subset\_tree(phylo\_obj, list\_to\_keep)

Given a list of tips and a tree, remove all other tips and resulting
redundant nodes to produce a new smaller tree.

#### prune\_single\_desc\_nodes(node)

Follow a tree from the bottom up, pruning any nodes with only one
descendant

#### find\_new\_root(node)

Search up tree from root and make new root at first divergence

#### make\_None\_list\_array(xdim, ydim)

Make a list of lists ("array") with the specified dimensions

#### get\_PD\_to\_mrca(node, mrca, PD)

Add up the phylogenetic distance from a node to the specified ancestor
(mrca). Find mrca with find\_1st\_match.

#### find\_1st\_match(list1, list2)

Find the first match in two ordered lists.

#### get\_ancestors\_list(node, anc\_list)

Get the list of ancestors of a given node

#### addup\_PD(node, PD)

Adds the branchlength of the current node to the total PD measure.

#### print\_tree\_outline\_format(phylo\_obj)

Prints the tree out in "outline" format (daughter clades are indented,
etc.)

#### print\_Node(node, rank)

Prints the node in question, and recursively all daughter nodes,
maintaining rank as it goes.

#### lagrange\_disclaimer()

Just prints lagrange citation etc. in code using lagrange libraries.

#### Code

-   [Code fulfilling these tasks is uploaded
    here](http://github.com/nmatzke/biopython/commits/Geography), along
    with an example script and data files to run.

### June, week 4: Functions to summarize taxon diversity in regions, given a phylogeny and a list of taxa and the regions they are in.

(note: I have scripts doing all of these functions already, so the work
is integrating them into a Biopython module, testing them, etc.)

Priority for this week:

Following up on suggestions to make the code more standard, with the
priority of figuring out how I can revise the current BioPython
phylogeny class, to resemble the better version in lagrange, so that
there is a generic flexible phylogeny/newick parser that can be used
generally as well as by my BioGeography package specifically.

#### alphadiversity

alpha diversity of a region (number of taxa in the region)

#### betadiversity

beta diversity (Sorenson’s index) between two regions

#### alphaphylodistance

total branchlength of a phylogeny of taxa within a region

#### phylosor

phylogenetic Sorenson’s index between two regions

#### meanphylodistance

average distance between all tips on a region’s phylogeny

#### meanminphylodistance

average distance to nearest neighbor for tips on a region’s phylogeny

#### netrelatednessindex

standardized index of mean phylodistance

#### nearesttaxonindex

standardized index of mean minimum phylodistance

### July, week 1: lagrange input/output handling (Task 6)

(note: lagrange requires a number of input files, e.g. hypothesized
histories of connectivity; the only inputs suitable for automation in
this project are the species ranges and phylogeny

#### make\_lagrange\_species\_range\_inputs

convert list of taxa/ranges to input format:
<http://www.reelab.net/lagrange/configurator/index>

#### check\_input\_lagrange\_tree

checks if input phylogeny meets the requirements for lagrange, i.e. has
ultrametric branchlengths, tips end at time 0, tip names are in the
species/ranges input file

#### parse\_lagrange\_output

take the output file from lagrange and get ages and estimated regions
for each node

### July, weeks 2-3: Devise algorithm for representing estimated node histories (location of nodes in categorical regions) as latitude/longitude points, necessary for input into geographic display files.

-   Regarding where to put reconstructed nodes, or tips that where the
    only location information is region. Within regions, dealing with
    linking already geo-located tips, spatial averaging can be used as
    currently happens with GeoPhyloBuilder. If there is only one node in
    a region the centroid or something similar could be used (i.e. the
    "root" of the polygon skeleton would deal even with weird
    concave polygons).
-   If there are multiple ancestral nodes or region-only tips in a
    region, they need to be spread out inside the polygon, or lines will
    just be drawn on top of each other. This can be done by putting the
    most ancient node at the root of the polygon skeleton/medial axis,
    and then spreading out the daughter nodes along the skeleton/medial
    axis of the polygon.

#### get\_polygon\_skeleton

this is a standard operation:
<http://en.wikipedia.org/wiki/Straight_skeleton>

#### assign\_node\_locations\_in\_region

within a region’s polygon, given a list of nodes, their relationship,
and ages, spread the nodes out along the middle 50% of the longest axis
of the polygon skeleton, with the oldest node in the middle

#### assign\_node\_locations\_between\_regions

connect the nodes that are linked to branches that cross between regions
(for this initial project, just the great circle lines)

### July, week 4 and August, week 1: Write functions for converting the output from the above into graphical display formats, e.g. shapefiles for ArcGIS, KML files for Google Earth.

#### write\_history\_to\_shapefile

write the biogeographic history to a shapefile

#### write\_history\_to\_KML

write the biogeographic history to a KML file for input into Google
Earth

### August, week 2: Beta testing

Make the series of functions available, along with suggested input
files; have others run on various platforms, with various levels of
expertise (e.g. Evolutionary Biogeography Discussion Group at U.C.
Berkeley). Also get final feedback from mentors and advisors.

### August, week 3: Wrapup

Assemble documentation, FAQ, project results writeup for
Phyloinformatics Summer of Code.
