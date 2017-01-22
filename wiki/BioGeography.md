---
title: Gathering and processing biogeographical data with the BioGeography module.
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
Kidd](https://www.nescent.org/science/awards_summary.php-id=59.html). The source code is in the
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

All classes and functions have been documented with standard docstrings.
Code is available at the most recent GitHub commit here:
<http://github.com/nmatzke/biopython/commits/Geography>

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
[searched manually](http://www.gbif.org/occurrence), and results
downloaded (see examples on GBIF website) in various formats:
spreadsheet, Google Earth KML, or the XML DarwinCore format.

GBIF can also be accessed via an API. Bio.Geography can process manually
downloaded DarwinCore results, or access GBIF directly.

### Background: organization of Bio.Geography

It is useful to understand the overall organization of classes in
Bio.Geography. There are four classes within the GbifXml module:

-   **GbifSearchResults** -- Contains the methods for conducting a GBIF
    search, as well as attributes storing the results in different
    objects, depending on their stage of processing.
    -   Also contains summary statistics on the search (e.g., number of
        records found).
    -   The three objects which store results in different forms are
        **GbifDarwincoreXmlString**, **GbifXmlTree**, and a list of
        individual **GbifObservationRecord** objects.
    -   Method print\_records for printing all contained records
        to screen.
-   **GbifDarwincoreXmlString** -- Contains the raw text returned
    by GBIF. If output to a file, this would be a standard XML file
    adhering to the DarwinCore standard. Inherits from the standard
    python String class.
-   **GbifXmlTree** -- Contains the ElementTree object which results
    from parsing GBIF's XML results. Also a number of methods for
    searching the ElementTree and finding matching elements, finding a
    certain element when it is contained within a certain larger
    element, etc.
-   **GbifObservationRecord** -- Contains the attributes which may be
    found within a certain record, e.g. taxon, genus, species, latitude,
    longitude, etc., as well as functions for classifying a record into
    a certain geographical area, printing a record to screen, etc.
-   **TreeSum** -- Contains functions and attributes for summarizing a
    phylogenetic tree, subsetting it, writing to screen, and calculating
    summary statistics.

### Parsing a local (manually downloaded) GBIF DarwinCore XML file

For one-off uses of GBIF, you may find it easiest to just download
occurrence data in spreadsheet format (for analysis) or KML (for
mapping). But for analyses of many groups, or for repeatedly updating an
analysis as new data is added to GBIF, automation is desirable.

A manual search conducted on the GBIF website can return results in the
form of an XML file adhering to the
[DarwinCore](http://en.wikipedia.org/wiki/Darwin_Core) data standard. An
example file can be found in Biopython's Tests/Geography directory, with
the name *utric\_search\_v2.xml*. This file contains over 1000
occurrence records for *Utricularia*, a genus of carnivorous plant.

Save the utric\_search\_v2.xml file in your working directory (or
download a similar file from GBIF). Here are suggested steps to parse
the file with Bio.Geography's GbifXml module. First, import the
necessary classes and functions, and specify the filename of the input
file.

    from Bio.Geography.GbifXml import GbifXmlTree, GbifSearchResults

    from Bio.Geography.GeneralUtils import fix_ASCII_file

    xml_fn = 'utric_search_v2.xml'

Second, in order to display results to screen in python, we need to
convert the file to plain ASCII (GBIF results contain all many of
unusual characters from different languages, and no standardization of
slanted quotes and the like; this can cause crashes when attempting to
print to screen in python or ipython).

    xml_fn_new = fix_ASCII_file(xml_fn)

This creates a new file with the string "\_fixed.xml" added to the
filename.

Next, we will parse the XML file into an ElementTree (a python object
which contains the data from the XML file as a nested series of lists
and dictionaries).

    from xml.etree import ElementTree as ET
    xmltree = ET.parse(xml_fn_new)

We can then store the element tree as an object of Class GbifXmlTree:

    gbif_recs_xmltree = GbifXmlTree(xmltree)

Then, with the xmltree stored, we parse it into individual records
(stored in individual objects of class GbifObservationRecord), which are
then stored as a group in an object of class GbifSearchResults.

    recs = GbifSearchResults(gbif_recs_xmltree)
    recs.extract_occurrences_from_gbif_xmltree(recs.gbif_recs_xmltree)

The list of individual observation records can be accessed at
recs.obs\_recs\_list. This will display the references to the first five
records:

    recs.obs_recs_list[0:4]

To get the data for the first individual record:

    rec = recs.obs_recs_list[0]

    dir(rec)

rec.lat will return the latitude, rec.long the longitude, etc. Certain
data attributes are not found in all GBIF records; if they are missing,
the field in question will contain "None".

To print all of the records in a tab-delimited table format:

    recs.print_records()

### Checking how many matching records are hosted by GBIF

Before we go through the trouble of downloading thousands of records, we
may wish to know how many there are in GBIF first. The user must set up
a dictionary containing the fields and search terms as keys and items,
respectively. I.e.,

    from GbifXml import GbifXmlTree, GbifSearchResults
    params = {'format': 'darwin', 'scientificname': 'Genlisea*'}

"'format': 'darwin'" specifies that GBIF should return the results in
DarwinCore format.

'scientificname' specifies the genus name to search on. Adding an '\*'
after the name will return anything that begins with "Genlisea".

The full list of search terms can be found on GBIF's [Occurrence record
data service](http://data.gbif.org/tutorial/services), which is linked
from the [Using data from the GBIF
portal](http://data.gbif.org/tutorial/services).

Once you have specified your search parameters, initiate a new
GbifSearchResults object and run get\_numhits to get the number of hits:

    params = {'format': 'darwin', 'scientificname': 'Genlisea*'}
    recs = GbifSearchResults()
    numhits = recs.get_numhits(params)

As of August 2009, 169 matching records existed in GBIF matching
"Genlisea\*"

For constrast, run the same search *without* the asterisk ('\*'):

    params = {'format': 'darwin', 'scientificname': 'Genlisea'}
    numhits = recs.get_numhits(params)

We only get ~10 results -- presumably records of specimens only
identified down to genus and no further.

### Downloading an individual record

Individual records can be downloaded by key. To download an individual
record:

    rec = recs.obs_recs_list[0]
    key = rec.gbifkey
    # (or manually)
    # key = 175067484
    xmlrec = recs.get_record(key)
    print xmlrec

If you want to print the xmlrec ElementTree object, store xmlrec in a
GbifXmlTree object and run print\_xmltree:

    GbifXmlTree(xmlrec).print_xmltree()

### Summary statistics for phylogenetic trees with TreeSum

Biogeographical regions are often characterized by alpha and
beta-diversity statistics: basically, these are indices of the number of
species found within or between regions. Given a phylogeny for organisms
in a region, phylogenetic alpha- and beta-diversity statistics can be
calculated. This has been implemented in a thorough way in the [phylocom
package](http://www.phylodiversity.net/phylocom/) by Webb et al., but
for some purposes it is useful to calculate the statistics directly in
python.

Here, we need to start with a Newick tree string:

    trstr2 = "(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436, t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171, t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662, t1:0.130208)F:0.0318288)D:0.0273876);"

    to2 = Tree(trstr2)

Then, we create a tree summary object:

    ts = TreeSum(to2)

The function test\_Tree will run the metrics (MPD = Mean Phylogenetic
Distance, NRI = Net Relatedness Index, MNPD = Mean Nearest Neighbor
Phylogenetic Distance, NTI = Nearest Taxon Index, PD = total
Phylogenetic distance) and output to screen:

    ts.test_Tree()

By subsetting a tree to taxa only existing within a region, statistics
can be calculated by region.

### Downloading and processing large numbers of records

GBIF only allows a maximum of 1000 observation records to be downloaded
at a time (10,000 for KML records). To get more, we need to download and
process them in stages.

Again we will set up our parameters dictionary, and also an "inc"
variable to specify the number of records to download per server
request.

    params = {'format': 'darwin', 'scientificname': 'Genlisea*'}
    inc = 100
    recs3 = GbifSearchResults()
    gbif_xmltree_list = recs3.get_all_records_by_increment(params, inc)

As with biopython's interactions with NCBI servers, the
GbifSearchResults module keeps track of when the last GBIF request was
made, and requires a 3-second wait before a new request.

Each server request returns an XML string; these are parsed into
GbifXmlTree objects, and a list of the returned GbifXmlTree objects is
returned to gbif\_xmltree\_list. The individual records have also been
parsed:

    recs3.print_records()

### Classifying records into geographical regions

Biogeographical analyses will often require that you determine what
area(s) a taxon lives in. Areas are not always obviously delineated, and
analysts may wish to try several different possible sets of areas and
see how this influences their analysis.

Below, we set up a polygon containing the latitude/longitude coordinates
for the Northern Hemisphere, and then set the "area" attribute for each
matching record to "NorthernHemisphere":

    ul = (-180, 90)
    ur = (180, 90)
    ll = (-180, 0)
    lr = (180, 0)
    poly = [ul, ur, ll, lr]
    polyname = "NorthernHemisphere"

    recs3.print_records()

This process can be repeated for all polygons of interest until all GBIF
records have been classified (except for GBIF records which lacked
lat/long data in the first place, which sometimes happens).

GeogUtils also contains open access libraries for processing
shapefile/dbf files -- these are standard GIS file formats, and various
publicly-accessible shapefiles might serve as sources for polygons.

Warning: the point-in-polygon operation will fail dramatically if your
polygon crosses the International Dateline. The best solution in this
case is to split any polygons crossing the dateline into two polygons,
one on each side of the line.

### General notes

GBIF search results often contain non-ASCII characters (e.g.
international placenames) and other confusing items, e.g., web links in
angle brackets, which can be misinterpreted as unmatched XML tags if a
GBIF search result is read to ASCII and then an attempt is made to parse
it.

In general, the Geography module will handle things fine if the results
are being processed in the background; but to print results to screen, a
series of functions from GeneralUtils are used to convert a string to
plain ASCII. This avoids crashes e.g. when printing data to screen.
Therefore, these printed-to-screen results may slightly alter the
content of the original search results.
