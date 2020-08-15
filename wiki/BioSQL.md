---
title: Managing local biological databases with the BioSQL module.
permalink: wiki/BioSQL
layout: wiki
tags:
 - Wiki Documentation
---

[BioSQL](http://www.biosql.org/wiki/Main_Page) is a joint effort between
the [OBF](http://open-bio.org/) projects (BioPerl, BioJava etc) to
support a shared database schema for storing sequence data. In theory,
you could load a GenBank file into the database with BioPerl, then using
Biopython extract this from the database as a record object with
features - and get more or less the same thing as if you had loaded
the GenBank file directly as a [SeqRecord](SeqRecord "wikilink") using
[SeqIO](SeqIO "wikilink").

We have some existing documentation
([HTML](http://biopython.org/DIST/docs/biosql/python_biosql_basic.html),
[PDF](http://biopython.org/DIST/docs/biosql/python_biosql_basic.pdf))
for the Biopython interfaces to BioSQL, covering installing Python
database adaptors and basic usage of BioSQL. This is a little old, and I
am hoping to use this wiki page to update the above documentation in
future.

The following text applies to Biopython 1.45 or later (and won't work
with Biopython 1.44).

Installation
============

This is fairly complicated - partly because there are so many options.
For example, you can use a range of different SQL database packages
(we'll focus on MySQL), you can have the database on your own computer
(the assumption here) or on a separate server, and of course there are
usernames and passwords associated with database. And finally the
details will also vary depending on your operating system.

This text is based in part on the [BioSQL schema INSTALL
instructions](https://github.com/biosql/biosql/blob/master/INSTALL)
which also covers alternatives to MySQL.

Installing Required Software
----------------------------

You will need to install some database software plus the associated
python library so that Biopython can "talk" to the database. In this
example we'll talk about the most common choice, MySQL. How you do this
will also depend on your operating system, for example on a Debian or
Ubuntu Linux machine try this:

``` bash
sudo apt-get install mysql-common mysql-server python-mysqldb
```

A password is required for logon. Please refer to [MySQL documentation](https://dev.mysql.com/doc/refman/8.0/en/linux-installation-debian.html).

It will also be important to have perl (to run some of the setup
scripts). Again, on a Debian or Ubuntu Linux machine try this:

``` bash
sudo apt-get install perl
```

You may find perl is already installed.

For Windows users, see [BioSQL on Windows](BioSQL_Windows "wikilink").

Downloading the BioSQL Schema & Scripts
---------------------------------------

Once the software is installed, your next task is to setup a database
and import the BioSQL schema (i.e. setup the relevant tables within the
database). See [BioSQL downloads](http://www.biosql.org/wiki/Downloads)
-- you'll need to unzip the archive.

Alternatively to get the very latest BioSQL, check out their git
repository. Or, navigate to the relevant schema file for your database
and download just that, e.g.
[biosqldb-mysql.sql](https://raw.github.com/biosql/biosql/master/sql/biosqldb-mysql.sql)
for MySQL. You will also want the NCBI Taxonomy loading perl script,
[load\_ncbi\_taxonomy.pl](https://raw.github.com/biosql/biosql/master/scripts/load_ncbi_taxonomy.pl).

Creating the empty database
---------------------------

### MySQL

The following command line should create a new database on your own
computer called *bioseqdb*, belonging to the *root* user account:

``` bash
mysqladmin -u root -p create biosqldb
```

We can then tell MySQL to load the BioSQL scheme we downloaded above.
Change to the scripts subdirectory from the unzipped BioSQL download,
then:

``` bash
mysql -u root -p bioseqdb < biosqldb-mysql.sql
```

You can have a quick play using the mysql command line tool, for
example:

``` bash
mysql --user=root -p bioseqdb -e "show tables"
```

giving:

```
+----------------------------+
| Tables_in_bioseqdb         |
+----------------------------+
| biodatabase                |
| bioentry                   |
| bioentry_dbxref            |
| bioentry_path              |
| bioentry_qualifier_value   |
| bioentry_reference         |
| bioentry_relationship      |
| biosequence                |
| comment                    |
| dbxref                     |
| dbxref_qualifier_value     |
| location                   |
| location_qualifier_value   |
| ontology                   |
| reference                  |
| seqfeature                 |
| seqfeature_dbxref          |
| seqfeature_path            |
| seqfeature_qualifier_value |
| seqfeature_relationship    |
| taxon                      |
| taxon_name                 |
| term                       |
| term_dbxref                |
| term_path                  |
| term_relationship          |
| term_relationship_term     |
| term_synonym               |
+----------------------------+
```

Or, to look inside a table:

``` bash
mysql --user=root -p bioseqdb -e "select * from bioentry;"
```

This should return no rows as the table is empty.

### PostgreSQL

**IMPORTANT NOTE FOR POSTRESQL USERS**: Before loading the
biosqldb-pg.sql schema into Postgres you must delete the two RULES named
**rule\_bioentry\_i1** and **rule\_bioentry\_i2**; lines 771-791 in
biosqldb-pg.sql BioSQL version 1.0.1

First you need to set up user permissions, if you are not sure how to do
this, try:

``` bash
su - postgres
createuser <your user name>
```

Then, assuming you are logged-in as <your user name> and Postgres is
running on the local machine, you should be able to do the following:

``` bash
createdb biosqldb
psql biosqldb < biosqldb-pg.sql
```

Run *psql* and type enter *\\d <ENTER>* to see all the entities created.

NCBI Taxonomy
-------------

The BioSQL package includes a perl script under
scripts/load\_ncbi\_taxonomy.pl to download and update the taxonomy
tables. The script should be able to download the files it needs from
the [NCBI taxonomy FTP site](ftp://ftp.ncbi.nih.gov/pub/taxonomy/)
automatically.

Prior to Biopython 1.49, if you wanted to work with the NCBI taxonomy
database it was good idea to pre-load the NCBI taxonomy before you start
trying to load sequences into the database. This isn't so important with
Biopython 1.49 onwards, where you can instead opt to have the
information needed downloaded as needed from Entrez.

To update the NCBI taxonomy, change to the scripts subdirectory from the
unzipped BioSQL download, then:

``` bash
./load_ncbi_taxonomy.pl --dbname bioseqdb --driver mysql --dbuser root --download true
```

For PostgreSQL you need to have the perl DBD-Pg module installed - Using
CPAN: "perl -MCPAN -e 'install DBI'; perl -MCPAN -e 'install DBD::Pg'".
Substitute *Pg* for *mysql* in the above command.

There is about 10MB to fetch, so it can take a little while (and doesn't
give any feedback while this happens). If you are worried, open a file
browser window and check to see it is downloading a file called
taxdump.tar.gz to the taxdata subdirectory.

You should see this output at the command prompt - be warned that some
of these steps do take a while (especially *rebuilding nested set
left/right values*):

```
Loading NCBI taxon database in taxdata:
        ... retrieving all taxon nodes in the database
        ... reading in taxon nodes from nodes.dmp
        ... insert / update / delete taxon nodes
        ... (committing nodes)
        ... rebuilding nested set left/right values
        ... reading in taxon names from names.dmp
        ... deleting old taxon names
        ... inserting new taxon names
        ... cleaning up
Done.
```

This might be a good point for a tea break - I didn't time this but it
was over ten minutes.

One the initial tables have been populated, re-running the script is
much faster. You can run this script again to update the taxonomy
tables, which the NCBI do add to regularly. You may want to setup a
scheduled job to do this automatically (say once a fortnight).

P.S. It is a particularly good idea to do update the taxonomy if you
will be working with the left/right values in the taxon table (see also
[BioSQL enhancement request GitHub issue 14
(Redmine 2493)](https://github.com/biosql/biosql/issues/14)).
Biopython 1.67 onwards will do this when writing to the taxon tables.

Note Biopython ignores these optional fields when loading or retrieving
sequences - instead using just the parent link. See
<http://www.oreillynet.com/pub/a/network/2002/11/27/bioconf.html> for
more about how this alternative tree representation works.

Running the unit tests
----------------------

Because there are so many ways you could have setup your BioSQL
database, you have to tell the unit test a few bits of information by
editing the file Tests/setup\_BioSQL.py and filling in the following
fields:

``` python
DBDRIVER = "MySQLdb"
DBTYPE = "mysql"
```

and a little lower down,

``` python
DBHOST = "localhost"
DBUSER = "root"
DBPASSWD = "your-password"
TESTDB = "biosql_test"
```

Change these to match your setup. You can then run the BioSQL unit tests
as normal, e.g.

``` bash
python run_tests.py test_BioSQL test_BioSQL_SeqIO
```

For PostgreSQL, use:

``` python
DBDRIVER = "psycopg2"
DBTYPE = "pg"
```

Creating a (sub) database
=========================

BioSQL lets us define named "sub" databases or "namespaces" within the
single SQL database (which we called *bioseqdb* earlier). For this
example, lets create a one for some orchid sequences:

``` python
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(
    driver="MySQLdb",
    user="root",
    passwd="your-password",
    host="localhost",
    db="bioseqdb",
)
db = server.new_database("orchids", description="Just for testing")
server.commit()  # On Biopython 1.49 or older, server.adaptor.commit()
```

(If you are using PostgreSQL rather than MySQL, just change the driver
argument to "psycopg2" instead. The same applies to the other examples
in this document)

The *commit* call tells the database to save the changes so far (commit
the SQL transaction). It is up to you to decide when to commit the SQL
transaction(s), and/or rollback changes, rather than having Biopython
try and decide for you and risk getting it wrong. See *Explicit is
better than implicit* ([The Zen of
Python](http://www.python.org/dev/peps/pep-0020/)).

There should now be a single row in the *biodatabase* table for our new
orchid namespace. You can check this at the command line:

MySQL:

``` bash
mysql --user=root -p bioseqdb -e "select * from biodatabase;"
```

PostgreSQL:

``` bash
psql -c "SELECT * FROM biodatabase;" bioseqdb
```

Which should give something like this (assuming you haven't done any
other testing yet):

```
+----------------+---------+-----------+------------------+
| biodatabase_id | name    | authority | description      |
+----------------+---------+-----------+------------------+
|              1 | orchids | NULL      | Just for testing |
+----------------+---------+-----------+------------------+
```

Now that we have setup an *orchids* namespace within our *biosqldb*
MySQL database, lets add some sequences to it.

Loading Sequences into a database
=================================

When loading sequences into a BioSQL database with Biopython we have to
provide annotated [SeqRecord](SeqRecord "wikilink") objects. This gives
us another excuse to use the [SeqIO](SeqIO "wikilink") module! A quick
recap on reading in sequences as SeqRecords, based on one of the orchid
examples in the Biopython Tutorial:

``` python
from Bio import Entrez
from Bio import SeqIO

handle = Entrez.efetch(
    db="nuccore", id="6273291,6273290,6273289", rettype="gb", retmode="text"
)
for seq_record in SeqIO.parse(handle, "genbank"):
    print seq_record.id, seq_record.description[:50] + "..."
    print "Sequence length %i," % len(seq_record.seq),
    print "from: %s" % seq_record.annotations["source"]
handle.close()
```

The expected output is below, note we have three records with a total of
nine features:

```
AF191665.1 Opuntia marenae rpl16 gene; chloroplast gene for c...
Sequence length 902, 3 features, from: chloroplast Opuntia marenae
AF191664.1 Opuntia clavata rpl16 gene; chloroplast gene for c...
Sequence length 899, 3 features, from: chloroplast Grusonia clavata
AF191663.1 Opuntia bradtiana rpl16 gene; chloroplast gene for...
Sequence length 899, 3 features, from: chloroplast Opuntia bradtianaa
```

Now, instead of printing things on screen, let's add these three records
to a new (empty) *orchid* database:

``` python
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(
    driver="MySQLdb",
    user="root",
    passwd="your-password",
    host="localhost",
    db="bioseqdb",
)
db = server["orchids"]
handle = Entrez.efetch(
    db="nuccore", id="6273291,6273290,6273289", rettype="gb", retmode="text"
)
count = db.load(SeqIO.parse(handle, "genbank"))
print "Loaded %i records" % count
server.commit()  # On Biopython 1.49 or older, server.adaptor.commit()
```

Again, you must explicitly call *commit* to record the SQL transaction
which is otherwise left pending.

The *db.load()* function should have returned the number of records
loaded (three in this example), and again have a look in the database
and you should see new rows in several tables.

The *bioentry* and *biosequence* tables should have three new rows:

``` bash
mysql --user=root bioseqdb -e "select * from bioentry;"
mysql --user=root bioseqdb -e "select * from biosequence;"
```

The should also be nine new features:

``` bash
mysql --user=root -p bioseqdb -e "select * from seqfeature;"
```

Next, we'll try and load these three records back from the database.

Extracting Sequences from the database
======================================

This continues from the previous example, where we loaded three records
into an *orchids* database (namespace):

``` python
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(
    driver="MySQLdb",
    user="root",
    passwd="your-password",
    host="localhost",
    db="bioseqdb",
)
db = server["orchids"]
for identifier in ["6273291", "6273290", "6273289"]:
    seq_record = db.lookup(gi=identifier)
    print seq_record.id, seq_record.description[:50] + "..."
    print "Sequence length %i," % len(seq_record.seq)
```

Giving:

```
AF191665.1 Opuntia marenae rpl16 gene; chloroplast gene for c...
Sequence length 902
AF191664.1 Opuntia clavata rpl16 gene; chloroplast gene for c...
Sequence length 899
AF191663.1 Opuntia bradtiana rpl16 gene; chloroplast gene for...
Sequence length 899
```

The objects you get back from BioSQL act like a
[SeqRecord](SeqRecord "wikilink") object with a [Seq](Seq "wikilink")
object as the sequence, but they are not exactly the same. You actually
get their BioSQL database equivalent, a **DBSeqRecord** object with a
**DBSeq** object for the sequence. These will only load the sequence and
annotation from the database on demand.

There are other ways to pull out records - the 'db' object here acts
somewhat like a dictionary (including supporting deleting entries from
their key). The python-dictionary-keys are actually the
database-primary-keys used inside the database for the bioentry table.
e.g.

``` python
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(
    driver="MySQLdb",
    user="root",
    passwd="your-password",
    host="localhost",
    db="bioseqdb",
)
db = server["orchids"]
print "This database contains %i records" % len(db)
for key, record in db.iteritems():
    print "Key %r maps to a sequence record with id %s" % (key, record.id)
```

Deleting a (sub) database
=========================

As mentioned above, BioSQL lets us define named "sub" databases (aka
namespaces) within the single SQL database (which we called *bioseqdb*).
In the previous example, we created a sub-database for some orchid
sequences. The following code will delete the *orchid* database (and all
the records in it):

``` python
from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(
    driver="MySQLdb",
    user="root",
    passwd="your-password",
    host="localhost",
    db="bioseqdb",
)
server.remove_database("orchids")
server.commit()  # On Biopython 1.49 or older, server.adaptor.commit()
```

Again, you must explicitly finialise the SQL transaction with a commit
to apply this change.

There should now be one less row in the *biodatabase* table, check this
at the command line:

``` bash
mysql --user=root -p bioseqdb -e "select * from biodatabase;"
```

You can also check that the three orchid sequences have gone from the
other tables.

How is the data stored
======================

If you need or want to go direct to the data, bypassing the Biopython
methods to retrieve records, then it would help to know more about how
the underlying tables are used. For this, we refer you to the BioSQL
documentation, starting with their [Schema
Overview](http://www.biosql.org/wiki/Schema_Overview) and the page on
[Annotation Mapping](http://www.biosql.org/wiki/Annotation_Mapping).

MySQL Tips and Tricks
=====================

If you are getting timeout errors, check to see if your SQL server has
any orphaned threads.

``` bash
mysql --user=root -p bioseqdb -e "SHOW INNODB STATUS\G" | grep "thread id"
```

And if there are, assuming you are the only person using this database,
you might try killing them off using the thread id given by the above
command:

``` bash
mysql --user=root -p bioseqdb -e "KILL 123;"
```

Use at your own risk!

Jython
======

Jython is supported for both PostgreSQL and MySQL from version 1.62.
Usage is almost fully transparent: You access the databases using the
exact same configuration parameters as specified above. The only issue
to take into consideration is that you will need the correct JBDC driver
in your CLASSPATH. For both databases you should use the official JDBC
drivers (available [here for PosgreSQL](http://jdbc.postgresql.org/) and
[here for MySQL](http://dev.mysql.com/downloads/connector/j/))
