---
title: BioSQL_Windows
permalink: wiki/BioSQL_Windows
layout: wiki
---

These are Windows specific notes on installing
[BioSQL](BioSQL "wikilink")

### MySQL

Download and install MySQL 5.1 from
<http://dev.mysql.com/downloads/mysql/5.1.html#win32>

Download and install MySQLdb (the python interface) from
<http://sourceforge.net/projects/mysql-python>

Create a new database, and load the schema etc as described on
[BioSQL](BioSQL "wikilink")

``` bash
mysqladmin -u root create bioseqdb
mysql -p -u root bioseqdb < biosqldb-mysql.sql
```