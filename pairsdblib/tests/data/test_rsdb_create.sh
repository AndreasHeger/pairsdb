#!/bin/bash

DATABASE=test

mysql -B -N -e "SELECT query_nid, sbjct_nid, evalue, query_from, query_to, query_ali, sbjct_from, sbjct_to, sbjct_ali, score, percent_identity
	FROM pairsdb_90x90 ORDER by query_nid, sbjct_nid, evalue" \
	${DATABASE} > pairsdb_90x90.links

mysql -B -N -e "SELECT nrdb.nid, length FROM nrdb, nrdb90 where nrdb.nid = nrdb90.nid" ${DATABASE} > nrdb.sizes

## try the database way
../../src/rsdb_create -v 2 --user test --database test --prefix-groupies=db_groupies --prefix-rsdb=db_rsdb --socket /var/lib/mysql/mysql.sock 30 40 50 60 70 80 >& rsdb_create_fromdb.log

## try the filesystem way
../../src/rsdb_create -v 2 -n pairsdb_90x90.links --prefix-groupies=fi_groupies --prefix-rsdb=fi_rsdb -q nrdb.sizes 30 40 50 60 70 80 >& rsdb_create_fromfile.log

## test with pre-built index
../../src/rsdb_index -n pairsdb_90x90.links -i rsdb.index >& rsdb_index.log
../../src/rsdb_create -v 2 -n pairsdb_90x90.links -i rsdb.index --prefix-groupies=pb_groupies --prefix-rsdb=pb_rsdb -q nrdb.sizes 30 40 50 60 70 80 >& rsdb_create_prebuilt.log
