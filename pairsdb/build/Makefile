PARAM_DIR_INPUT=../input/
PARAM_DIR_SRC=../src/

#############################################################
#############################################################
#############################################################
## primary targets
#############################################################
all: nrdb.table nrdb90.table nrdb90_masks.table \
	nrdb.index nrdb.sary \
	pairsdb_100x90.table \
	pairsdb_90x90.table \
	pairsdb_90x40.table \
	pairsdb_100x40.table \
	scop.log \
	cath.log \
	pdb.log \
	interpro.log

clean:
	rm -rf *.table* *.log nrdb* sary* indexed*

#############################################################
#############################################################
#############################################################
## secondary targets
#############################################################
EXTRA_INPUT_FILES=$(wildcard *.fasta)
EXTRA_INPUT=$(EXTRA_INPUT_FILES:%=% generic)
nrdb.table: 
	python $(PARAM_DIR_SRC)fasta2nrdb.py \
		--verbose=1 \
		--tempdir=. \
		--input-filename-taxonomy=$(PARAM_DIR_INPUT)names.dmp \
		--input-filename-nrdb=$(PARAM_DIR_INPUT)nrdb.table \
		--output-filename-nrdb=$@ \
		$(PARAM_DIR_INPUT)pdbfinder.fasta generic \
		$(PARAM_DIR_INPUT)ncbi.fasta generic \
		$(PARAM_DIR_INPUT)refseq.fasta generic \
		$(PARAM_DIR_INPUT)swiss.fasta uniprot \
		$(PARAM_DIR_INPUT)trembl.fasta uniprot \
		$(PARAM_DIR_INPUT)wormpep.fasta generic \
		$(EXTRA_INPUT) \
	> $@.log

#############################################################
#############################################################
#############################################################
## build taxonomy table
#############################################################
taxonomy.table:
	python $(PARAM_DIR_SRC)Taxonomy.py \
		--action=build-pairsdb \
		--input-filename-names=$(PARAM_DIR_INPUT)names.dmp \
		--input-filename-nodes=$(PARAM_DIR_INPUT)nodes.dmp \
		--input-filename-division=$(PARAM_DIR_INPUT)division.dmp \
		--log=$@.log \
	> $@

#############################################################
#############################################################
#############################################################
## nrdb90
#############################################################
nrdb90.table: nrdb.table
	python $(PARAM_DIR_SRC)nrdb2rsdb.py \
		--verbose=1 \
		--tempdir=. \
		--level=90 \
		--filename-rsdb-table=nrdb90.table \
		--filename-rsdb-fasta=nrdb90.fasta \
		--filename-pairsdb-table=pairsdb_100x90.table \
		--memory=2000 \
	nrdb.fasta > $@.log

#############################################################
#############################################################
#############################################################
## 
#############################################################
nrdb.index: nrdb.table
	python $(PARAM_DIR_SRC)IndexedFasta.py \
		indexed_nrdb nrdb.fasta > $@

#############################################################
#############################################################
#############################################################
## 
#############################################################
pairsdb_100x90.table: nrdb90.table nrdb.index
	python $(PARAM_DIR_SRC)align_pairs.py \
		--filename-fasta=indexed_nrdb \
		--output-format=pairsdb \
		--output-filename-unaligned=$@.unaligned \
		--output-filename-notfound=$@.notfound \
		--log=$@.log \
        < pairsdb.table > $@

#############################################################
#############################################################
#############################################################
## compute masks
#############################################################

## variables for tmhmm
DIR_TMHMM_DATA=/usr/local/data/tmhmm
DIR_COILS_DATA=/usr/local/data/coils

##############################################################
## calculate various masks
nrdb90_masks_bias: nrdb90.table
	perl $(PARAM_DIR_SRC)biasdb.pl -rdb < nrdb90.fasta |\
	perl $(PARAM_DIR_SRC)bias_parser.pl 1 > $@

nrdb90_masks_tmhmm: nrdb90.table 
	perl -p -e "if (!/^>/) { s/[BJUO]/X/g; }" < nrdb90.fasta |\
	decodeanhmm -f $(DIR_TMHMM_DATA)/TMHMM1.0.options $(DIR_TMHMM_DATA)/TMHMM1.0.model |\
	perl $(PARAM_DIR_SRC)tmhmm_parser.pl 2 > $@

nrdb90_masks_coils: nrdb90.table
	export COILSDIR=$(DIR_COILS_DATA); \
	coils -f < nrdb90.fasta |\
	perl $(PARAM_DIR_SRC)mask_parser.pl 3 > $@

nrdb90_masks_short: nrdb90.table
	export PERL5LIB=$(PARAM_DIR_SRC); \
	perl $(PARAM_DIR_SRC)run_short.pl < nrdb90.fasta > $@

nrdb90_masks.table: nrdb90_masks_coils \
	nrdb90_masks_short \
	nrdb90_masks_tmhmm \
	nrdb90_masks_bias
	cat $^ > $@

#############################################################
#############################################################
#############################################################
## apply masks to nrdb90
#############################################################
nrdb90_masked.fasta: nrdb90_masks.table
	python $(PARAM_DIR_SRC)mask_sequences.py \
		--filename-masks=$< \
		--masks=all \
	< nrdb90.fasta > $@	

#############################################################
#############################################################
#############################################################
## build all-versus-all graph
#############################################################
pairsdb_90x90.table: nrdb90_masked.fasta
	formatdb -p T -n nrdb90_masked -i nrdb90_masked.fasta
	blastall -p blastp -d nrdb90_masked -i nrdb90_masked.fasta -e 1.0 -z 65000000 -b 100000 -v 100000 -F F -a 1|\
	perl $(PARAM_DIR_SRC)psiblastparser.pl --ends --table --log |\
	sort -T. -k1n -k2n -k3n > $@


#############################################################
#############################################################
#############################################################
## 
#############################################################
## Alignments from nrdb90 to nrdb40 are extracted from the
## blast graph. This uses a place holder script, use rsdb_create
## if possible.
pairsdb_90x40.table: pairsdb_90x90.table
	python $(PARAM_DIR_SRC)rsdb_create.py \
		--filename-graph=$< \
		--filename-fasta=indexed_nrdb \
	> $@.log

#############################################################
#############################################################
#############################################################
## 
#############################################################
pairsdb_100x40.table: pairsdb_100x90.table pairsdb_90x40.table
	python $(PARAM_DIR_SRC)map_alignments.py \
		--filename-fasta=indexed_nrdb \
		--allow-partial-map \
		--input-format=pairsdb \
		--output-format=pairsdb \
		--output-filename-unaligned=$@.unaligned \
		--output-filename-notfound=$@.notfound \
		--log=$@.log \
		pairsdb_100x90.table \
		pairsdb_90x40.table \
		> $@

#############################################################
#############################################################
#############################################################
## 
#############################################################
nrdb.sary:
	python $(PARAM_DIR_SRC)SaryFasta.py \
		sary_nrdb nrdb.fasta > $@

#############################################################
#############################################################
#############################################################
## Annotations: scop
#############################################################
scop.log: nrdb.sary pairsdb_100x40.table
	rm -f nrdb*_scop*
	python $(PARAM_DIR_SRC)Scop.py \
		--input-filename-sary=sary_nrdb \
		--input-filename-fasta=indexed_nrdb \
		--input-filename-annotations=$(PARAM_DIR_INPUT)scop.fasta \
		--input-filename-descriptions=$(PARAM_DIR_INPUT)scop.hierarchy \
		--input-filename-references=references.table \
	> $@


#############################################################
#############################################################
#############################################################
## Annotations: cath
#############################################################
cath.log: nrdb.sary pairsdb_100x40.table
	rm -f nrdb*_scop*
	python $(PARAM_DIR_SRC)Cath.py \
		--input-filename-sary=sary_nrdb \
		--input-filename-fasta=indexed_nrdb \
		--input-filename-annotations=$(PARAM_DIR_INPUT)cath.fasta \
		--input-filename-descriptions=$(PARAM_DIR_INPUT)cath.hierarchy \
		--input-filename-references=references.table \
	> $@


#############################################################
#############################################################
#############################################################
## Annotations: pdb
#############################################################
pdb.log: nrdb.sary pairsdb_100x40.table
	rm -f nrdb*_pdb*
	python $(PARAM_DIR_SRC)Pdb.py \
		--input-filename-fasta=indexed_nrdb \
		--input-filename-nrdb=nrdb.table \
		--input-filename-references=references.table \
	> $@

#############################################################
#############################################################
#############################################################
## Annotations: interpro
#############################################################
interpro.log: nrdb.sary pairsdb_100x40.table
	rm -f nrdb*_interpro*
	python $(PARAM_DIR_SRC)Interpro.py \
		--input-filename-fasta=indexed_nrdb \
		--input-filename-nrdb=nrdb.table \
		--input-filename-annotations=$(PARAM_DIR_INPUT)match_complete.xml.gz \
		--input-filename-descriptions=$(PARAM_DIR_INPUT)entry.list \
		--input-filename-references=references.table \
	> $@



