#!/usr/bin/env python

import os


#Sets up the relative paths to all database files
def init(vrtg_main):
	global vrtg_path
	vrtg_path = os.path.dirname(vrtg_main)

	global codon_tbl_path
	global blast_db_path
	global hmmer_db_path
	global acc2tax_refseq_path
	global acc2tax_vogs_path
	global tree_path
	global prot_tbl_path
	global acc2fun_path
	codon_tbl_path = 		os.path.join(vrtg_path, "db", "translation_tables", "codon_table.txt")
	blast_db_path = 		os.path.join(vrtg_path, "db", "refseq_blast", "blast_db")
	hmmer_db_path = 		os.path.join(vrtg_path, "db", "pvogs_hmm", "hmmer_db")
	acc2tax_refseq_path = 	os.path.join(vrtg_path, "db", "translation_tables", "acc2taxid_viral_refseq.txt")
	acc2tax_vogs_path = 	os.path.join(vrtg_path, "db", "translation_tables", "acc2taxid_vogs.txt")
	tree_path = 			os.path.join(vrtg_path, "db", "tax_tree", "viral_tax_tree.txt")
	prot_tbl_path = 		os.path.join(vrtg_path, "db", "pvogs_prot_tbl", "protein_table.txt")
	acc2fun_path = 			os.path.join(vrtg_path, "db", "translation_tables", "acc2func.txt")