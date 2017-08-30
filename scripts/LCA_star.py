#!usr/bin/python

import os
import re
import sys
from math import log
from copy import copy

#Import classes
from scripts.classes.class_node import Node
from scripts.classes.class_lca import LCA


#Build the accesion to taxid conversion table
def build_acc_taxid():

	taxids = dict()
	path = os.path.join("db", "translation_tables", "acc2taxid.txt")
	with open(path, "r") as file:
		for line in file:
			if line[0] != "#":
				line = line.strip().split()
				acc = line[0]
				taxid = line[1]
				if line[0] not in taxids.keys():
					taxids[acc] = taxid
	return(taxids)


#Build the VOG to taxid conversion table
def build_vog_taxid(taxids):

	vog_to_taxids = dict()
	path = os.path.join("db", "pVOGs", "protein_table", "protein_table.txt")
	with open(path, "r") as file:
		for line in file:
			if line[0] != "#":
				data = line.strip().split(":")
				vog_id = data[0]
				acc = data[1]
				taxid = taxids[data[1]]
				if vog_id not in vog_to_taxids.keys():
					vog_to_taxids[vog_id] = list()
				vog_to_taxids[vog_id].append(taxid)
	return(vog_to_taxids)


#Build the NCBI tax tree (for viruses) based on a custom file.
def build_ncbi_tree():

	#Read the file
	path = os.path.join("db", "tax_tree", "viral_tax_tree.txt")
	path = os.path.normpath(path)
	file = open(path, "r")
	data = file.readlines()
	tree = dict()

	#Store every line as a node
	for line in data:
		if line[0] != "#":
			node = Node(line)
			node.parse()
			if node.id != "NA":
				tree[node.id] = node

	#Point to the corresponding parent and child nodes
	to_pop = list()
	for id in tree.keys():

		#point to parents
		parent_id = tree[id].parent_id
		if parent_id == "1": #This is the "root" (of all evil)
			continue
		try:
			tree[id].parent = tree[parent_id]
		except KeyError:
			to_pop.append(id) #Non-viral node, remove.
			continue

		#point parent to node
		child_id = id
		tree[parent_id].children[child_id] = tree[child_id]

	for id in to_pop:
		tree.pop(id)

	#update "no rank" assignments based on annotations of other children.
	for id in tree.keys():
		if tree[id].level == "no rank":
			if tree[id].parent.level == "superkingdom":
				tree[id].level = "genome type"
			elif tree[id].parent.level == "family":
				if "subfamily" in tree[id].parent.get_child_levels():
					tree[id].level = "subfamily"
				elif "genus" in tree[id].parent.get_child_levels():
					tree[id].level = "genus"
				elif "species" in tree[id].parent.get_child_levels():
					tree[id].level = "species"
			elif tree[id].parent.level == "genus":
				if "species" in tree[id].parent.get_child_levels():
					tree[id].level = "subgenus"
		tree[id].calc_index()

	return(tree)


#Do the LCA* calculations and write to output file.
def calc_write_lca(lca_obj, out_file, tree, n_orfs = False):

	lca_taxid = lca_obj.calc_lca(n_orfs)
	p_value = lca_obj.calc_p()
	lca_obj.clean_tree()
	write_lca(lca_obj.id, lca_taxid, p_value, out_file, tree)
	return(lca_taxid)


#Calculate the LCA of every predicted ORF based on BLAST out.
def get_lca_blast(output_file, taxids, tree, alpha):

	contig_taxids = dict() 	#stores assigned ORF LCA taxids per contig id
	with open(output_file + ".blast", "r") as blast_out:
		with open(output_file + ".orf.tax", "w") as tax_out:
			
			#Write file header
			tax_out.write("%s\t%s\t%s\t%s\n" % ("#orf_id", "LCA_taxid", \
												"p_value", \
												"\t".join(tree["10239"].tax_levels)))
			
			#Loop through BLAST out and calculate LCA for every predicted ORF
			orf_lca = LCA(None, tree, alpha)
			for line in blast_out:

				#Get basic data
				data = line.split("\t")
				query_id = data[0]
				query_data = data[0].split("|")
				target_ids = data[1].split("|")
				contig_id = "|".join(query_data[0:len(query_data) - 1])
				orf_id = query_data[-1]
				acc = target_ids[1].split("-")[0]
				taxid = taxids[acc]

				#Create contig_taxids list
				if contig_id not in contig_taxids.keys():
					contig_taxids[contig_id] = list()

				#Analyse and remove the previous LCA object if a new ORF is 
				#encountered in BLAST out.
				if query_id != orf_lca.id:
					if orf_lca.id: 		#Skip the empty LCA object
						lca_taxid = calc_write_lca(orf_lca, tax_out, tree)
						contig_taxids[contig_id].append(lca_taxid)
					orf_lca = LCA(query_id, tree, alpha)

				#Add taxids to LCA object
				if taxid != "NA":
					orf_lca.taxids.append(taxid)
			
			#Last ORF
			if orf_lca.id: 		#Skip the empty LCA object
				lca_taxid = calc_write_lca(orf_lca, tax_out, tree)
				contig_taxids[contig_id].append(lca_taxid)

			return(contig_taxids)


#Calculates the LCA of every contig based on LCAs of predicted ORFs.
def get_lca_contig(output_file, contig_taxids, tree, alpha, orf_count):

	with open(output_file + ".con.tax", "w") as tax_out:
		tax_out.write("%s\t%s\t%s\t%s\n" % ("#contig_id", "NCBI_TaxID", \
												"p_value", \
												"\t".join(tree["10239"].tax_levels)))
		for id in contig_taxids.keys():
			n_orfs = orf_count[id]
			contig_lca = LCA(id, tree, alpha)
			contig_lca.taxids = contig_taxids[id]
			calc_write_lca(contig_lca, tax_out, tree, n_orfs)


#Calculate the LCA of every predicted ORF based on HMMer out.
def get_lca_hmmer(output_file, tree, alpha, vog_to_taxids, to_write):

	contig_taxids = dict() 	#stores assigned ORF LCA taxids per contig id
	with open(output_file + ".hmmer", "r") as hmmer_out:
		with open(output_file + ".orf.tax", to_write) as tax_out:
			
			if to_write == "w":
				#Write file header
				tax_out.write("%s\t%s\t%s\t%s\n" % ("#orf_id", "LCA_taxid", \
												"p_value", \
												"\t".join(tree["10239"].tax_levels)))
					
			#Loop through BLAST out and calculate LCA for every predicted ORF
			orf_lca = LCA(None, tree, alpha)
			for line in hmmer_out:

				if line[0] != "#":
					#Get basic data
					data = line.split()
					query_id = data[2]
					query_data = data[2].split("|")
					target_id = data[0]
					contig_id = "|".join(query_data[0:len(query_data) - 1])
					orf_id = query_data[-1]
					taxids = vog_to_taxids[target_id]

					#Create contig_taxids list
					if contig_id not in contig_taxids.keys():
						contig_taxids[contig_id] = list()

					#Analyse and remove the previous LCA object if a new ORF is 
					#encountered in HMMer out.
					if query_id != orf_lca.id:
						if orf_lca.id: 		#Skip the empty LCA object
							lca_taxid = calc_write_lca(orf_lca, tax_out, tree)
							contig_taxids[contig_id].append(lca_taxid)
						orf_lca = LCA(query_id, tree, alpha)

					#Add taxids to LCA object
					for taxid in taxids:
						if taxid != "NA":
							orf_lca.taxids.append(taxid)
			
			#Last ORF
			if orf_lca.id:
				lca_taxid = calc_write_lca(orf_lca, tax_out, tree)
				contig_taxids[contig_id].append(lca_taxid)
			return(contig_taxids)


#Calculates the LCA-star based on BLAST and HMMer output for ORFs and contigs.
def lca_star(output_file, alpha, soft, orf_alg, cont_alg, orf_count):

	#Build required translation tables and NCBI taxon tree
	tree = build_ncbi_tree()
	taxids = build_acc_taxid()
	vog_to_taxids = build_vog_taxid(taxids)

	#Set alphas
	orf_alpha = alpha
	cont_alpha = alpha
	if orf_alg == "lca":
		orf_alpha = 0.99999 #LCA looks for a parent node shared by all children
	if cont_alg == "lca":
		cont_alpha = 0.99999

	#Calculate ORF LCAs
	if soft == "blast" or soft == "both":
		contig_taxids_blast = get_lca_blast(output_file, taxids, tree, orf_alpha)
		contig_taxids = contig_taxids_blast
	if soft == "hmmer" or soft == "both":
		to_write = "a" #Append to file if BLAST output has been processed
		if soft == "hmmer":
			to_write = "w" #Write to file if only using HMMer
		contig_taxids_hmmer = get_lca_hmmer(output_file, tree, orf_alpha, vog_to_taxids, to_write)
		contig_taxids = contig_taxids_hmmer
	if soft == "both":
		contig_taxids = merge_dicts(contig_taxids_blast, contig_taxids_hmmer)

	#Calculate contig LCAs
	get_lca_contig(output_file, contig_taxids, tree, cont_alpha, orf_count)


#Merge two dictionaries
def merge_dicts(a, b):
	c = dict()
	for id in a.keys():
		c[id] = copy(a[id])
		if id in b.keys():
			c[id] += b[id]
	for id in b.keys():
		if id not in c.keys():
			c[id] = copy(b[id])
	return(c)


#Write information about the LCA and its ancestry to file
def write_lca(query_id, lca_taxid, p_value, file, tree):

	tax = list()
	lineage = tree[lca_taxid].return_parents()
	for level in tree["10239"].tax_levels:
		if level in lineage.keys():
			tax.append(lineage[level])
		else:
			tax.append("NA")
	tax = "\t".join(tax)
	file.write("%s\t%s\t%.2E\t%s\n" % (query_id, lca_taxid, \
										p_value, tax))


#The main function is only used for testing.
def main():

	os.chdir("..")
	output_file = "output/lambda_virus"
	alpha = 0.5
	soft = "both"
	orf_alg = "lca_star"
	cont_alg = "lca_star"
	orf_count = {'gi|9626243|ref|NC_001416.1|__Enterobacteria__phage__lambda,__complete__genome': 60}
	lca_star(output_file, alpha, soft, orf_alg, cont_alg, orf_count)


if __name__ == '__main__':
	main()