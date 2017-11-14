#!usr/bin/python

import os
import re
import sys


#Stores linked annotations
class Annotation:

	def __init__(self, acc, name, EC, KO, GO):

		self.acc = acc
		self.name = name
		self.EC = EC
		self.KO = KO
		self.GO = GO.replace(" ", "")

	def add_ann(self, ann_obj):

		for attr in ["acc", "name", "EC", "KO", "GO"]:
			if getattr(ann_obj, attr) != "":
				if getattr(self, attr) != "":
					setattr(self, attr, "%s;%s" % (getattr(self, attr), \
													getattr(ann_obj, attr)))
				else:
					setattr(self, attr, getattr(ann_obj, attr))


#Build the annotation table that converts NCBI Accession to
#functional annotations (GO, KEGG and EC)
def build_ann_table():

	#Read annotation table and store in dict based on NCBI accession
	ann_tab = dict()
	path = os.path.join("db", "translation_tables", "acc2func.txt")
	with open(path, "r") as file:
		for line in file:
			if line[0] != "#":
				data = line.split("\t")
				acc = data[1]
				name = data[2]
				EC = data[3]
				KO = data[4]
				GO = data[5]
				ann_tab[acc] = Annotation(acc, name, EC, KO, GO)
	return(ann_tab)


#Build the table that converts VOG ids to NCBI accessions
def build_vog_table():

	#A subset of species contains a ":" in their name, but ":" is also a split-
	#character. 
	def sub_colon(line):

		colon = re.search("phi[0-9]+.?.?:[0-9]", line)
		if colon:
			line = line.replace(colon.group(0), \
			colon.group(0).replace(":", "="))
		return(line)

	#Read protein table and store translation from VOG id to NCBI accessions
	vog_tab = dict()
	path = os.path.join("db", "pVOGs", "protein_table", "protein_table.txt")
	with open(path, "r") as file:
		for line in file:
			if line[0] != "#":
				line = sub_colon(line)
				data = line.split(":")
				id = data[0]
				acc = data[3].split("|")[0].split("-")[1]
				if id not in vog_tab.keys():
					vog_tab[id] = list()
				vog_tab[id].append(acc)
	return(vog_tab)


#Writes mapped functions based on BLAST output
def get_funcs_blast(output_file, orfF, ann_tab):

	blast_path = output_file + ".blast"
	func_path = output_file + ".orf.fun"
	prev_query = None
	with open(blast_path, "r") as in_file:
		with open(func_path, "w") as out_file:
			#Write header
			out_file.write("#ORF\tNCBI_accessions\tgene_names\tEC_numbers\tKO_ids\tGO_ids\n")
			for line in in_file:
				data = line.split("\t")
				query = data[0]
				acc = data[1]
				try:
					ann = ann_tab[acc]
				except KeyError:
					ann = Annotation(acc, "hypothetical protein", "", "", "")
				if query != prev_query:
					if prev_query != None: #Skip empty
						write_funcs(prev_query, total_ann, out_file)
					prev_query = query
					total_ann = Annotation("", "", "", "", "")
					total_ann.add_ann(ann)
				if orfF == "all":
					total_ann.add_ann(ann)
			#Write the last one to file
			write_funcs(prev_query, total_ann, out_file)


#Writes mapped functions based on HMMer output
def get_funcs_hmmer(output_file, orfF, ann_tab, vog_tab, to_write):

	hmmer_path = output_file + ".hmmer"
	func_path = output_file + ".orf.fun"
	prev_query = None
	with open(hmmer_path, "r") as in_file:
		with open(func_path, to_write) as out_file:
			#Write header
			if to_write == "w":
				out_file.write("#ORF\tNCBI_accessions\tgene_names\tEC_numbers\tKO_ids\tGO_ids\n")
			for line in in_file:
				if line[0] != "#":
					data = line.split()
					query = data[2]
					vog_id = data[0]
					if query != prev_query:
						if prev_query != None: #Skip empty
							write_funcs(prev_query, total_ann, out_file)
						prev_query = query
						total_ann = Annotation("", "", "", "", "")
						for acc in vog_tab[vog_id]:
							try:
								ann = ann_tab[acc]
							except KeyError:
								ann = Annotation(acc, "hypothetical protein", "", "", "")
							total_ann.add_ann(ann)
					
					if orfF == "all":
						for acc in vog_tab[vog_id]:
							try:
								ann = ann_tab[acc]
							except KeyError:
								ann = Annotation(acc, "hypothetical protein", "", "", "")
							total_ann.add_ann(ann)

			#Write the last one to file
			try:
				write_funcs(prev_query, total_ann, out_file)
			except UnboundLocalError: #this means there were no HMMer annotations
				pass


#Get functional annotation output tables for BLAST and HMMer output
def get_functions(output_file, orfF, soft):

	ann_tab = build_ann_table()
	vog_tab = build_vog_table()
	if soft == "blast" or soft == "both":
		get_funcs_blast(output_file, orfF, ann_tab)
	if soft == "hmmer" or soft == "both":
		to_write = "w"
		if soft == "both":
			to_write = "a"
		get_funcs_hmmer(output_file, orfF, ann_tab, vog_tab, to_write)


def write_funcs(query, ann_obj, out_file):

	out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (query, ann_obj.acc, ann_obj.name, \
											ann_obj.EC, ann_obj.KO, ann_obj.GO))


#The main function is only used for testing.
def main():

	os.chdir("..")
	orfF = "best_hit"
	output_file = "output/lambda_virus"
	soft = "both"
	get_functions(output_file, orfF, soft)


if __name__ == '__main__':
	main()