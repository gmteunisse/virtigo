#!/usr/bin/env python

#This script downloads the databases that Virtigo requires to run.

#Import libraries
import argparse
import os
import sys
from shutil import rmtree, move
import tarfile
import gzip
import subprocess
from time import sleep
import re
try:
    from urllib.request import urlopen # Python 3
except ImportError:
    from urllib2 import urlopen # Python 2

#Import settings
import modules.settings as settings


#Append new indices to acc2tax table
def append_acc2tax(node_merges, acc2tax, db):

	if db == "RefSeq":
		path = settings.acc2tax_refseq_path
	elif db == "pVOGs":
		path = settings.acc2tax_vogs_path

	sys.stdout.write("Added %d %s entries\n" % (len(acc2tax), db))
	merged_nodes = node_merges.keys()
	lines = list()
	with open(path, "r") as in_file:
		in_file.readline()
		for line in in_file:
			data = line.strip().split()
			acc = data[0]
			taxid = data[1]
			if taxid in merged_nodes:
				taxid = node_merges[taxid]
				line = "%s\t%s\n" % (acc, taxid)
			lines.append(line)
	with open(path, "w") as out_file:
		out_file.write("#Accession\tTaxid\n")
		for line in lines:
			out_file.write(line)
		for acc in acc2tax.keys():
			out_file.write("%s\t%s\n" % (acc, acc2tax[acc]))


#Check whether database files exists
def check_dbs():

	databases = {"RefSeq Viral BLAST DB" : settings.blast_db_path + ".phr", \
				"pVOGs HMMER DB" : settings.hmmer_db_path + ".h3f", \
				"pVOGs protein table" : settings.prot_tbl_path, \
				"NCBI taxonomy tree" : settings.tree_path}
	missing_dbs = list()

	sys.stdout.write("\nChecking databases...\n")
	for db in databases.keys():
		if not os.path.exists(databases[db]):
			missing_dbs.append(db)

	missing = "\n".join(missing_dbs)

	if missing_dbs:
		exit_message = "Error: the following database files are missing:\n" \
						"%s\nPlease initialize the databases by running " \
						"Virtigo with the -i option." % missing
		sys.exit(exit_message)
	else:
		sys.stdout.write("Databases okay.\n")


#Creates directories to store a file in
def create_dir(path):

	if not os.path.exists(os.path.dirname(path)):
		os.makedirs(os.path.dirname(path))


#Downloads and extracts a file to a path
def download_file(db, path, url, data):

	create_dir(path)
	response = urlopen(url)
	if ".tar.gz" in url:
		file_type = ".tar.gz"
	elif ".gz" in url:
		file_type = ".gz"
	else:
		file_type = str()

	file = path + file_type

	#Get the file size for progress report
	size = float(response.info().getheader("Content-Length").strip()) / 1024**2
	chunk_size = int(8192)
	downloaded = 0

	#Download the file in chunks and report the progress
	with open(file, "wb") as f:
		while 1:
			chunk = response.read(chunk_size)
			downloaded += float(chunk_size) / float(1024**2)
			if not chunk:
				sys.stdout.write("Downloaded %0.2f%% of %.2f MB.\n" % \
							(100, size))
				break
			sys.stdout.write("Downloaded %0.2f%% of %.2f MB.\r" % \
							(float(downloaded / size) * 100, size))
			f.write(chunk)

	sys.stdout.write("Extracting...\n")
	return(file)
	#Extract file if it is in a .gz archive
	if file_type == ".tar.gz":
		data += extract_tar(file, path)
	#Extract .gz files and return content
	if file_type == ".gz":
		data += extract_gz(file, path)
	return(data)


#Use NCBI Eutils to obtain the missing indices
def download_acc2tax(missing, db):

	pattern = "<Link>\n.*<Id>(.*)</Id>"
	acc2tax = dict()

	if db == "refseq":
		eutils_db = "protein"
	elif db == "pvogs":
		eutils_db = "nuccore"

	for acc in missing:
		url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?" \
				"dbfrom=%s&db=taxonomy&id=%s" % (eutils_db, acc)
		response = urlopen(url)
		chunk = response.read()
		taxid = re.findall(pattern, chunk)
		if len(taxid) < 1:
			taxid = ["10239"] #Viral root
		elif len(taxid) > 1:
			taxid = taxid_in_tree(taxid)
		acc2tax[acc] = taxid[0]

	return(acc2tax)
		

#Extract the file according to its type
def extract_file(db, file, path, data):

	if db == "RefSeq Viral":
		data += extract_refseq(file, path)
	elif db == "pVOGs HMMs":
		data += extract_hmms(file, path)
	elif db == "NCBI taxonomy tree":
		extract_tree(file, path)
	return(data)

#Extract .tar.gz and remove file, and only keep relevant files
def extract_tree(file, path):

	tmp_path = path + "_tmp"
	tar = tarfile.open(file)
	tar_members = tar.getnames()
	tar.extractall(path = tmp_path)
	tar.close()
	os.remove(file)
	relevant = ["nodes.dmp", "names.dmp", "merged.dmp"]
	for member in tar_members:
		if member in relevant:
			move(os.path.join(tmp_path, member), os.path.dirname(path))
	rmtree(tmp_path)


#Get contents from .gz file and remove .gz file
def extract_refseq(file, path):

	with gzip.open(file) as f:
		content = f.read()
	os.remove(file)
	return(content)


#Extract tar.gz, remove .tar.gz and extracted files and return contents
def extract_hmms(file, path):

	tmp_path = path + "_tmp"
	tar = tarfile.open(file)
	tar_members = tar.getnames()
	tar.extractall(path = tmp_path)
	tar.close()
	os.remove(file)
	data = str()
	for member in tar_members:
		if "." in member:		#Only read files, not directories
			data += open(os.path.join(tmp_path, member), "r").read()
	rmtree(tmp_path)
	return(data)


#Download the databases and create files for blast and hmmer
def init_dbs():

	#Set paths and URLS
	databases = {"RefSeq Viral" : settings.blast_db_path, \
				"pVOGs HMMs" : settings.hmmer_db_path, \
				"pVOGs protein table" : settings.prot_tbl_path, \
				"NCBI taxonomy tree" : settings.tree_path}

	urls = {"RefSeq Viral" : \
			["ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz",
			"ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz"],
	"pVOGs protein table" : ["http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/" \
	"VOGProteinTable.txt"],
	"pVOGs HMMs" : ["http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/AllvogHMMprofiles.tar.gz"], \
	"NCBI taxonomy tree" : ["ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"]}

	#Download and extract files
	for db in databases.keys():
		sys.stdout.write("\nDownloading and extracting %s: \n" % (db))
		data = str() #Variable to store data spread over multiple files
		path = databases[db]
		for url in urls[db]:
			file = download_file(db, path, url, data)
			data = extract_file(db, file, path, data)
			sys.stdout.write("Done\n")

		if data:
			with open(path, "w") as out_file:
				out_file.write(data)

	#Create the viral_tax_tree.txt file by removing all non-viral nodes
	sys.stdout.write("\nCreating taxonomy tree...\n")
	make_tree_file()

	#Add missing acc2taxids
	sys.stdout.write("Done\nAdding missing NCBI Accession-to-TaxID entries...\n")
	update_acc2tax()

	#Create HMMER and BLAST DBs
	sys.stdout.write("Done\nCreating HMMER DB...\n")
	make_hmmer_db()
	sys.stdout.write("Done\nCreating BLAST DB...\n")
	make_blast_db()
	sys.stdout.write("Done\n")

	#Remove files
	sys.stdout.write("Cleaning up...\n")
	os.remove(settings.hmmer_db_path)
	os.remove(settings.blast_db_path)
	nodes_file = os.path.join(os.path.dirname(settings.tree_path), "nodes.dmp")
	names_file = os.path.join(os.path.dirname(settings.tree_path), "names.dmp")
	merged_file = os.path.join(os.path.dirname(settings.tree_path), "merged.dmp")
	os.remove(nodes_file)
	os.remove(names_file)
	os.remove(merged_file)
	sys.stdout.write("Done\n")


def make_blast_db():

	path = settings.blast_db_path
	call = 	"makeblastdb -in %s -title blast_db -out %s -dbtype prot" % (path, path)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if error:
		sys.stderr.write("BLAST error:\n%s" % (error.decode("ascii")))


def make_hmmer_db():

	path = settings.hmmer_db_path
	call = 	"hmmpress -f %s" % (path)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if error:
		sys.stderr.write("HMMer error:\n%s" % (error.decode("ascii")))

#Creates a file that contains only viral nodes of the NCBI Tax Tree
def make_tree_file():

	class Node:
		def __init__(self):
			self.tax_id = int()
			self.parent_tax_id = int()
			self.rank = str()
			self.name = str()

	#Subset to viral nodes
	nodes = dict()
	path = os.path.join(os.path.dirname(settings.tree_path), "nodes.dmp")
	with open(path, "r") as nodes_file:
		for line in nodes_file:
			data = line.split("\t|\t")
			div = int(data[4])
			if div in [3, 9, 11]:	#3: phages; 9: viruses; 11: environment
				viral_node = Node()
				viral_node.tax_id = int(data[0])
				viral_node.parent_tax_id = int(data[1])
				viral_node.rank = data[2]
				nodes[viral_node.tax_id] = viral_node

	#Get names of viral nodes
	path = os.path.join(os.path.dirname(settings.tree_path), "names.dmp")
	ids = list(nodes.keys())
	ids = sorted(ids)
	min_id = min(ids)
	max_id = max(ids)

	with open(path, "r") as names_file:
		for line in names_file:
			data = line.split("\t|\t")
			tax_id = int(data[0])
			if tax_id < min_id: 	#Skip lines with too low an ID
				continue
			if tax_id == ids[0]:
				if "scientific name" in data[3].lower():
					nodes[tax_id].name = data[1]
					ids.remove(tax_id)
			if len(ids) < 1 or tax_id > max_id:
					break

	#Write to file
	ids = list(nodes.keys())
	ids = sorted(ids)
	with open(settings.tree_path, "w") as out_file:
		out_file.write("tax_id\tparent_tax_id\trank\tname\n")
		for id in ids:
			out_file.write("%d\t%d\t%s\t%s\n" % (nodes[id].tax_id, \
												nodes[id].parent_tax_id, \
												nodes[id].rank, \
												nodes[id].name))


#Compares all accessions in the downloaded Viral RefSeq with those found
#in acc2tax file and returns the missing entries
def missing_acc2tax_pvogs():

	#Get all the accession numbers in acc2tax
	accs_tax = list()
	with open(settings.acc2tax_vogs_path, "r") as acc2tax:
		acc2tax.readline() #Skip first
		for line in acc2tax:
			data = line.split("\t")
			accs_tax.append(data[0])

	#Create dict to replace segmentphis with accessions
	phis = {"segmentsphi12" :	"NC_004175.1", \
			"segmentsphi13" :	"NC_004171.1", \
			"segmentsphi2954" :	"NC_012093.1", \
			"segmentsphi6":		"NC_003715.1", \
			"segmentsphi8" :	"NC_003299.1", \
			"segmentsphiNN" :	"KJ957164.1"}

	#Loop throug pvogs and store accessions
	accs_pvogs = list()
	with open(settings.prot_tbl_path, "r") as pvogs:
		for line in pvogs:
			data = line.split(":")
			acc = data[1]
			if acc in phis.keys():
				acc = phis[acc]
			accs_pvogs.append(acc)

	#Sort for fast searching
	accs_tax = sorted(accs_tax)
	accs_pvogs = sorted(list(set(accs_pvogs)))

	#Loop through accesssions and find missing
	missing = list()
	if len(accs_tax) < 1:
		missing = accs_pvogs
		return(missing)

	for acc in accs_pvogs:
		if acc == accs_tax[0]:
			accs_tax.remove(acc)
		else:
			missing.append(acc)

	return(missing)


#Compares all accessions in the downloaded Viral RefSeq with those found
#in acc2tax file and returns the missing entries
def missing_acc2tax_refseq():

	#Get all the accession numbers in acc2tax
	accs_tax = list()
	with open(settings.acc2tax_refseq_path, "r") as acc2tax:
		acc2tax.readline() #Skip first
		for line in acc2tax:
			data = line.split("\t")
			accs_tax.append(data[0])

	#Loop throug RefSeq and store accessions
	accs_refseq = list()
	with open(settings.blast_db_path, "r") as refseq:
		for line in refseq:
			if line[0] == ">":
				acc = line.split()[0][1:]
				accs_refseq.append(acc)

	#Sort for fast searching
	accs_tax = sorted(accs_tax)
	accs_refseq = sorted(accs_refseq)

	#Loop through accesssions and find missing
	missing = list()
	for acc in accs_refseq:
		if acc == accs_tax[0]:
			accs_tax.remove(acc)
		else:
			missing.append(acc)

	return(missing)

#Returns a dict of merged nodes in the NCBI Tax tree
def read_merged_nodes():

	merged_file = os.path.join(os.path.dirname(settings.tree_path), "merged.dmp")
	node_merges = dict()
	with open(merged_file) as in_file:
		for line in in_file:
			data = line.split("|")
			old_id = data[0].strip()
			new_id = data[1].strip()
			node_merges[old_id] = new_id
	return(node_merges)


#When multiple TaxIDs are returned, pick the one
#that has a node in the tree
def taxid_in_tree(taxids):

	with open(settings.tree_path) as tree:
		tree.readline() #Skip header
		for line in tree:
			data = line.split()
			id = data[0]
			if id in taxids:
				taxid = [str(id)]
				return(taxid)
	taxid = ["10239"] 	#Nothing found, so assign viral root
	return(taxid)



#Check which accession2tax are missing from the file and add them through
#NCBI Entrez calls
def update_acc2tax():

	missing_refseq = missing_acc2tax_refseq()
	missing_pvogs = missing_acc2tax_pvogs()
	acc2tax_refseq = download_acc2tax(missing_refseq, "refseq")
	acc2tax_pvogs = download_acc2tax(missing_pvogs, "pvogs")
	node_merges = read_merged_nodes()
	append_acc2tax(node_merges, acc2tax_refseq, "RefSeq")
	append_acc2tax(node_merges, acc2tax_pvogs, "pVOGs")


def main():

	init_dbs()


if __name__ == '__main__':
	main()
