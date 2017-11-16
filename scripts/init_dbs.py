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
import scripts.settings as settings


#Append new indices to acc2tax table
def append_acc2tax(acc2tax):

	sys.stdout.write("Added %d entries\n" % len(acc2tax))
	with open(settings.acc2tax_refseq_path, "a") as out_file:
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
def download_acc2tax(missing):

	#store accessions as strings of 250 entries
	id_calls= list()
	ids = list()
	count = 0
	for acc in missing:
		ids.append(acc)
		count += 1
		if count == 200:
			id_calls.append("&id=".join(ids))
			ids = list()
			count = 0
	if count != 200:
		id_calls.append("&id=".join(ids))

	#Loop through id_calls and fetch xml output
	xmls = str()
	for call in id_calls:
		url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?" \
				"dbfrom=protein&db=taxonomy&id=%s" % call
		response = urlopen(url)
		chunk = response.read()
		xmls += chunk

	#Get all the taxids
	pattern = "<Link>\n.*<Id>(.*)</Id>"
	taxids = re.findall(pattern, xmls)

	#Store in acc2tax dict
	acc2tax = dict()
	for i in range(len(missing)):
		acc = missing[i]
		taxid = taxids[i]
		acc2tax[acc] = taxid

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
	relevant = ["nodes.dmp", "names.dmp"]
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


#Compares all accessions in the downloaded Viral RefSeq with those found
#in acc2tax file and returns the missing entries
def find_missing_acc2tax():

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
		print("\nDownloading and extracting %s: " % (db))
		data = str() #Variable to store data spread over multiple files
		path = databases[db]
		for url in urls[db]:
			file = download_file(db, path, url, data)
			data = extract_file(db, file, path, data)
			sys.stdout.write("Done\n")

		if data:
			with open(path, "w") as out_file:
				out_file.write(data)

	#Create the viral_tax_tree.txt file according to the accession numbers in RefSeq Viral
	sys.stdout.write("\nCreating taxonomy tree...\n")
	make_tree_file()

	#Add mossing acc2taxids
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
	os.remove(nodes_file)
	os.remove(names_file)
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

#Creates a file that contains all information about the Viral Tax Tree
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
			if div in [3, 9]:	#Divs 3 and 9 represent phages and viruses
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

#Check which accession2tax are missing from the file and add them through
#NCBI Entrez calls
def update_acc2tax():

	missing = find_missing_acc2tax()
	acc2tax = download_acc2tax(missing)
	append_acc2tax(acc2tax)



def main():

	init_dbs()


if __name__ == '__main__':
	main()
