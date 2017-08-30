#!/usr/bin/python

#This script downloads the pVOGs database and creates the correct 
#directory structure to hold those databases. By defaultatabases will be 
#download to the working directory, unless specified otherwise by the user
#using the appropriate flag.

#Import libraries
import argparse
import os
from shutil import rmtree, move
import tarfile
import subprocess
try:
    from urllib.request import urlopen # Python 3
except ImportError:
    from urllib2 import urlopen # Python 2


def get_user_input():

	parser = argparse.ArgumentParser(description = "This script creates " \
		"the directories for the pVOGs databases and download them.",
		epilog = "By default, the databases will be downloaded to the " \
		"working directory. Users can specify a different directory " \
		"using the appropriate flag.\nThe databases are quite large, " \
		"so it is recommended not run this script unless databases need " \
		"to be updated. See -h for more.")
	parser.add_argument("-d", metavar = "directory", dest = "directory", \
		help = "User specified directory for the databases.", \
		default = os.getcwd())
	parser.add_argument("-o", default = False, \
		action = "store_true", dest = "overwrite", help = "Overwrite the "\
		"current database. Specifying this flag will set overwrite to TRUE." \
		" <<<WARNING>>> The pVOGs databases may take some time to download " \
		"and extract.")
	user_args = parser.parse_args()
	path = os.path.join(user_args.directory, "db", "pVOGs")
	overwrite = user_args.overwrite
	return (path, overwrite)


#Check whether directories exists and, if not, create directories
def check_directories(path, overwrite):

	DATABASES = ["members", "proteins", "protein_table", 'VOGs', \
	"alignments", "HMMs"]
	missing_dbs = list()

	print("\nChecking database directories...\n")
	for db in DATABASES:
		directory = os.path.join(path, db)
		if check_directory(directory, overwrite):
			missing_dbs.append(db)
	print("Directories ready.\n\n%s\n" %(60*"="))
	return missing_dbs


def check_directory(directory, overwrite):

	db_missing = False

	try:
		os.makedirs(directory)
		print("Created %s\n" % (directory))
		db_missing = True
	except OSError:
		if overwrite:
			rmtree(directory)
			print("Overwriting directory '%s'.\n" %(directory))
			os.makedirs(directory)
			db_missing = True
		else:
			print("Directory '%s' already exists. This db " \
			"will not be downloaded. See -h for help.\n" % (directory))
	return db_missing

#Extract database, remove .gz file and move extracted files to correct db dir.
def extract_database(db, file, path):

	db_path = os.path.join(path, db)
	tar = tarfile.open(file)
	tar_members = tar.getnames()
	tar.extractall(path = db_path)
	tar.close()
	os.remove(file)
	for member in tar_members:
		if "." in member:		#Only move files, not directories
			move(os.path.join(db_path, member), db_path)
		else:
			old_dir = os.path.join(db_path, member)
	rmtree(old_dir)


def download_database(db, path, URLS):

	response = urlopen(URLS[db])
	file_type = URLS[db].split(".")[-1]
	file = os.path.join(path, db, db + "." + file_type)

	print("Downloading and extracting %s db..." % (db))

	with open(file, "wb") as f:
		chunk = response.read()
		f.write(chunk)
	if file_type == "gz":
		extract_database(db, file, path)

	print("Done.\n")


def get_databases(path, to_download):

	URLS = {"members" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/AllFamilyMembers.tsv",
	"proteins" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/AllFamilyProteinList.tsv",
	"protein_table" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/" \
	"VOGProteinTable.txt",
	"VOGs" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/Allvogtables.tar.gz",
	"alignments" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/Allvogsalignments.tar.gz",
	"HMMs" : "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All" \
	"/AllvogHMMprofiles.tar.gz"}
	
	for db in to_download:
		download_database(db, path, URLS)


def main():

	(path, overwrite) = get_user_input()
	to_download = check_directories(path, overwrite)
	get_databases(path, to_download)


if __name__ == '__main__':
	main()
