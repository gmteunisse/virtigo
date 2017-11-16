#!/usr/bin/env python

#Import libraries
import argparse
import multiprocessing
import subprocess
import os
import sys
import distutils


#Import modules
from scripts.init_dbs import check_dbs, init_dbs
from scripts.predict_orfs import predict_orfs, update_contig_ids
from scripts.extract_orfs import extract_orfs
from scripts.BLASTp import blast_search, parse_blast_result
from scripts.HMMer import hmmer_search, parse_hmmer_result, hmmer_best_hit
from scripts.map_reads import map_reads
from scripts.LCA_star import lca_star
from scripts.func_annotation import get_functions
import scripts.settings as settings

#Define global variable: path to virtigo software
settings.init(__file__)


#Check whether a bin exists in the '$PATH' environmental variable. Courtesy of
#https://stackoverflow.com/questions/377017/
#test-if-executable-exists-in-python/377028#377028
def which(program):
		
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	for path in os.environ["PATH"].split(os.pathsep):
		path = path.strip('"')
		exe_file = os.path.join(path, program)
		if is_exe(exe_file):
			return exe_file
	return None


#check whether all software has been installed and is callable.
def check_software():

	software = ["blastp", "hmmscan"]
	for tool in software:
		avail = which(tool)
		if not avail:
			exit_message = "Error: no %s installation found. Make sure to " \
			"install %s and add the path to your $PATH variable." % (tool, tool)
			sys.exit(exit_message)

	software = ["samtools", "genomeCoverageBed"]
	for tool in software:
		avail = which(tool)
		if not avail:
			sys.stderr.write = "Error: no %s installation found. Read " \
			"mapping and abundance estimation will not function. To " \
			"prevent errors, install %s and add the path to your $PATH " \
			"variable." % (tool, tool)

#Check whether required input has been provided
def check_input(user_args):
	
	if not user_args.contigs:
		exit_message = "\nError: no contigs file provided. See -h.\n"
		sys.exit(exit_message)
	if not user_args.output:
		exit_message = "\nError: no output filename provided. See -h.\n"
		sys.exit(exit_message)

	#Convert to abspaths
	user_args.contigs = os.path.abspath(user_args.contigs)
	user_args.output = os.path.abspath(user_args.output)
	if user_args.reads:
		user_args.reads = os.path.abspath(user_args.reads)
	return(user_args)


#Get the OS so that the correct bins are called
def get_bins():

	os_type = sys.platform.lower()
	if "darwin" in os_type:		#OSX
		bin_path = os.path.join(settings.vrtg_path, "bin", "osx")
	elif "linux" in os_type:
		bin_path = os.path.join(settings.vrtg_path, "bin", "linux")
	else:
		print("OS not recognized, assuming linux.")
		bin_path = os.path.join(settings.vrtg_path, "bin", "linux")
	bin_path = os.path.abspath(bin_path)
	return(bin_path, os_type)


#If no number of threads has been specified, use the maximum number of threads
#available.
def get_n_threads(args):

	if args.threads == 0:
		args.threads = multiprocessing.cpu_count()
	return(args)


#Obtain the user arguments
def get_user_args():

	parser = argparse.ArgumentParser(description = "This script annotates " \
		"a file of presumed viral contiguous sequences.",
		epilog = "This script will only run on Linux and OSX systems. " \
		"Required dependencies: " \
		"\n\t BEDtools (tested with v2.26);" \
		"\n\t BLAST+ (tested with v2.6.0);" \
		"\n\t Samtools (tested with v1.5);" \
		"\n\t HMMer (tested with v3.1);" \
		"\n\t SciPy (tested with v0.18.1).")

	parser.add_argument("-a", metavar = "alpha", dest = "alpha", \
		help = "Majority criterium for LCA* algorithm (default = 0.5).", \
		default = 0.5, type = float)

	parser.add_argument("-c", metavar = "input_contigs", dest = "contigs", \
		help = "Required. Path to (multi-)FASTA file containing contigs. " \
		"Entire header will be used as the contig ID. Sequences need to be "\
		"nucleotide sequences.")

	parser.add_argument("--contA", metavar = "cont_tax", dest = "cont_alg", \
		help = "Algorithm to use to determine taxonomy of contigs. " \
		"Use one of 'lca' or 'lca_star' (default = lca_star), " \
		"which respectively assign the naive LCA as in MEGAN and the LCA*.", \
		default = "lca_star", type = str, \
		choices = ["lca", "lca_star"])

	parser.add_argument("-e", metavar = "e_value", dest = "e_value", \
		help = "E-value threshold for BLASTp and HMMer (default = 0.00001).", \
		default = 0.00001, type = float)

	parser.add_argument("--err", dest = "err", \
		help = "Redirect stderr to <output>.err.", default = False,
		action = "store_true")

	parser.add_argument("-i", dest = "init", \
		help = "Initializes Virtigo by downloading the required databases " \
		"and generating searchable files. An active internet connection and " \
		"~2 GB of disk space are required." , default = False,
		action = "store_true")

	parser.add_argument("-n", metavar = "n_threads", dest = "threads", \
		help = "Number of threads to use for BLASTp and HMMer. " \
		"Default (0) will use all detected threads.", \
		default = 0, type = int)

	parser.add_argument("-o", metavar = "output", dest = "output", \
		help = "Required. Path and filename (without extension) for output " \
		"files.")

	parser.add_argument("--orfA", metavar = "orf_tax", dest = "orf_alg", \
		help = "Algorithm to use to determine taxonomy of predicted ORFs. " \
		"Use one of 'best_hit', 'lca' or 'lca_star' (default = lca_star), " \
		"which respectively assign the BLAST/HMMer best-hit, the naive LCA as " \
		"in MEGAN and the LCA*. The HMMer best-hit is the LCA* for the " \
		"highest scoring VOG alignment.", \
		default = "lca_star", type = str, \
		choices = ["best_hit", "lca", "lca_star"])

	parser.add_argument("--orfF", metavar = "orf_func", dest = "orf_func", \
		help = "Algorithm to use to determine function of predicted ORFs. " \
		"Use one of 'best_hit' or 'all' (default = 'best_hit'), " \
		"which respectively assign the annotation of the BLAST/HMMer " \
		" best-hit or a list of all annotations of significant BLAST/HMMer " \
		"hits.", \
		default = "best_hit", type = str, \
		choices = ["best_hit", "all"])

	parser.add_argument("-r", metavar = "input_reads", dest = "reads", \
		help = "Path to fastq file containing (QC-ed) reads. If specified, " \
		"the number of mapped bases per contig and predicted ORF will be " \
		"calculated. If using paired-end reads, first merge them using PEAR " \
		"(https://sco.h-its.org/exelixis/web/software/pear)", \
		default = None)

	parser.add_argument("-s", metavar = "align_soft", dest = "soft", \
		help = "Software to align predicted ORFs to database with. Use one " \
		"of 'blast', 'hmmer' or 'both' (default = 'both'). In terms of " \
		"computation time and sensitivity, 'blast' < 'hmmer' < 'both'.", \
		default = "both", type = str, \
		choices = ["blast", "hmmer", "both"])

	user_args = parser.parse_args()
	return(user_args)


def main():

	#Get input
	args = get_user_args()
	args = get_n_threads(args)
	if args.err:
		sys.stderr = open(args.output + ".err", "w")

	#Get OS type to call appropriate bins
	bin_path, os_type = get_bins()

	#Check whether required software is callable
	check_software()

	#Check whether Virtigo is run in initialization mode
	if args.init:
		init_dbs()
		sys.exit(0)

	#Check whether databases are present
	check_dbs()

	#Check input
	args = check_input(args)

	#Predict ORFs and extract the translated sequences from contigs
	print("\nUsing %d threads on a %s based OS." % (args.threads, os_type))
	print("\nPredicting ORFs with MetaGeneAnnotator...")
	args.contigs, tmp = update_contig_ids(args.contigs)
	mga_out = predict_orfs(bin_path, args.contigs, args.output)
	n, orf_count = extract_orfs(mga_out, args.contigs, args.output)
	print("Predicted %d ORFs." % n)


	#Do a BLASTp mapping against pVOGs sequences
	if args.soft == "blast" or args.soft == "both":
		print("\nBLASTp searching %d predicted ORFs against pVOGs " \
			"sequences..." % n)
		blast_search(args.output, args.e_value, args.threads, args.orf_alg)
		m, hmmer_in_tmp = parse_blast_result(args.output)
		hmmer_in = hmmer_in_tmp.name
		print("Mapped %d of %d ORFs using BLASTp." % (m, n))
		if args.soft == "blast":
			del(hmmer_in)


	#Do a HMMer mapping against pVOGs HMMs
	if args.soft == "hmmer":
		hmmer_in = "%s.fasta" % args.output
		m = 0
	if args.soft == "hmmer" or args.soft == "both":
		print("\nHMMer searching %d predicted ORFs against pVOGs " \
			"HMMs..." % (n-m))
		hmmer_search(hmmer_in, args.output, args.e_value, args.threads)
		o = parse_hmmer_result(args.output)
		print("Mapped %d of %d ORFs using HMMer." % (o, n-m))
		del(hmmer_in)
		if args.orf_alg == "best_hit":
			hmmer_best_hit(args.output)


	#Calculate the ORF LCAs based on BLAST results
	print("\nObtaining taxonomic annotations for predicted ORFs and contigs...")
	lca_star(args.output, args.alpha, args.soft, args.orf_alg, \
		args.cont_alg, orf_count)
	print("Done.")


	#Output functional annotation tables
	print("\nObtaining functional annotations for predicted ORFs...")
	get_functions(args.output, args.orf_func, args.soft)
	print("Done.")


	#Map reads to assembled contigs to estimate coverage and abundance
	if args.reads:
		print("\nMapping reads to contigs...")
		map_reads(bin_path, args.contigs, args.reads, args.output, \
			args.threads, mga_out)
		del(mga_out)
		print("Done.")


if __name__ == '__main__':
	main()
