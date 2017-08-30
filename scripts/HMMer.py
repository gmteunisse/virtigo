#!/usr/bin/python

#Import libraries
import subprocess
import sys
import os
from tempfile import NamedTemporaryFile


#Do a HMMer search vs pVOGs HMMs
def hmmer_search(hmmer_in, output_file, e_value, threads):

	call = 	"hmmscan --incE %f -E %f " \
	"--tblout %s.hmmer " \
	"--cpu %d " \
	"db/pVOGs/hmmer_db/hmmer_db " \
	"%s" % (e_value, e_value, output_file, threads, hmmer_in)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if error:
		exit_message = "HMMer error:\n%s" % (error.decode("ascii"))
		sys.exit(exit_message)

#Count the number of ORFs that has been annotated
def parse_hmmer_result(output_file):

	#Get all gene_ids that have been annotated by HMMer
	hmmer_file = "%s.hmmer" % output_file
	annotated = set()
	with open(hmmer_file, "r") as hmmer:
		for line in hmmer:
			if line[0] != "#":
				data = line.split()
				annotated.add(data[2])

	#Get number of annotated ORFs
	o = len(annotated)
	return(o)


#Reduces HMMer output to only the best hits per gene
def hmmer_best_hit(output_file):

	hmmer_file = "%s.hmmer" % output_file
	tmp = NamedTemporaryFile()
	annotated = set()
	with open(hmmer_file, "r") as hmmer:
		for line in hmmer:
			if line[0] != "#":
				data = line.split()
				if data[2] not in annotated:
					annotated.add(data[2])
				else:
					continue
			b_line = str.encode(line)
			tmp.write(b_line)
	tmp.flush()
	tmp = open(tmp.name).read()
	with open(hmmer_file, "w") as hmmer:
		hmmer.write(tmp)
				

#Main is just used for testing
def main():
	
	output_file = "output/test.fasta"
	e_value = 0.00001
	hmmer_search(output_file, e_value)


if __name__ == '__main__':
	main()
