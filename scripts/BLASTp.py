#!/usr/bin/python

#Import libraries
import subprocess
import sys
from tempfile import NamedTemporaryFile


#Do a BLASTp search vs pVOGs sequences
def blast_search(output_file, e_value, threads, method):

	blast_opt = ""
	if method == "best_hit":
		blast_opt = "-max_target_seqs 1"

	call = 	"blastp -query %s.fasta " \
			"-outfmt 6 " \
			"-evalue %f " \
			"-num_threads %d " \
			"-db db/pVOGs/blast_db/pVOGs_db " \
			"%s " \
			"-out %s.blast "  %  (output_file, e_value, threads, blast_opt,
									output_file)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if output:
		print("BLASTp:\n%s" % (output.decode("ascii")))
	if error:
		exit_message = "BLASTp error:\n%s" % (error.decode("ascii"))
		sys.exit(exit_message)


#Create a temporary file with sequences for HMMer searching
def parse_blast_result(output_file):

	#Get all gene_ids that have been annotated by BLAST
	blast_file = "%s.blast" % output_file
	annotated = set()
	with open(blast_file, "r") as blast:
		for line in blast:
			data = line.split("\t")
			annotated.add(data[0])

	#Get number of annotated ORFs
	m = len(annotated)
	
	#Write all remaining fastas to a temporary file for hmmer
	fasta_file = "%s.fasta" % output_file
	hmmer_in = NamedTemporaryFile()
	with open(fasta_file, "r") as fasta:
		for line in fasta:
			if line[0] == ">":
				gene_id = line[1:].strip()
				if gene_id in annotated:
					write = False
				else:
					write = True
			if write:
				b_line = str.encode(line)
				hmmer_in.write(b_line)
	
	#Make temp file readable and return
	hmmer_in.flush()
	return(m, hmmer_in)
		

def main():
	
	output_file = "output/test"
	e_value = 0.001
	blast_search(output_file, e_value)
	parse_blast_result(output_file)


if __name__ == '__main__':
	main()
