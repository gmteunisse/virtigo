#!/usr/bin/env python

#Import libraries
import sys
import os
import subprocess
import re
from tempfile import NamedTemporaryFile

#Import settings
import modules.settings as settings

#Calls MetaGenomeAnnotator to predict genes and stores predictions in file
def predict_orfs(bin_path, input_file, output_file):

	#Predict ORFs using MGA
	call = "%s/mga %s -m" %  (bin_path, input_file)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	
	#Catch errors
	if error:
		exit_message = "MetaGeneAnnotator error:\n%s" % (error.decode("ascii"))
		sys.stderr.write(exit_message)

	output = output.decode("ascii")
	#Write to file
	write_mga(output, output_file)
	return(output)


#Replaces spaces in fasta headers of contigs with double underscores. 
def update_contig_ids(input_file):

	tmp = NamedTemporaryFile()
	printed = False
	with open(input_file, "r") as file:
		for line in file:
			if line[0] == ">":
				if " " in line and not printed:
					print("Replacing whitespace in contig headers with " \
						"double underscores ('__').")
					printed = True
				#line = line.replace(" ", "__")
				line = re.sub("[ \t]+", '__', line)
			b_line = str.encode(line)
			tmp.write(b_line)
	tmp.flush()
	return(tmp.name, tmp)


#Writes the MGA output to a text file
def write_mga(output, output_file):

	path = "%s.mga" % (output_file)
	file = open(path, "w")
	file.write(output)
	file.close()


def main():
	
	bin_path = bin_path = os.path.join("..", "bin", "osx")
	input_file = os.path.join("..", "data", "contigs", "filtered_VirS1.fna")
	output_file = os.path.join("..", "output", "test")
	output = predict_orfs(bin_path, input_file, output_file)


if __name__ == '__main__':
	main()
