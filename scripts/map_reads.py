#!/usr/bin/env python

#Import libraries
import sys
import os
import subprocess
from tempfile import mkdtemp
from shutil import rmtree
from scripts.extract_orfs import ORF, parse_mga_out
from collections import OrderedDict

#Import settings
import scripts.settings as settings

class Contig_coverage:

	def __init__(self, contig_name):

		self.id = contig_name
		self.length = int()				#Contig length
		self.mapped = int() 			#Number of mapped bases
		self.coverage = OrderedDict()	#Coverage (item) per segment (key)


	def add_coverage(self, end, coverage):

		self.coverage[end] = coverage


	#Calculate the number of mapped bases for this contig
	def calc_map(self):

		i = 0
		for j in self.coverage.keys():
			length = j - i
			self.mapped += length * self.coverage[j]
			i = j
		self.length = i


	#Calculate the number of bases mapped to a segment between start and stop
	def calc_segment_map(self, start, stop):

		i = start - 1	#Correct for offset
		stop = stop - 1
		mapped = 0
		first = True

		for j in self.coverage.keys():
			if i >= j:	#Start of ORF not in segment
				continue
			elif i < j:	#ORF in segment
				if stop >= j: 	#ORF spans full segment
					length = j - i
					mapped += length * self.coverage[j]
					i = j
				else: 			#ORF spans until ORF stop
					length = stop - i
					mapped += length * self.coverage[j]
					break
		return(int(mapped))

#Calls Bowtie V2.1 to create an index of contig sequences
def bowtie_build(bin_path, contigs):

	#Create a temporary directory to store index files
	cwd = os.getcwd()
	tmp_dir = mkdtemp(dir = cwd, prefix = "bowtie_tmp_")

	#Change dir, because bowtie can't be run from a different dir.
	os.chdir(os.path.join(bin_path, "bowtie"))

	#Create Bowtie index files for the contigs file
	call = "bowtie2-build %s %s/tmp" % \
		(contigs, tmp_dir)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()

	#Catch errors
	if error:
		warning_message = "Bowtie error/output:\n%s" % (error.decode("ascii"))
		sys.stderr.write(warning_message)
	return(tmp_dir, cwd)


#Calls Bowtie V2.1 to map reads against contig index files
def bowtie_map(reads, tmp_dir, n_threads):
	
	call = 	"bowtie2 " \
			"-p %d " \
			"-x %s/tmp " \
			"-U %s " \
			"-S %s/tmp.sam" % \
		(n_threads, tmp_dir, reads, tmp_dir)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()

	#Catch errors
	if error:
		warning_message = "Bowtie error/output:\n%s" % (error.decode("ascii"))
		sys.stderr.write(warning_message)


#Use BEDtools v2.26 to obtain the coverage of each contig
def calc_cov(tmp_dir):

	call = "genomeCoverageBed -ibam %s/tmp.bam -bga" % (tmp_dir)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()

	#Catch errors
	if error:
		error_message = "BEDtools error:\n%s" % error.decode("ascii")
		sys.stderr.write(error_message)

	output = output.decode("ascii")

	#Store all coverage data in useable objects
	coverage = dict()
	for line in output.splitlines():
		data = line.split("\t")
		contig = data[0]
		end = int(data[2])
		cov = int(data[3])
		if contig not in coverage.keys():
			coverage[contig] = Contig_coverage(contig)
		coverage[contig].add_coverage(end, cov)

	return(coverage)


#Calculate the number of mapped based per contig and per gene
def calc_bases_mapped(coverage, mga_out):
	
	#Calculate the number of bases mapped per contig and per predicted ORF
	orfs = parse_mga_out(mga_out)
	for contig in orfs:
		#Contig mapping
		coverage[contig].calc_map()
		for orf in orfs[contig]:
			#ORF mapping
			orf.aa_length = orf.stop - orf.start
			orf.mapped = coverage[contig].calc_segment_map(orf.start, \
																orf.stop)
	return(coverage, orfs)

#Maps reads to contigs, removes PCR duplicates and calculates the number of
#mapped bases per contig and predicted ORF.
def map_reads(bin_path, contigs, reads, output_file, n_threads, mga_out):

	#Align reads to contigs with bowtie
	tmp_dir, cwd = bowtie_build(bin_path, contigs) #Also switches dir
	bowtie_map(reads, tmp_dir, n_threads)
	os.chdir(cwd) #Move back to original dir

	#Create SAM file
	samtools_build(contigs, tmp_dir)
	#remove_duplicates(bin_path, tmp_dir, output_file)

	#Calculate coverage and mapping per contig and predicted ORF
	coverage = calc_cov(tmp_dir)
	coverage, orfs = calc_bases_mapped(coverage, mga_out)

	#Write calculated mapping to file
	write_mapping(output_file, coverage, orfs)

	#Remove tmp_dir and its contents
	tmp_dir = tmp_dir.split("/")[-1]
	rmtree(tmp_dir)


#This function is not yet properly implemented
'''#Remove PCR duplicates using picard v2.10 and samtools v1.5
def remove_duplicates(bin_path, tmp_dir, output_file):

	call = "bash scripts/picard_remove_duplicates.sh %s %s/tmp %s" % \
																(bin_path,
																tmp_dir, \
																output_file)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if error:
		error_message = "Samtools / picard error:\n%s" % error.decode("ascii")
		sys.stderr.write(error_message)'''


#Uses samtools v1.5 to create a sorted and indexed BAM-file
def samtools_build(contigs, tmp_dir):

	call = "bash scripts/samtools_build_bam.sh %s %s/tmp" % (contigs, tmp_dir)
	process = subprocess.Popen(call.split(), stdout = subprocess.PIPE, \
		stderr = subprocess.PIPE)
	output, error = process.communicate()
	if error:
		error_message = "Samtools error:\n%s" % error.decode("ascii")
		sys.stderr.write(error_message)

#Write mapped bases per contig and predicted ORFs to files
def write_mapping(output_file, coverage, orfs):

	#Open files
	contig_file = open(output_file + ".con.tbl", "w")
	orf_file = open(output_file + ".orf.tbl", "w")

	#Add headers
	contig_file.write("#contig_id\tcontig_length\tmapped_bases\n")
	orf_file.write("#contig_id\tORF_id\tORF_length\tmapped_bases\n")

	#Write mapped bases
	for contig in orfs:
		contig_file.write("%s\t%d\t%d\n" % (coverage[contig].id, \
											coverage[contig].length, \
											coverage[contig].mapped))
		for orf in orfs[contig]:
			orf_file.write("%s\t%s\t%d\t%d\n" % (coverage[contig].id, \
											orf.gene_id, \
											orf.aa_length, \
											orf.mapped))
	contig_file.close()
	orf_file.close()


#The main function is only used for testing.
def main():

	n_threads = 8
	bin_path = os.path.join("..", "bin", "osx")
	bin_path = os.path.abspath(bin_path)
	contigs = os.path.join("..", "data", "contigs", "reference", "lambda_virus.fa.tmp")
	contigs = os.path.abspath(contigs)
	reads = os.path.join("..", "data", "reads", "reads", "reads_1.fq")
	reads = os.path.abspath(reads)
	output_file = os.path.join("..", "output", "lambda_virus")
	output_file = os.path.abspath(output_file)
	mga_out = open(os.path.join("..", "output", "lambda_virus.mga")).read()

	tmp_dir, cwd = bowtie_build(bin_path, contigs)
	bowtie_map(reads, tmp_dir, n_threads)
	#Move back to original dir
	os.chdir(cwd)
	samtools_build(contigs, tmp_dir)
	remove_duplicates(bin_path, tmp_dir, output_file)
	coverage = calc_cov(tmp_dir)
	coverage, orfs = calc_bases_mapped(coverage, mga_out)
	write_mapping(output_file, coverage, orfs)
	tmp_dir = tmp_dir.split("/")[-1]
	#Remove tmp_dir and its contents
	rmtree(tmp_dir)


if __name__ == '__main__':
	main()
