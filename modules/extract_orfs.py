#!/usr/bin/env python

#Import libraries
import sys
import os

#Import settings
import modules.settings as settings


#Class that stores data for an MGA predicted ORF.
class ORF:

	def __init__(self, contig_id):

		self.contig = contig_id
		self.gene_id = str()
		self.start = int()
		self.stop = int()
		self.strand = str()
		self.frame = int()
		self.complete = str()
		self.org = str()
		self.nt_seq = str()
		self.aa_seq = str()
		self.length = int()
		self.mapped = int()


	#Get the complementary sequence in case of anti-sense strand prediction
	def get_complement(self, nt_seq):

		compl_seq = str()
		for nt in nt_seq:				#Get complement
			if nt == "A":
				compl_seq += "T"
			elif nt == "T":
				compl_seq += "A"
			elif nt == "G":
				compl_seq += "C"
			elif nt == "C":
				compl_seq += "G"
		compl_seq = compl_seq[::-1] 	#Reverse the sequence
		return(compl_seq)


	#Gets the nucleotide sequence for this predicted gene from a contig seq.
	def get_nt_seq(self, seq):

		if self.strand == "+":
			self.nt_seq = seq[self.start-1:self.stop-1].upper()
		elif self.strand == "-":
			self.nt_seq = seq[self.start:self.stop].upper()
			self.nt_seq = self.get_complement(self.nt_seq)


	#Parses an MGA output line
	def parse(self, line):

		data = line.strip().split("\t")
		if len(data) > 1:
			self.gene_id = data[0]
			self.start = int(data[1])
			self.stop = int(data[2])
			self.strand = data[3]
			self.frame = int(data[4])
			self.complete = data[5]
			self.org = data[7]


	#Translate an nt sequence to an amino acid sequence
	def transl_seq(self, codon_table):

		for i in range(self.frame, len(self.nt_seq), 3):
			codon = self.nt_seq[i:i+3]
			if len(codon) == 3:
				aa = codon_table[codon]
				if aa == "*":			#Stop codon, so stop translating
					break
				else:
					self.aa_seq += aa
			else:
				break


#Count the total number of predicted orfs 
def count_orfs(orfs):

	n = 0
	orf_count = dict()
	for contig in orfs:
		n += len(orfs[contig])
		orf_count[contig] = len(orfs[contig])
	return(n, orf_count)


#Extracts the predicted ORF sequences, count them and translates them.
def extract_orfs(mga_out, input_file, output_file):

	orfs = parse_mga_out(mga_out)
	n_genes, orf_count = count_orfs(orfs)
	translate_seqs(orfs, input_file, output_file)
	return(n_genes, orf_count)


#Extract the nucleotide sequence of each predicted ORF
def extract_seqs(contig_id, seq, orfs):
	
	codon_table = read_codon_tab()
	for gene in orfs[contig_id]:
		gene.get_nt_seq(seq)
		gene.transl_seq(codon_table)


#Store the MGA output into gene objects per contig. Each MGA annotation
#Consists of three lines starting with "#", followed by lines with predicted
#ORFs. 
def parse_mga_out(mga_out):

	orfs = dict()
	count = 0
	for line in mga_out.splitlines():

		if line[0] == "#":		#Contig header line
			if count == 0:
				id = line[1:].strip()
				orfs[id] = list()
			count += 1
			if count == 3:		#Ready for the next contig
				count = 0
		elif line[0] != "#":	#ORF line
			orf = ORF(id)
			orf.parse(line)
			orfs[id].append(orf)
	return(orfs)


#Reads a codon translation table.
def read_codon_tab():

	codon_table = dict()
	path = os.path.join(settings.codon_tbl_path)
	for line in open(path, "r"):
		data = line.strip().split("\t")
		codon_table[data[0]] = data[1]
	return(codon_table)


#Extract the corresponding nucleotide sequences from the input multi FASTA,
#translate them and store in file.
def translate_seqs(orfs, input_file, output_file):

	#Open a file for writing
	file_path = "%s.fasta" % (output_file)
	file = open(file_path, "w")

	#Read the contig file, extract the ORF sequences and write them to file
	with open(input_file) as contig_file:
		first_line = True
		for line in contig_file:
			if line[0] == ">":							#Header line
				if not first_line:						#Extract and write seq
					extract_seqs(contig_id, seq, orfs)
					write_seqs(contig_id, orfs, file)
					orfs.pop(contig_id)					#Remove to save mem
				first_line = False
				contig_id = line[1:].strip()
				seq = str()								#New empty sequence
			else:										#Sequence line
				seq += line.strip()
	#Extract and write away the last sequences
	extract_seqs(contig_id, seq, orfs)
	write_seqs(contig_id, orfs, file)
	file.close()


#Writes the translated ORFs to file.
def write_seqs(contig_id, orfs, file):

	for gene in orfs[contig_id]:
		file.write(">%s|%s\n" % (gene.contig, gene.gene_id))
		for i in range(0, len(gene.aa_seq), 60):
			file.write("%s\n" % (gene.aa_seq[i:i+60]))


def main():
	
	input_file = os.path.join("..", "data", "contigs", "test.fna")
	output_file = output_file = os.path.join("..", "output", "test")
	mga_out = open(os.path.join("..", "output", "test.mga"))
	extract_orfs(mga_out, input_file, output_file)


if __name__ == '__main__':
	main()
