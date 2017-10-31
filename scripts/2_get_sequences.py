#!usr/bin/python

import os

class Fasta_Sequence:

	def __init__(self):

		self.header = str()
		self.sequence = ""


def create_sequence_db():

	base_path = os.getcwd()
	sequence_db_path = os.path.join(base_path, "db", "pVOGs", "sequences")
	try:
		os.mkdir(sequence_db_path)
	except OSError:
		print("Directory already exists, files may be overwritten.")
	return(sequence_db_path)


def get_file_names():

	base_path = os.getcwd()
	alignment_db_path = os.path.join(base_path, "db", "pVOGs", "alignments")
	file_names = os.listdir(alignment_db_path)
	return(file_names, alignment_db_path)


def remove_gaps(files, aln_path):

	vogs = dict()

	for file in files:

		file_path = os.path.join(aln_path, file)
		data = open(file_path, "r")
		lines = data.readlines()
		seqs = list()
		i = -1

		for line in lines:
			if line[0] == ">":
				i = i + 1
				seqs.append(Fasta_Sequence())
				seqs[i].header = line.rstrip()
			else:
				line = line.replace("-", "") #Remove gap characters ("-")
				seqs[i].sequence = "".join((seqs[i].sequence, line.rstrip()))

		vog_id = file.replace(".aln", "")
		vogs[vog_id] = seqs
	return(vogs)


def write_seqs(vogs, seq_path):

	for vog_id in vogs:

		file_path = os.path.join(seq_path, ".".join((vog_id, "fasta")))
		file = open(file_path, "w")
		for seq in vogs[vog_id]:
			file.write("%s\n%s\n" % (seq.header, seq.sequence))
		file.close()


def main():

	files, aln_path = get_file_names()
	seq_path = create_sequence_db()
	vogs = remove_gaps(files, aln_path)
	write_seqs(vogs, seq_path)


if __name__ == '__main__':
	main()