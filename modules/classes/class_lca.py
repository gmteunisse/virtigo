#!usr/bin/python

import os
import re
import sys
from math import log
try:
	from scipy.stats import chisquare
	scipy = True
except ImportError:
	sys.stderr.write("Error: no SciPy installation found. Cannot compute " \
	"p-value for LCA* assignment, will report 0 instead.\n")
	scipy = False

#Import classes
from modules.classes.class_node import Node


class LCA:

	def __init__(self, id, tree, alpha):

		self.id = id 			#ORF or contig id
		self.taxids = list()	#Taxids assigned
		self.to_visit = list()	#Taxids + all parent taxids
		self.visited = list() 	#Nodes to use in LCA* calculation
		self.total = int()		#Total number of taxids assigned
		self.tree = tree 		#NCBI tax tree for viruses

		self.alpha = alpha 		#Majority treshold

		self.H = dict()			#Entropy (H) per taxid
		self.S = dict()			#S = sum n per taxid
		self.L = dict()			#L = sum n log n per taxid

		self.lca_taxid = str()	#Stores taxid of LCA*

	#Calculate the LCA* based on available taxids.
	def calc_lca(self, total = False):

		#Add counts and "color" the visited nodes
		if not total:
			self.total = len(self.taxids)
		else:
			self.total = total
		for id in self.taxids:
			try:
				self.to_visit += self.tree[id].add_count(1)
			#This means that the node_id does not exist. This should not
			#occur, but sometimes it happens inexplicably.
			except KeyError:
				self.to_visit += self.tree["10239"].add_count(1)
				sys.stderr.write("Error: non-existent node with TaxID %s" % id)
		self.to_visit = set(self.to_visit)	#Reduce to union

		candidate = ["10239", 10000000.00] 	#Root with infinite entropy
		stack = ["10239"]					#List of taxids to consider
		while len(stack) > 0:

			id = stack.pop()

			if id in self.visited:

				c = list()	#Children with counts in their downstream lineage
				for c_id in self.tree[id].children.keys():
					if self.tree[c_id].visited:
						c.append(c_id)
				
				self.H[id] = 0 	#Entropy = 0.

				count = self.tree[id].count
				if count > 0: #Prevent logarithmic errors
					self.S[id] = float(count)
					self.L[id] = float(count) * log(float(count))
				else:
					self.S[id] = 0
					self.L[id] = 0

				for c_id in c:
					self.S[id] += self.S[c_id]
					self.L[id] += self.L[c_id]

				try:
					self.H[id] = -(self.L[id]/self.S[id] - log(self.S[id]))
				except:
					break
				
				if self.S[id] > self.total * self.alpha:
					if candidate[1] > self.H[id]:
						candidate = [id, self.H[id]]
			else:
				self.visited.append(id)
				stack.append(id)
				c = list()
				for c_id in self.tree[id].children.keys():
					if self.tree[c_id].visited:
						c.append(c_id)
						stack.append(c_id)
		self.lca_taxid = candidate[0]
		return(self.lca_taxid)


	#Calculate a chi-square test statistic and p-value for a collapsed
	#distribution of counts.
	def calc_p(self):

		if scipy:
			#Collapse children of LCA* and remove from distribution
			to_remove = self.tree[self.lca_taxid].clean_children()
			for id in to_remove:
				try:
					self.taxids.remove(id)
				except ValueError:
					pass
			self.taxids.append(self.lca_taxid)

			#Get distribution of counts
			counts = list()
			for id in self.taxids:
				counts.append(self.tree[id].count)
			unknown = self.total - sum(counts)	#Add unannotated ORFs
			if unknown > 0:
				counts.append(unknown)

			#Perform chi-square
			if len(counts) > 1:
				p_value = chisquare(counts)[1]
			else:
				p_value = float(0) 	#LCA* == LCA, thus no p-value can be calculated
			return(p_value)
		else:
			return(float(0))


	#Revert the NCBI taxonomic tree back to its clean state, i.e. all counts
	#back to 0.
	def clean_tree(self):

		for id in self.to_visit:
			self.tree[id].count = 0
			self.tree[id].count_collapsed = 0
			self.tree[id].visited = False


	#Print all nodes in the tree that have counts assigned to them.
	def print_tree(self, out_file):

		with open(out_file, "a") as out_file:
			out_file.write(">Contig:\t%s\nNumber of predicted ORFs:\t%s\n\n" % (self.id, self.total))
			header = "\t" + "%s\t" * 9 + "%s\n\n"
			out_file.write(header % ("s_king", "gen_typ", "order", "family", \
			 				"subfam", "genus", "sub_gen", "species", "sub_sp", "no rank"))
			to_write = str()
			tree = self.tree["10239"].print_tree(depth = 10, lca = True) #10239 is root
			if tree:
				for subtree in tree:
					for line in subtree:
						to_write += line
			to_write += "\n"
			out_file.write(to_write)

