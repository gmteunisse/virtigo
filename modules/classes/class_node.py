#!/usr/bin/python
import re
import sys
from collections import OrderedDict

class Node:

	def __init__(self, line):

		self.line = line
		self.level = str()
		self.name = str()
		self.id = str()
		self.parent_id = str()
		self.parent = None
		self.children = dict()

		self.tax_levels = ["superkingdom", "genome type", "order", "family", \
		 "subfamily", "genus", "subgenus", "species", "subspecies", "no rank"]
		self.index = int()

		self.count = 0
		self.count_collapsed = 0
		self.visited = False


	#Add a count to a node and color its parents. Also add a collapsed tree
	#count: the sum of the counts of a node and all its children.
	def add_count(self, count):

		self.count_collapsed += 1
		self.count += count
		self.visited = True
		if self.level != "superkingdom":
			visited = self.parent.add_count(0) #'True' marks parents as visisted
		else:
			visited = list()
		visited.append(self.id)
		return(visited)


	#Calculate the index of a node based on it's taxonomic level
	def calc_index(self):

		self.index = self.tax_levels.index(self.level)
	

	#Find the node that corresponds to a taxid. Once found, obtain 
	#annotation of parent nodes and return as a dict.
	def find_node(self, tax, level):
		tax_index = self.tax_levels.index(level)
		if tax_index == self.index + 1:
			if tax in self.children.keys():
				ancestry = self.return_parents()
				return(ancestry)
		elif tax_index > self.index + 1:
			for key in self.children.keys():
				ancestry = self.children[key].find_node(tax, level)
				if ancestry:
					return(ancestry)
		else:
			exit_message = "Did not find node. This should not happen."
			sys.exit(exit_message)


	#Returns levels of children
	def get_child_levels(self):

		levels = list()
		for key in self.children.keys():
			levels.append(self.children[key].level)
		return(levels)


	#Parse a line describing an NCBI taxon node
	def parse(self):

		data = self.line.split("\t")
		self.id = data[0].strip()
		self.parent_id = data[1].strip()
		self.level = data[2].strip()
		if self.level not in self.tax_levels:
			self.level = "no rank"
		self.name = data[3].strip()		


	#Print the current node name and recursively print child node names in a 
	#tree-like structure.
	def print_tree(self, depth, lca = False, sisters = False):

		if self.index < depth:
			if self.visited:
				tree = list()
				active_children = 0
				for key in self.children.keys():
					if self.children[key].visited:
						active_children += 1
				if self.parent:
					length = self.index - self.parent.index
				else:
					length = 1
				for key in self.children.keys():
					if self.children[key].visited:
						active_children -= 1
						if active_children > 0:
							pass_sisters = True
						else:
							pass_sisters = False
						subtree = self.children[key].print_tree(depth, lca, pass_sisters)
						if not subtree:
							continue
						for line in subtree:
							if sisters:
								line = "|" + " " * 7 + " " * 8 * (length - 1)+ line
							else:
								line = " " * 8 * length + line
							tree.append(line)
				branch = "+" + "-" * 7 + ("-" * 8 * (length - 1))
				line = "%s%s (%d; %d)\n" % (branch, self.name, \
										self.count, self.count_collapsed)
				tree.insert(0, line)
				return(tree)


	'''#Print the current node name and recursively print child node names in a 
	#tree-like structure.
	def print_tree(self, depth, lca = False, prev_ind = 0):
		to_write = str()
		if not lca:
			if self.index < depth:
				indent = "\t" * self.index
				print("%s%s" % (indent, self.name))
				for key in self.children.keys():
					self.children[key].print_tree(depth)
		else:
			if self.index < depth:
				if self.visited:
					indent = "\t" * prev_ind * 2
					indent += "+" + "-"* 7 * (self.index - prev_ind)
					to_write = "%s%s (%d; %d)\n" % (indent, self.name, \
											self.count, self.count_collapsed)
					for key in self.children.keys():
						to_write += self.children[key].print_tree(depth, lca, self.index)
					return(to_write)

		return(to_write)'''


	#Recursively reset the visited status of child nodes to their unvisited
	#state
	def clean_children(self):

		self.count = self.count_collapsed #Children are set to 0, so collapse
		children = list()
		if len(self.children) != 0:
			for id in self.children.keys():
				if self.children[id].visited:
					self.children[id].count = 0
					self.children[id].count_collapsed = 0
					self.children[id].visited = False
					children.append(id)
					children += self.children[id].clean_children()
		return(children)



	#Return annotation of parent nodes, unless at root.
	def return_parents(self):

		if self.level == "superkingdom":
			ancestry = OrderedDict()
		else:
			ancestry = self.parent.return_parents()
		ancestry[self.level] = self.name
		return(ancestry)

