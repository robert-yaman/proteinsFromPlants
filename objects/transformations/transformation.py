import abc

from Bio.Alphabet import generic_protein
from Bio.Seq import Seq

from sequenceSet import SequenceSet

class Transformation:
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def transform(self, input_aas):
		pass

	@abc.abstractmethod
	def cost(self):
		pass

	@abc.abstractmethod
	def name(self):
		pass

	def parse_code(self, code):
		"""Convenience method that parses a peptide filter."""
		code_list = code.split("-")
		# Disjunctive positive conditions for area around possible cleavage 
		# site.
		positive_list = []
		# Conjunctive negative conditions.
		negative_list = []
		for item in code_list:
			split_item = item.split("^")
			positive_list.append(split_item[0])
			if len(split_item) == 1:
				negative_list.append("")
			else:
				negative_list.append(split_item[1])

		return positive_list, negative_list


	def _sequence(self, sequence):
		return Seq(sequence, generic_protein)
