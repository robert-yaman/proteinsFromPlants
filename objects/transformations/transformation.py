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

	def _enzymeCleave(self, code, sequence_set):
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
		# On each sequence, indicates indices that should be cleaved at the
		# C-terminal.
		cleaveage_sites = {}
		sequences = sequence_set.sequences()
		for sequence in sequences:
			for i, cleaveage_site in enumerate(sequence):
				# Iterating over amino acids which we consider cleaving on the
				# C-terminal. Don't consider the last AA since we can't cleave 
				# there.
				if i == len(sequences) - 1:
					continue
				cleave = True
				for p in range(len(positive_list)):
					if cleave == False:
						continue
					p_relation_to_i = p - 3
					p_index_on_sequence = i + p_relation_to_i
					if (p_index_on_sequence < 0 or 
						p_index_on_sequence >= len(sequence)):
						if positive_list[p]:
							# We require an amino acid to be here, but we are
							# outisde of the sequence bounds. Don't cleave.
							cleave = False
							continue

					if (positive_list[p] and 
						not sequence[p_index_on_sequence] in positive_list[p]):
						cleave = False

					if (negative_list[p] and 
						sequence[p_index_on_sequence] in negative_list[p]):
						cleave = False
				if cleave:
					if sequence in cleaveage_sites:
						cleaveage_sites[sequence].append(i)
					else:
						cleaveage_sites[sequence] = [i]

		# Perform actual cleaving.
		if len(cleaveage_sites) == 0:
			return sequence_set
		output_set = SequenceSet()
		for sequence in sequences:
			if not sequence in cleaveage_sites:
				output_set.addFrequency(sequence, 
					sequence_set.frequency(sequence))
			else:
				cleavages = cleaveage_sites[sequence]
				cleavages.sort()
				# Add 1 for more convenient indexing
				cleavages = [c + 1 for c in cleavages]
				cleavages_with_padding = [0] + cleavages + [len(sequence)]
				for new_seq_i in range(len(cleavages_with_padding) - 1):
					new_sequence = sequence[cleavages_with_padding[new_seq_i] : cleavages_with_padding[new_seq_i + 1]]
					output_set.addFrequency(new_sequence, 
						sequence_set.frequency(sequence))

		return output_set

	def _sequence(self, sequence):
		return Seq(sequence, generic_protein)
