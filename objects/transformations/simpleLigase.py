import csv

import transformation
from sequenceSet import SequenceSet

def transformations(params):
	with open('ligases.csv', 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		return [SimpleLigase(row[0], row[1], int(row[2])) for row in reader]

class SimpleLigase(transformation.Transformation):
	"""General transformation for basic ligases."""
	def __init__(self, name, code, cost):
		self.name_val = name
		self.code = code
		self.cost_val = cost

	def transform(self, input_aas):
		return self._enzymeBind(self.code, input_aas)

	def cost(self):
		return self.cost_val

	def name(self):
		return self.name_val

	def _enzymeBind(self, code, input_aas):
		positive_list, negative_list = self.parse_code(code)
		if (not len(positive_list) % 2 == 0 
			and len(positive_list) == len(negative_list)):
			print "ERROR: condition lists have incorrect lengths: %d, %d." %(
				positive_list, negative_list)

		# Conditions on the FRONT of the protein i.e. N-terminus (end of the 
		# condition list).
		positive_front = positive_list[len(positive_list) / 2 :]
		negative_front = negative_list[len(negative_list) / 2 :]
		# Conditions on the BACK of the protein i.e. C-terminus (beginning of 
		# the condition list).
		positive_back = positive_list[0 : len(positive_list) / 2]
		negative_back = negative_list[0 : len(negative_list) / 2]

		def matches_positive_front(sequence, pattern):
			matches = True
			for i, condition in enumerate(pattern):
				if len(sequence) <= i:
					if condition:
						matches = False
						break
					else:
						continue
				if condition and not sequence[i] in condition:
					matches = False
					break
			return matches

		def matches_positive_back(sequence, pattern):
			matches = True
			for i, condition in enumerate(pattern):
				seq_index = -len(positive_back) + i
				if len(sequence) < abs(seq_index):
					if condition:
						matches = False
						break
					else:
						continue
				if condition and not sequence[seq_index] in condition:
					matches = False
					break
			return matches

		def avoids_negative_front(sequence, pattern):
			avoids = True
			for i, condition in enumerate(pattern):
				if len(sequence) <= i:
					continue
				if condition and sequence[i] in condition:
					avoids = False
					break
			return avoids

		def avoids_negative_back(sequence, pattern):
			avoids = True
			for i, condition in enumerate(pattern):
				seq_index = -len(positive_back) + i
				if len(sequence) < abs(seq_index):
					continue
				if condition and sequence[seq_index] in condition:
					avoids = False
					break
			return avoids

		sequences = input_aas.sequences()
		output_sequences = SequenceSet()
		front_sequences = []
		back_sequences = []
		for sequence in sequences:
			matched = False
			if (matches_positive_front(sequence, positive_front) and 
				avoids_negative_front(sequence, negative_front)):
				front_sequences.append(sequence)
				matched = True
			if (matches_positive_back(sequence, positive_back) and
				avoids_negative_back(sequence, negative_back)):
				back_sequences.append(sequence)
				matched = True
			if not matched:
				output_sequences.setFrequency(sequence, 
					input_aas.frequency(sequence))

		# If there is nothing to bind, add proteins back in.
		if len(front_sequences) == 0:
			for back_sequence in back_sequences:
				output_sequences.setFrequency(back_sequence, 
					input_aas.frequency(back_sequence))
		if len(back_sequences) == 0:
			for front_sequence in front_sequences:
				output_sequences.setFrequency(front_sequence, 
					input_aas.frequency(front_sequence))

		# We assume that there will be an equal distribution among all the 
		# pairs of proteins that can bind with each other. For now, the number
		# of merged sequences is the lower of the two subtrates. Will create a
		# more nuanced picture later.
		for front_sequence in front_sequences:
			for back_sequence in back_sequences:
				combined_sequence = back_sequence + front_sequence
				combined_frequency = min(input_aas.frequency(
					front_sequence), input_aas.frequency(back_sequence))
				output_sequences.setFrequency(combined_sequence, 
					combined_frequency)

		return output_sequences