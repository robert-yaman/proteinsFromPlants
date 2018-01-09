from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

class SequenceSet:
	def __init__(self, frequency_dict=None):
		if frequency_dict == None:
			frequency_dict = {}
		# Enforce uniqueness on sequences by using a dict.
		self.frequency_dict = frequency_dict

	def setFrequency(self, sequence, frequency):
		if not type(sequence) == Seq:
			print "ERROR: sequence should be Seq object: " + str(type(sequence))
			raise Exception
		self.frequency_dict[str(sequence)] = frequency

	def frequency(self, sequence):
		if not type(sequence) == Seq:
			print "ERROR: sequence should be Seq object: " + str(type(sequence))
			raise Exception
		# Use str(sequence) for backward compatibility.
		key = str(sequence)
		if not key in self.frequency_dict:
			return 0
		else:
			return self.frequency_dict[key]

	def sequences(self):
		return [Seq(seq, generic_protein) for seq in self.frequency_dict.keys()]