from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

class SequenceSet:
	def __init__(self, frequency_dict=None):
		if frequency_dict == None:
			frequency_dict = {}
		# Enforce uniqueness on sequences by using a dict.
		self.frequency_dict = frequency_dict

	def __len__(self):
		return len(self.frequency_dict)

	def __hash__(self):
		sorted_keys = self.frequency_dict.keys()
		sorted_keys.sort()
		hash_string = ""
		for key in sorted_keys:
			hash_string += "%s%d" % (key, self.frequency_dict[key])
		return hash(hash_string)

	def __eq__(self, other):
		# Use internal struct for efficiency here. It's fine since it's an 
		# internal method.
		return self.frequency_dict == other.frequency_dict

	def setFrequency(self, sequence, frequency):
		if not type(sequence) == Seq:
			print "ERROR: sequence should be Seq object: %s" %(
				str(type(sequence)))
			raise Exception
		self.frequency_dict[str(sequence)] = frequency
		# Don't store 0s for equality comparisons.
		if int(frequency) == 0:
			del self.frequency_dict[str(sequence)]

	def frequency(self, sequence):
		if not type(sequence) == Seq:
			print "ERROR: sequence should be Seq object: %s" %(
				str(type(sequence)))
			raise Exception
		# Use str(sequence) for backward compatibility.
		key = str(sequence)
		if not key in self.frequency_dict:
			return 0
		else:
			return self.frequency_dict[key]

	def addFrequency(self, sequence, frequency):
		self.setFrequency(sequence, self.frequency(sequence) + frequency)

	def sequences(self):
		return [Seq(seq, 
			generic_protein) for seq in self.frequency_dict.keys()]

	def contains(self, seq):
		return self.frequency(seq) > 0

	def isExactly(self, seq):
		return len(self.sequences()) == 1 and self.contains(seq)
