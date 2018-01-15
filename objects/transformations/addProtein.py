import csv

import transformation
from sequenceSet import SequenceSet

def transformations(params):
	# |start_seqs| is a list of lists (starter_sequence, starter_cost).
	start_seqs = params['start_seqs']
	return [AddProtein(starter[0], starter[1]) for starter in start_seqs]

class AddProtein(transformation.Transformation):
	def __init__(self, protein, cost):
		self.protein = protein
		self.cost_val = cost

	def transform(self, input_aas):
		output_aas = SequenceSet()
		for sequence in input_aas.sequences():
			output_aas.setFrequency(sequence, input_aas.frequency(sequence))
		output_aas.setFrequency(self.protein, 
			output_aas.frequency(self.protein) + 1)
		return output_aas

	def cost(self):
		return self.cost_val

	def name(self):
		prefix = self.protein[0:4]
		if len(self.protein) > 4:
			prefix += "..."
		return "Add %s" % prefix
