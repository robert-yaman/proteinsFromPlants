import transformation
from cost import Cost
from sequenceSet import SequenceSet

def transformations(params):
	present_sequences = params['sequence_set'].sequences()
	return [Purification(sequence) for sequence in present_sequences]

class Purification(transformation.Transformation):
	def __init__(self, target_sequence):
		self.target_sequence = target_sequence

	def transform(self, input_aas):
		output_aas = SequenceSet()
		for sequence in input_aas.sequences():
			if sequence == self.target_sequence:
				output_aas.setFrequency(sequence, 
					input_aas.frequency(sequence))
				return output_aas
		print "ERROR: Target protein not present in input aas."
		return output_aas

	def cost(self):
		return Cost(50)

	def name(self):
		return "Purify %s" % self.target_sequence
