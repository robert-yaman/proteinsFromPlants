import transformation
from cost import Cost
from sequenceSet import SequenceSet

def transformations(params):
	start_proteins = params['start_seqs']
	return [AddProtein(protein) for protein in start_proteins]

class AddProtein(transformation.Transformation):
	def __init__(self, protein):
		self.protein = protein

	def transform(self, input_aas):
		output_aas = SequenceSet()
		for sequence in input_aas.sequences():
			output_aas.setFrequency(sequence, input_aas.frequency(sequence))
		output_aas.setFrequency(self.protein, 
			output_aas.frequency(self.protein) + 1)
		return output_aas

	def cost(self):
		return 25

	def name(self):
		prefix = self.protein[0:4]
		if len(self.protein) > 4:
			prefix += "..."
		return "Add %s" % prefix
