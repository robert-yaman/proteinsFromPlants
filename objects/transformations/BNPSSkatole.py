import transformation
from cost import Cost
from sequenceSet import SequenceSet

def transformations(params):
	return [BNPSSkatole()]

class BNPSSkatole(transformation.Transformation):
	"""Source: http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions"""
	def transform(self, input_aas):
		return self._enzymeCleave('---W----', input_aas)

	def cost(self):
		return 12

	def name(self):
		return "BNPS-Skatole"
