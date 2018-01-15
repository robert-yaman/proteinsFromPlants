import transformation
from sequenceSet import SequenceSet

def transformations(params):
	return [ArgCProteinase()]

class ArgCProteinase(transformation.Transformation):
	"""Source: http://enzyme.expasy.org/EC/3.4.21.35"""
	def transform(self, input_aas):
		return self._enzymeCleave('---R----', input_aas)

	def cost(self):
		return 1

	def name(self):
		return "Arg-C Proteinase"
