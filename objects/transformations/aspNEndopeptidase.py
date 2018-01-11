import transformation
from cost import Cost
from sequenceSet import SequenceSet

def transformations(params):
	return [AspNEndopeptidase()]

class AspNEndopeptidase(transformation.Transformation):
	"""Source: http://enzyme.expasy.org/EC/3.4.24.33
	TODO: check this one - there is some discrepancy in expasy.
	"""
	def transform(self, input_aas):
		return self._enzymeCleave('----DQ---', input_aas)

	def cost(self):
		return 1

	def name(self):
		return "Asp-N Endopeptidase"
