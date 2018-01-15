import transformation
from sequenceSet import SequenceSet

def transformations(params):
	return [Caspase1()]

class Caspase1(transformation.Transformation):
	"""Source: http://www.uniprot.org/uniprot/P29466"""
	def transform(self, input_aas):
		return self._enzymeCleave('FWYL--HAT-D-^PEDQKR---', input_aas)

	def cost(self):
		return 1

	def name(self):
		return "Caspase 1"
