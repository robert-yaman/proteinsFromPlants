import transformation
from sequenceSet import SequenceSet
from simpleCleaver import SimpleCleaver

def transformations(params):
	return [BNPSSkatole()]

class BNPSSkatole(SimpleCleaver):
	"""Source: http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions"""
	def __init__(self):
		return super(BNPSSkatole, self).__init__("BNPS-Skatole", '---W----', 1)