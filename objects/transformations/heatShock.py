from Bio.Alphabet import generic_protein
from Bio.Seq import Seq

import transformation
from cost import Cost
from sequenceSet import SequenceSet

class HeatShock(transformation.Transformation):
	def transform(self, input_aas):
		returnAAs = SequenceSet()
		for seq in input_aas.sequences():
			for aa in seq:
				aa_seq = self._sequence(aa)
				returnAAs.setFrequency(aa_seq, 
					returnAAs.frequency(aa_seq) + input_aas.frequency(seq))
		return returnAAs

	def cost(self):
		return Cost(100)

	def name(self):
		return "Heat Shock"

	def _sequence(self, sequence):
		return Seq(sequence, generic_protein)