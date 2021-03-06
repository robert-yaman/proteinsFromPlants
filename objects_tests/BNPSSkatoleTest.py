import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.BNPSSkatole import BNPSSkatole
from sequenceSet import SequenceSet

class BNPSSkatoleTest(unittest.TestCase):
	def setUp(self):
		self.transformation = BNPSSkatole()

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFWCT"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATFW")), 1)
		self.assertEqual(output_set.frequency(self._sequence("CT")), 1)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
