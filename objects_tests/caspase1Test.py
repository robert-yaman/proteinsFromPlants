import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.caspase1 import Caspase1
from sequenceSet import SequenceSet

class Caspase1Test(unittest.TestCase):
	def setUp(self):
		self.transformation = Caspase1()

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAADTW"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAAD")), 1)
		self.assertEqual(output_set.frequency(self._sequence("TW")), 1)

	def test_missing_initial_condition(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATAADTW"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATAADTW")), 1)
		self.assertEqual(output_set.frequency(self._sequence("TW")), 0)

	def test_negative_condition(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAADEW"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAADEW")), 1)
		self.assertEqual(output_set.frequency(self._sequence("EW")), 0)


	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
