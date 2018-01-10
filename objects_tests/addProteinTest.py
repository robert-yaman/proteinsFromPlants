import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.addProtein import AddProtein
from sequenceSet import SequenceSet

class AddProteinTest(unittest.TestCase):
	def setUp(self):
		self.add_protein = AddProtein(Seq("FTF", generic_protein))

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("A"), 2)
		output_set = self.add_protein.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("FTF")), 1)
		self.assertEqual(output_set.frequency(self._sequence("A")), 2)

	def test_adjusts_freq(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("A"), 2)
		inputs_set.setFrequency(self._sequence("FTF"), 1)
		output_set = self.add_protein.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("FTF")), 2)
		self.assertEqual(output_set.frequency(self._sequence("A")), 2)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
