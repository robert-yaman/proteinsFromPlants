import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.purification import Purification
from sequenceSet import SequenceSet

class purificationTest(unittest.TestCase):
	def setUp(self):
		target_seq = Seq("A", generic_protein)
		self.purification = Purification(target_seq)

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("FTF"), 1)
		inputs_set.setFrequency(self._sequence("A"), 2)
		output_set = self.purification.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("FTF")), 0)
		self.assertEqual(output_set.frequency(self._sequence("A")), 2)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
