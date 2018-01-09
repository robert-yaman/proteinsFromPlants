import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.heatShock import HeatShock
from sequenceSet import SequenceSet

class HeatShockTest(unittest.TestCase):
	def setUp(self):
		self.heatShock = HeatShock()

	def test_transform(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("FTF"), 1)
		output_set = self.heatShock.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("F")), 2)
		self.assertEqual(output_set.frequency(self._sequence("T")), 1)

	def test_different_sequences(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("FTF"), 1)
		inputs_set.setFrequency(self._sequence("AFA"), 1)
		output_set = self.heatShock.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("F")), 3)
		self.assertEqual(output_set.frequency(self._sequence("A")), 2)
		self.assertEqual(output_set.frequency(self._sequence("T")), 1)

	def test_double_sequences(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("FTF"), 1)
		inputs_set.setFrequency(self._sequence("AFA"), 2)
		output_set = self.heatShock.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("F")), 4)
		self.assertEqual(output_set.frequency(self._sequence("A")), 4)
		self.assertEqual(output_set.frequency(self._sequence("T")), 1)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
