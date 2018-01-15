import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.simpleLigase import SimpleLigase
from sequenceSet import SequenceSet

class SimpleLigaseTest(unittest.TestCase):
	def test_base_case(self):
		self.transformation = SimpleLigase('test ligase', '---W-D---', 0)
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAW"), 1)
		inputs_set.setFrequency(self._sequence("D"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAWD")), 1)

	def test_open_back_condition_case(self):
		self.transformation = SimpleLigase('test ligase', '---W----', 0)
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAW"), 1)
		inputs_set.setFrequency(self._sequence("D"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAWD")), 1)
		self.assertEqual(output_set.frequency(self._sequence("AFAWAFAW")), 1)

	def test_negative_condition(self):
		self.transformation = SimpleLigase('test ligase', '---W-^DA---', 0)
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAW"), 1)
		inputs_set.setFrequency(self._sequence("D"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAW")), 1)
		self.assertEqual(output_set.frequency(self._sequence("D")), 1)

	def test_negative_back_condition(self):
		self.transformation = SimpleLigase('test ligase', '--^A-W-D---', 0)
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("AFAW"), 1)
		inputs_set.setFrequency(self._sequence("D"), 1)
		output_set = self.transformation.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("AFAW")), 1)
		self.assertEqual(output_set.frequency(self._sequence("D")), 1)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
