import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.argCProteinase import ArgCProteinase
from sequenceSet import SequenceSet

class ArgCProteinaseTest(unittest.TestCase):
	def setUp(self):
		self.argC = ArgCProteinase()

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFRCT"), 1)
		output_set = self.argC.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATFR")), 1)
		self.assertEqual(output_set.frequency(self._sequence("CT")), 1)

	def test_two_cleaves(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFRCTRY"), 1)
		output_set = self.argC.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATFR")), 1)
		self.assertEqual(output_set.frequency(self._sequence("CTR")), 1)
		self.assertEqual(output_set.frequency(self._sequence("Y")), 1)

	def test_two_close_cleaves(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFRRY"), 1)
		output_set = self.argC.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATFR")), 1)
		self.assertEqual(output_set.frequency(self._sequence("R")), 1)
		self.assertEqual(output_set.frequency(self._sequence("Y")), 1)

	def test_doesnt_cleave_at_end(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("TTQ"), 1)
		output_set = self.argC.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("TTQ")), 1)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
