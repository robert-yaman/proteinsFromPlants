import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from transformations.aspNEndopeptidase import AspNEndopeptidase
from sequenceSet import SequenceSet

class AspNEndopeptidaseTest(unittest.TestCase):
	def setUp(self):
		self.aspN = AspNEndopeptidase()

	def test_base_case(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFDCT"), 1)
		output_set = self.aspN.transform(inputs_set)
		print output_set.sequences()
		self.assertEqual(output_set.frequency(self._sequence("ATF")), 1)
		self.assertEqual(output_set.frequency(self._sequence("DCT")), 1)

	def test_both_cleave_sites(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("ATFDCQT"), 1)
		output_set = self.aspN.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("ATF")), 1)
		self.assertEqual(output_set.frequency(self._sequence("DC")), 1)
		self.assertEqual(output_set.frequency(self._sequence("QT")), 1)

	def test_cleaves_at_end(self):
		inputs_set = SequenceSet()
		inputs_set.setFrequency(self._sequence("TTQ"), 1)
		output_set = self.aspN.transform(inputs_set)
		self.assertEqual(output_set.frequency(self._sequence("TT")), 1)
		self.assertEqual(output_set.frequency(self._sequence("Q")), 1)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
