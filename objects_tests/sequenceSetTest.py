import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '../objects/'))
from sequenceSet import SequenceSet

class SequenceSetTest(unittest.TestCase):
	def test_equality(self):
		s1 = SequenceSet()
		s2 = SequenceSet()
		s3 = SequenceSet()
		s4 = SequenceSet()

		s1.setFrequency(self._sequence("A"), 1)
		s2.setFrequency(self._sequence("A"), 1)
		s3.setFrequency(self._sequence("A"), 2)
		s4.setFrequency(self._sequence("A"), 1)
		s4.setFrequency(self._sequence("F"), 1)

		self.assertEqual(hash(s1), hash(s2))
		self.assertNotEqual(hash(s1), hash(s3))
		self.assertNotEqual(hash(s1), hash(s4))
		self.assertTrue(s1 == s2)
		self.assertFalse(s1 == s3)
		self.assertFalse(s1 == s4)

	def _sequence(self, seq):
		return Seq(seq, generic_protein)

if __name__ == '__main__':
	unittest.main()
