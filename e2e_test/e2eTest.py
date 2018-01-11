import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], '..'))
sys.path.insert(1, os.path.join(sys.path[0], '../objects'))

import fakeTransformation

from proteinFromPlants import findTransformations
from objects.transformations import addProtein
from objects.transformations import heatShock
from objects.transformations import purification

class E2eTest(unittest.TestCase):
	def setUp(self):
		self.work_list = fakeTransformation.WORK_LIST
		while len(self.work_list) > 0:
			self.work_list.pop()

		self.transformations = [
			addProtein,
			fakeTransformation,
			heatShock,
			purification,
		]

	def test_base_case(self):
		start_seqs = [Seq("AFT", generic_protein)]
		target_seq = Seq("A", generic_protein)
		output = findTransformations(start_seqs, self.transformations, target_seq, 6)
		transformation_chain = output[0]

		# Add AFT -> Heat shock -> purification
		self.assertEqual(len(transformation_chain), 3)

		first_transformation = transformation_chain[0]
		self.assertEqual(first_transformation.name(), "Add AFT")

		second_transformation = transformation_chain[1]
		self.assertEqual(second_transformation.name(), "Heat Shock")

		third_transformation = transformation_chain[2]
		self.assertEqual(third_transformation.name(), "Purify A")

	def test_impossible(self):
		start_seqs = [Seq("AFT", generic_protein)]
		target_seq = Seq("L", generic_protein)
		output = findTransformations(start_seqs, self.transformations, target_seq, 6)

		self.assertEqual(output[0], None)

if __name__ == '__main__':
	unittest.main()
