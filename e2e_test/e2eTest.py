import os
import sys
import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from sets import Set

sys.path.insert(1, os.path.join(sys.path[0], '..'))
sys.path.insert(1, os.path.join(sys.path[0], '../objects'))

import fakeTransformation

from proteinFromPlants import findTransformations
from proteinFromPlants import total_cost
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
		start_seqs = [[Seq("AFT", generic_protein), 1]]
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
		start_seqs = [[Seq("AFT", generic_protein), 1]]
		target_seq = Seq("L", generic_protein)
		output = findTransformations(start_seqs, self.transformations, target_seq, 6)

		self.assertEqual(output[0], None)

	def test_stops_search_when_cost_is_above_cheapest_success(self):
		# Since the fake transformation doesnt modify the sequence set, there 
		# shouldn't be two fake transformations with the same sequence set.
		start_seqs = [
			[Seq("AFT", generic_protein), 1], 
			[Seq("CFQ", generic_protein), 1],
		]
		target_seq = Seq("C", generic_protein)
		output = findTransformations(start_seqs, self.transformations, target_seq, 6)

		for item in self.work_list:
			params = item[1]
			self.assertTrue(params["cheapest_success"] >=
				total_cost(params["transformation_chain"]))


	def test_never_has_two_fake_transitions(self):
		# Since the fake transformation doesnt modify the sequence set, there 
		# shouldn't be two fake transformations with the same sequence set.
		start_seqs = [
			[Seq("AFT", generic_protein), 1], 
			[Seq("CFQ", generic_protein), 1],
		]
		target_seq = Seq("C", generic_protein)
		output = findTransformations(start_seqs, self.transformations, target_seq, 6)

		sequence_sets = [item[1]["sequence_set"] for item in self.work_list]
		self.assertEqual(len(sequence_sets), len(Set(sequence_sets)))


if __name__ == '__main__':
	unittest.main()
