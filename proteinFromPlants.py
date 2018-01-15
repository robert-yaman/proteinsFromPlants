import argparse
import csv
import importlib
import time
import os
import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], './objects'))

from sequenceSet import SequenceSet
import report

def findTransformations(start_seqs, transformation_modules, target_seq, 
	max_cost):
	lowest_score_per_set = {}
	# Initialize the cheapest cost to max_cost so we don't consider anything
	# more expensive.
	cheapest_success = [max_cost]

	def should_process_from_state(sequence_set, transformation_chain):
		if not sequence_set in lowest_score_per_set:
			return True

		return lowest_score_per_set[sequence_set] > total_cost(
			transformation_chain)

	def should_process_new_chain(transformation_chain):
		# Don't continue search if our cost already succeed the cheapest chain
		# to obtain the target.
		return total_cost(transformation_chain) < cheapest_success[0]

	def process(sequence_set, transformation_chain):
		# Returns transformation chain that yields the target seq, or None if
		# there is none.
		if total_cost(transformation_chain) > max_cost:
			print "ERROR: cost %d greater than max cost %d." % (
				total_cost(transformation_chain), max_cost)
		if not should_process_from_state(sequence_set, transformation_chain):
			return None
		lowest_score_per_set[sequence_set] = total_cost(transformation_chain)

		params = {
			"start_seqs" : start_seqs,
			"sequence_set" : sequence_set,
			# For tesing.
			"cheapest_success" : cheapest_success[0],
			"transformation_chain" : transformation_chain,
		}

		if sequence_set.isExactly(target_seq):
			total_cost_ = total_cost(transformation_chain)
			if total_cost_ < cheapest_success[0]:
				cheapest_success[0] = total_cost_
			return transformation_chain
		else:
			answer = None
			for transformation_module in transformation_modules:
				transformation_instances = (transformation_module.
					transformations(params))
				for transformation in transformation_instances:
					new_chain = transformation_chain + [transformation]
					if not should_process_new_chain(new_chain):
						continue
					child_answer = process(transformation.transform(
						sequence_set), new_chain)
					if not child_answer == None:
						if answer == None or len(child_answer) < len(answer):
							answer = child_answer
			return answer

	start_time = time.time()
	transformation_chain = process(SequenceSet(), [])
	time_elapsed = time.time() - start_time
	return transformation_chain, total_cost(transformation_chain), time_elapsed

def total_cost(transformation_chain):
	if transformation_chain == None:
		return 0
	else:
		return sum([transformation.cost() for transformation in transformation_chain])

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--transformations',
		help="Python module that holds the transformations.",
		required=True,
	)
	parser.add_argument(
		'--starters',
		help='CSV with the starter proteins. First column is sequence, second is cost per unit.',
		required=True,
	)
	parser.add_argument(
		'--target',
		help='TXT file with the target sequence.',
		required=True,
	)
	parser.add_argument(
		'--max-cost',
		default=7,
		type=int,
	)
	parser.add_argument(
		'--verbose',
		default=False,
		type=bool,
	)

	args = parser.parse_args()

	transformations = importlib.import_module(args.transformations)

	with open(args.target, 'r') as target_file:
		target_sequence = target_file.read()[0].rstrip()
		target_seq = Seq(target_sequence, generic_protein)

	start_seqs_path = args.starters
	with open(start_seqs_path, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		start_seqs = [[Seq(row[0], generic_protein), int(row[1]
			)] for row in reader]

	output = findTransformations(start_seqs, 
		transformations.TRANSFORMATIONS, target_seq, args.max_cost)
	
	print report.report(output[0], output[1], output[2], args.verbose)
