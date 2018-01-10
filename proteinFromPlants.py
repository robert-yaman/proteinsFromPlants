import argparse
import importlib
import time
import os
import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], './objects'))

from sequenceSet import SequenceSet
import report

def main(start_seqs, transformation_modules, target_seq, max_length):
	lowest_score_per_set = {}
	def stop_search(sequence_set, transformation_chain):
		if not sequence_set in lowest_score_per_set:
			return False
		return lowest_score_per_set[sequence_set] <= total_cost(
			transformation_chain)

	def process(sequence_set, transformation_chain):
		# Returns transformation chain that yields the target seq, or None if
		# there is none.
		if len(transformation_chain) > max_length:
			print "ERROR: length %d greater than max length %d." % (
				len(transformation_chain), max_length)
		
		if stop_search(sequence_set, transformation_chain):
			return None
		lowest_score_per_set[sequence_set] = total_cost(transformation_chain)

		params = {
			"start_seqs" : start_seqs,
			"sequence_set" : sequence_set,
		}

		if sequence_set.isExactly(target_seq):
			return transformation_chain
		elif len(transformation_chain) == max_length:
			return None
		else:
			answer = None
			for transformation_module in transformation_modules:
				transformation_instances = (transformation_module.
					transformations(params))
				for transformation in transformation_instances:
					new_chain = transformation_chain + [transformation]
					child_answer = process(transformation.transform(
						sequence_set), new_chain)
					if not child_answer == None:
						if answer == None or len(child_answer) < len(answer):
							answer = child_answer
			return answer

	start_time = time.time()
	transformation_chain = process(SequenceSet(), [])
	time_elapsed = time.time() - start_time
	print report.report(transformation_chain, 
		total_cost(transformation_chain), time_elapsed)

def total_cost(transformation_chain):
	if transformation_chain == None:
		return 0
	else:
		return sum([transformation.cost() for transformation in transformation_chain])

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--transformations',
		help="Module that holds the transformations.",
		required=True,
	)
	parser.add_argument(
		'--max-length',
		default=10,
		type=int,
	)

	args = parser.parse_args()
	transformations = importlib.import_module(args.transformations)
	# Load passed in module.

	start_seq = Seq("MASLPWSLTTSTAIANTTNISAFPPSPLFQRASHVPVARNRSRRFAPSKVSCNSANGDPNSDSTSDVRETSSGKLDRRNVLLGIGGLYGAAGGLGATKPLAFGAPIQAPDISKCGTATVPDGVTPTNCCPPVTTKIIDFQLPSSGSPMRTRPAAHLVSKEYLAKYKKAIELQKALPDDDPRSFKQQANVHCTYCQGAYDQVGYTDLELQVHASWLFLPFHRYYLYFNERILAKLIDDPTFALPYWAWDNPDGMYMPTIYASSPSSLYDEKRNAKHLPPTVIDLDYDGTEPTIPDDELKTDNLAIMYKQIVSGATTPKLFLGYPYRAGDAIDPGAGTLEHAPHNIVHKWTGLADKPSEDMGNFYTAGRDPIFFGHHANVDRMWNIWKTIGGKNRKDFTDTDWLDATFVFYDENKQLVKVKVSDCVDTSKLRYQYQDIPIPWLPKNTKAKAKTTTKSSKSGVAKAAELPKTTISSIGDFPKALNSVIRVEVPRPKKSRSKKEKEDEEEVLLIKGIELDRENFVKFDVYINDEDYSVSRPKNSEFAGSFVNVPHKHMKEMKTKTNLRFAINELLEDLGAEDDESVIVTIVPRAGGDDVTIGGIEIEFVSD", generic_protein)
	start_seq2 = Seq("YQPPSTNKNTKSQRRKGSTFEEHK", generic_protein)
	target_seq = Seq("F", generic_protein)
	main([start_seq, start_seq2], transformations.TRANSFORMATIONS, 
		target_seq, args.max_length)
