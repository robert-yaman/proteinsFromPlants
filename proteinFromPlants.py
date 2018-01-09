import argparse
import importlib
import os
import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

sys.path.insert(1, os.path.join(sys.path[0], './objects'))

from sequenceSet import SequenceSet

def main(start_seq, transformations, target_seq, max_length):
	def process(sequence_set, transformation_chain):
		# Returns transformation chain that yields the target seq, or None if
		# there is none.
		if len(transformation_chain) > max_length:
			print "ERROR: length %d greater than max length %d." % (
				len(transformation_chain), max_length)
		
		if sequence_set.contains(target_seq):
			return transformation_chain
		elif len(transformation_chain) == max_length:
			return None
		else:
			answer = None
			for transformation_type in transformations:
				transformation = transformation_type()
				new_chain = transformation_chain + [transformation]
				child_answer = process(transformation.transform(sequence_set), 
					new_chain)
				if not child_answer == None and len(child_answer) < len(answer):
					answer = child_answer
			return answer

	start_seq_set = SequenceSet()
	start_seq_set.setFrequency(start_seq, 1)
	print process(start_seq_set, [])






if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--transformations',
		help="Module that holds the transformations.",
		required=True,
	)
	parser.add_argument(
		'--max_length',
		default=10,
		type=int,
	)

	args = parser.parse_args()
	transformations = importlib.import_module(args.transformations)
	# Load passed in module.

	start_seq = Seq("MASLPWSLTTSTAIANTTNISAFPPSPLFQRASHVPVARNRSRRFAPSKVSCNSANGDPNSDSTSDVRETSSGKLDRRNVLLGIGGLYGAAGGLGATKPLAFGAPIQAPDISKCGTATVPDGVTPTNCCPPVTTKIIDFQLPSSGSPMRTRPAAHLVSKEYLAKYKKAIELQKALPDDDPRSFKQQANVHCTYCQGAYDQVGYTDLELQVHASWLFLPFHRYYLYFNERILAKLIDDPTFALPYWAWDNPDGMYMPTIYASSPSSLYDEKRNAKHLPPTVIDLDYDGTEPTIPDDELKTDNLAIMYKQIVSGATTPKLFLGYPYRAGDAIDPGAGTLEHAPHNIVHKWTGLADKPSEDMGNFYTAGRDPIFFGHHANVDRMWNIWKTIGGKNRKDFTDTDWLDATFVFYDENKQLVKVKVSDCVDTSKLRYQYQDIPIPWLPKNTKAKAKTTTKSSKSGVAKAAELPKTTISSIGDFPKALNSVIRVEVPRPKKSRSKKEKEDEEEVLLIKGIELDRENFVKFDVYINDEDYSVSRPKNSEFAGSFVNVPHKHMKEMKTKTNLRFAINELLEDLGAEDDESVIVTIVPRAGGDDVTIGGIEIEFVSD", generic_protein)
	target_seq = Seq("YQPPSTNKNTKSQRRKGSTFEEHK", generic_protein)
	main(start_seq, transformations.TRANSFORMATIONS, 
		target_seq, args.max_length)
