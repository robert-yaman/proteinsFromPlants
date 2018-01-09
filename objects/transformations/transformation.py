import abc

from Bio.Alphabet import generic_protein
from Bio.Seq import Seq

class Transformation:
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def transform(self, input_aas):
		pass

	@abc.abstractmethod
	def cost(self):
		pass

	@abc.abstractmethod
	def name(self):
		pass

	def _sequence(self, sequence):
		return Seq(sequence, generic_protein)