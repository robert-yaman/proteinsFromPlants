import abc

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