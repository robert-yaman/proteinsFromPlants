from objects.transformations import transformation

# Mutated to keep track of what work any fake transition has done.
WORK_LIST = []

def transformations(params):
	return [FakeTransformation(params)]

class FakeTransformation(transformation.Transformation):
	def __init__(self, params):
		self.params = params

	def transform(self, input_aas):
		WORK_LIST.append([input_aas, self.params])
		return input_aas

	def cost(self):
		return 1

	def name(self):
		return "Fake Transformation"
