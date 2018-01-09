def short_report(transformation_chain, total_cost):
	string = "TRANSFORMATION CHAIN:\n"
	for transformation in transformation_chain:
		string += "    - %s\n" % transformation.name()
	string += "TOTAL COST: %s" % str(total_cost)
	return string

def long_report(transformation_chain):
	return "TODO: long report."

def report(transformation_chain, total_cost):
	"""Generates a human readable report about the transformation chain."""
	short = short_report(transformation_chain, total_cost)
	long_ = long_report(transformation_chain)
	return short + "\n ########## \n" + long_

