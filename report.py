def short_report(transformation_chain, total_cost, time_elapsed):
	string = "TRANSFORMATION CHAIN:\n"
	for transformation in transformation_chain:
		string += "    - %s\n" % transformation.name()
	string += "TOTAL COST: %s\n" % str(total_cost)
	string += "TIME ELAPSED: %f" % time_elapsed
	return string

def long_report(transformation_chain):
	return "TODO: long report."

def report(transformation_chain, total_cost, time_elapsed):
	"""Generates a human readable report about the transformation chain."""
	short = short_report(transformation_chain, total_cost, time_elapsed)
	long_ = long_report(transformation_chain)
	return short + "\n ########## \n" + long_

