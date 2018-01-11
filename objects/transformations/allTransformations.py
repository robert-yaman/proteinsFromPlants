import addProtein
import argCProteinase
import BNPSSkatole
import caspase1
import aspNEndopeptidase
import heatShock
import simpleCleaver
import purification

TRANSFORMATIONS = [
	addProtein,
	# Now in simply cleavers:
	# argCProteinase,
	# aspNEndopeptidase,
	# BNPSSkatole,
	# caspase1,
	heatShock,
	simpleCleaver,
	purification,
]