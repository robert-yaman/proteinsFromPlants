# Proteins From Plants

This project implements an algorithm that finds a sequence of enzymes to a set
of starter proteins that yields a target protein. The intended application is
to find a cheap way to synthesize proteins important for animal muscle cell
development by using plant proteins.

Important concepts:

- A sequence set is a collection of peptide sequences, and their respective
frequencies (relative to the starting unit).
	- Defined in ./objects/sequenceSet.py
- A transformation is some mapping between sequence sets. It represents a step
in an experimental procedure, e.g. apply a particular ligase, or purify a
particular protein.
	- Base class defined in ./objects/transformations/transformation.py

# Usage

The base script is proteinFromPlants.py. Run this script passing in the 
following arguments:

- The path to a CSV file containing the starter proteins. The first column is
the sequence of the protein, and the second column is the cost per unit.
./testStarters.csv is an example file.
- A TXT file containing the sequence of the target protein. ./testTarget.txt is
an example file.
- A python module containing all of the possible transformations. The set of 
all transformations so far is objects.transformations.allTransformations. The
module should contain a python list named TRANSFORMATIONS contains all of the
transformation modules. Each transformation module should have a method 
`transformations` that returns a list of individual transformation objects. 
A `params` dictionary is passed into the `transformations` method that can
assist in transformation object creation.

The usage cleaving enzymes and ligases are defined in ./cleavers.csv and 
./ligases.csv respectively. The columns of the csvs are defined as follows:
- The name of the enzyme
- The enzyme's code. This code signifies the 8-length amino acid sequence that 
the enzyme will act on. The middle (4th) position is the cleaving position, or
the position at which the proteins will be binded. Each position can have 
either positive conditions or negative conditions (signified by '^'). If there
is a positive conditions, then the sequence at that position must be one of the
listed amino acids. If there is a negative condition, the sequence at that
position must not be any of the listed amino acids.
- The cost
- A link to more info about the enzyme
- Comments

An example usage can be run by:

`bash run.bash`

# Limitations

- Our algorithm won't find parallel workflows - e.g. apply an enzyme to a protein, then another enzyme to another protein, then combine them.
- Our algorithm only support purifying one protein at a time (I'm not sure if it's physically possible to purify multiple proteins ata time ).

# Todos

- Deal with exceptions in cleaving for Trypsin
- Efficiency: 
	- Try out non-deterministic algorithms for efficiency.
	- Parallel execution.
	- Mode that ignores frequency.
- Generate full report
- Use smart alg to figure out which node to process next. Because of early returning, this will speed up avg case.
	- Purification step only happens at the end?
		- Or, if our target is a subset of the sequence set, look at 
		purifications first.
	- If we don't have the total AAs we need, only perform "additive" 
	transformations. We can define this dynamically in init - if giving an 
	empty set yields a non-empty set.
- Improve pose-ligase sequence quantities.
