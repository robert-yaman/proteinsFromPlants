# Todos

- Deal with exceptions in cleaving for Trypsin
- Efficiency: 
	- Try out non-deterministic algorithms for efficiency.
	- Parralel execution.
	- Mode that ignores frequency.
- Generate full report
- Use smart alg to figure out which node to process next. Because of early returning, this will speed up avg case.
	- Purification step only happens at the end?
		- Or, if our target is a subset of the sequence set, look at 
		purifications first.
	- If we don't have the total AAs we need, only perform "addititive" 
	transformations. We can define this dynamically in init - if giving an 
	empty set yields a non-empty set.
- Make a more nuanced picture of the quantites after a ligase is applied.
- Dynamic inputs for script.

# Limitations

- Our algorithm won't find parallel workflows - e.g. apply this enzyme to this protein, then this other enzyme to another protein, then combine them.
- Our algorithm only supporte purifying one protein at a time (I'm not sure if it's physically possible to purify multiple proteins ata time ).