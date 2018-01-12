# Todos

- Deal with exceptions in cleaving for Trypsin
- Efficiency: 
	- Try out non-deterministic algorithms for efficiency.
	- Parralel execution.
- Generate full report
- Purification step only happens at the end? This will greatly help efficiency since purification blows up as proteins cleave.
- Use max cost rather than max length
- Use smart alg to figure out which node to process next. Because of early returning, this will speed up avg case.

# Limitations

- Our algorithm won't find parallel workflows - e.g. apply this enzyme to this protein, then this other enzyme to another protein, then combine them.
- Our algorithm only supporte purifying one protein at a time (I'm not sure if this is an issue).