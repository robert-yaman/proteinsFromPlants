# Todos

- Deal with exceptions in cleaving for Trypsin
- Efficiency: 
	- Try out non-deterministic algorithms for efficiency.
	- Parralel execution.
- Generate full report
- Purification step only happens at the end? This will greatly help efficiency since purification blows up as proteins cleave.
- Don't continue graph search if current cost is more than the lowest found cost!!!
	- then directly algs could speed up avg runtime
- Write E2E tests.
	- need a --verbose tag for this. Prints out each node we visit.
- Use max cost rather than max lenght

# Limitations

- Our algorithm won't find parallel workflows - e.g. apply this enzyme to this protein, then this other enzyme to another protein, then combine them.
- Our algorithm only supporte purifying one protein at a time (I'm not sure if this is an issue).