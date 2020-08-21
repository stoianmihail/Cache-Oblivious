# Cache-Oblivious
Cache-oblivious Data Structures, as described in MIT 6.851 Advanced Data Structures, Spring 2012

1) Static cache-oblivious search trees (see Solver for the main idea and BoostedSolver for a faster version) - It's indeed faster. For 2^26 ~ 67M elements one lookup costs:
Binary-Search: 1103 ns
Solver: 875 ns
BoostedSolver: 857 ns
