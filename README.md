# SubSieve
Implementation of the SubSieve algorithm, following the article
[Léo Ducas, Shortest Vector from Lattice Sieving: a Few Dimensions for Free]

It is written in C++ and python, relying on the [fpylll library](https://github.com/fplll/fpylll). It also requires the [ctypes](http://www.python.net/crew/theller/ctypes/) and [numpy](http://www.numpy.org/) package. The C++ component is used through middleware.py, that wrap the shared object binary SubSieveLib.so, using python ctypes. 

To use this SubSieve.py, you must:
- be in your fpylll environment (See item 7. of [this README](https://github.com/fplll/fpylll#getting-started))
- compile the C++ coponent using Make

To reproduce the benchmarks of the paper, simply run bench.py. Of course, results may differs depending on your hardware, or on the version of fpylll used.



