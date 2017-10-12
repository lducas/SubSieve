# SubSieve
Implementation of the SubSieve algorithm, following the article 
[Shortest Vector from Lattice Sieving: a Few Dimensions for Free](https://eprint.iacr.org/2017/999)

It is written in C++ and python, relying on the python [fpylll library](https://github.com/fplll/fpylll). It also requires the [ctypes](http://www.python.net/crew/theller/ctypes/) and [numpy](http://www.numpy.org/) python packages. The C++ component is used through middleware.py, that wrap the shared object binary SubSieveLib.so, using python ctypes. 

To use this SubSieve.py, you must first:
- compile the C++ component using Make
- be in your fpylll environment (See item 7. of [this howto](https://github.com/fplll/fpylll#getting-started))

The current makefile is very specific, and should be adjusted to your hardware and platforms. In particular, for OSX, you may refer to this [issue](https://github.com/lducas/SubSieve/issues/1).

To reproduce the benchmarks of the paper, simply run bench.py. Of course, results may differs depending on your hardware, or on the version of fpylll used.
