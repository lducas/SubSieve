all:
	g++ -fPIC -Ofast -march=native -std=c++11 -c SubSieveLib.cpp -o SubSieveLib.o
	g++ -shared -Ofast -march=native -std=c++11 SubSieveLib.o -o SubSieveLib.so
	g++ -Ofast -march=native -std=c++11 SubSieveLibTest.cpp -o SubSieveLibTest