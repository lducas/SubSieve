all:
	g++ -fPIC -Ofast -mavx2 -funroll-loops -std=c++11 -c SubSieveLib.cpp -o SubSieveLib.o
	g++ -shared -Ofast -mavx2 -funroll-loops -std=c++11 SubSieveLib.o -o SubSieveLib.so

