lsn5.1 : lsn5.1.o ../random/random.o ../random/random.h
	c++ -std=c++11 -Wall --pedantic lsn5.1.o ../random/random.o -o lsn5.1.out

../random/random.o : ../random/random.cpp ../random/random.h
	c++ -std=c++11 -Wall --pedantic -c ../random/random.cpp -o ../random/random.o

lsn5.1.o : lsn5.1.cpp
	c++ -std=c++11 -Wall --pedantic -c lsn5.1.cpp -o lsn5.1.o
