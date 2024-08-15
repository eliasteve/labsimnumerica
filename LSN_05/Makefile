CXX = c++
CFLAGS = -Wall --pedantic -std=c++11

lsn5.1 : lsn5.1.o ../random/random.o ../lib_NSL/misc.o ../lib_NSL/metropolis.o ../lib_NSL/point.o
	$(CXX) $(CFLAGS) $^ -o lsn5.1.out

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@

