CXX = c++
CFLAGS = -Wall --pedantic -std=c++11

lsn3.1 : lsn3.1.o ../random/random.o ../lib_NSL/misc.o
	$(CXX) $(CFLAGS) $^ -o lsn3.1.out

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@

