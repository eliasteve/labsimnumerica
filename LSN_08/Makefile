CXX = c++
CXXFLAGS = -std=c++11 -Wall --pedantic
LDFLAGS = -std=c++11 -Wall --pedantic 

lsn8.1: lsn8.1.o ../random/random.o ../lib_NSL/metropolis.o mcmc.o 
	$(CXX) $(LDFLAGS) $^ -o $@.out

lsn8.2: lsn8.2.o ../random/random.o ../lib_NSL/metropolis.o mcmc.o
	$(CXX) $(LDFLAGS) $^ -o $@.out

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o *.out ../random/*.o ../lib_NSL/*.o
