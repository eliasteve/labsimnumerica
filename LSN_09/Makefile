CXX = c++
CXXFLAGS = -std=c++11 -Wall --pedantic
LDFLAGS = -std=c++11 -Wall --pedantic -larmadillo

lsn9.1 : lsn9.1.o ../random/random.o ../lib_NSL/individual.o ../lib_NSL/optimization.o
	$(CXX) $(LDFLAGS) $^ -o $@.out

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
