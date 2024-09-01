CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall --pedantic
LDFLAGS = -std=c++11 -Wall --pedantic -larmadillo

lsn10.1: lsn10.1.o ../random/random.o ../lib_NSL/individual.o ../lib_NSL/optimization.o
	$(CXX) $(LDFLAGS) $^ -o $@.out

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
