CXX = c++
CFLAGS = -Wall --pedantic -std=c++11

lsn2.2 : lsn2.2.o ../random/random.o ../lib_NSL/lib_NSL.o
	$(CXX) $(CFLAGS) $^ -o lsn2.2.out

lsn2.1 : lsn2.1.o ../random/random.o ../lib_NSL/lib_NSL.o
	$(CXX) $(CFLAGS) $^ -o lsn2.1.out

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@

