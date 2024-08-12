CXX = c++
CFLAGS = -Wall --pedantic -std=c++11

lsn1.3 : lsn1.3.o ../random/random.o ../lib_NSL/misc.o
	$(CXX) $(CFLAGS) $^ -o lsn1.3.out

lsn1.2 : lsn1.2.o ../random/random.o 
	$(CXX) $(CFLAGS) $^ -o lsn1.2.out

lsn1.1 : lsn1.1.o ../random/random.o ../lib_NSL/misc.o
	$(CXX) $(CFLAGS) $^ -o lsn1.1.out

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@

