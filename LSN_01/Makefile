CXX = c++
CFLAGS = -Wall --pedantic -std=c++11

lsn1.3 : lsn1.3.o ../random/random.o ../lib_NSL/lib_NSL.o
	$(CXX) $(CFLAGS) $^ -o lsn1.3.out

lsn1.2 : lsn1.2.o ../random/random.o ../random/random.h
	$(CXX) $(CFLAGS) lsn1.2.o ../random/random.o -o lsn1.2.out

lsn1_01.1 : lsn1_01.1.o ../random/random.o ../random/random.h
	$(CXX) $(CFLAGS) lsn1_01.1.o ../random/random.o -o lsn1_01.1.out

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c $< -o $@

%.o: %.cpp 
	$(CXX) $(CFLAGS) -c $< -o $@

