#include <iostream>
#include <fstream>
#include <cstdlib> // EXIT_FAILURE
#include "../random/random.h"
#include "../lib_NSL/lib_NSL.h"

int throwNeedle(Random& rng, double L, double d);

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  double L = 1., d = 1.1;
  int M = 10000, nBlocks = 100, hits = 0;
  double blockValue = 0, meanAccumulator = 0, mean2Accumulator = 0;
  std::ofstream outf("pi.dat");

  if(!outf.is_open()) {
    std::cerr << "I/O error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < nBlocks; i++) {
    hits = 0;
    for (int j = 0; j < M; j++) {
      hits += throwNeedle(rng, L, d);
    }
    blockValue = 2*L*M/(d*hits);
    computeUpdateMeanAndError(
      i,
      blockValue,
      meanAccumulator,
      mean2Accumulator,
      outf    
    );
  }

  return 0;
}

int throwNeedle(Random& rng, double L, double d) {
  double p = rng.Rannyu(0, d); //Position of one end on the needle
  double x, y; //Relative coordinates of the other end of the needle
  double norm;

  //Select an angle uniformly between -pi and pi by sampling a point in
  //"the first two quadrant" and keeping it only if it's inside the unit
  //disk. This implementation could technically run in an infinite loop,
  //but this happens with probability zero, and is unavoidable provided
  //we always want to generate a point per function call.
  do {
    x = rng.Rannyu(), y = rng.Rannyu(-1, 1);
    norm = pow(x, 2) + pow(y, 2);
  } while (norm > 1);

  // Normalize to reflect the length of the needle
  x *= L/sqrt(norm);

  if(x+p >= d) return 1;
  else return 0;
}
