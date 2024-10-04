#include <iostream>
#include <fstream>
#include <cstdlib> // EXIT_FAILURE
#include "../random/random.h"
#include "../lib_NSL/misc.h"

int throwNeedle(Random& rng, double L, double d);

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  double L = 1.; //Length of the needle
  double d = 1.1; //Distance of the "strips"
  int M = 10000; //Throws per block
  int nBlocks = 100; //Number of blocks
  int hits = 0; //Number of times the needle crosses the line in a block
  double blockValue = 0; //Value of pi in every block
  //Accumulators for mean, mean^2 of the block values
  double meanAccumulator = 0, mean2Accumulator = 0;
  std::ofstream outf("pi.dat"); //Output file

  if(!outf.is_open()) {
    std::cerr << "I/O error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < nBlocks; i++) {
    hits = 0;
    for (int j = 0; j < M; j++) {
      hits += throwNeedle(rng, L, d);
    }
    //Value of pi from number of throws and relevant parameters
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

//Throw the needle once. Returns 1 if the needle "crosses the line",
//0 otherwise
int throwNeedle(Random& rng, double L, double d) {
  double p = rng.Rannyu(0, d); //Position of one end on the needle
  double x, y; //Relative coordinates of the other end of the needle
  //Norm of the randomly generated point (for sampling the angle)
  double norm;

  //Select an angle uniformly between -pi and pi by sampling a point in
  //"the first two quadrants" and keeping it only if it's inside the unit
  //disk. This implementation could technically run in an infinite loop,
  //but this happens with probability zero, and is unavoidable provided
  //we always want to generate a point per function call.
  do {
    x = rng.Rannyu(), y = rng.Rannyu(-1, 1);
    norm = pow(x, 2) + pow(y, 2);
  } while (norm > 1);

  // Normalize to reflect the length of the needle. No need to normalize
  // y as it's not used in the following
  x *= L/sqrt(norm);

  if(x+p >= d) return 1;
  else return 0;
}
