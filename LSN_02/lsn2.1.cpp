#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"
#include "../lib_NSL/misc.h"

//Integrand function
double integrand(double x) {
  return M_PI/2 * cos(M_PI * x / 2);
}

// Probability density function (PDF) used for importance sampling,
// P(x) = 3/2 * (1 - x^2) over the interval [0, 1]
double myPDF(double x) {
  return 3.0 / 2.0 * (1 - pow(x, 2));
}

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  // Number of blocks for the Monte Carlo integration
  int nBlocks = 100;
  int throwsPerBlock = 1000; // Number of throws per block
  double blockAccumulator = 0; //Accumulator for block value
  double meanAccum = 0; //Global accumulator for the mean
  double mean2Accum = 0; //Global accumulator for the mean^2
  double val = 0; //Stores the output by the RNG

  // Open a file to write the results for uniform sampling
  std::ofstream outf("lsn2.1_uniform.dat");
  if (!outf.is_open()) {
    std::cerr << "Error trying to open lsn2.1_uniform.dat. Terminating.";
    return -1;
  }

  // Perform Monte Carlo integration using uniform sampling
  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < throwsPerBlock; j++) {
      val = rng.Rannyu();  // Generate a random number uniformly between 0 and 1
      blockAccumulator += integrand(val);  // Evaluate the integrand at this point and accumulate
    }
    blockAccumulator /= throwsPerBlock;  // Average over the block
    // Update the running mean and error estimates, and write to the file
    computeUpdateMeanAndError(i, blockAccumulator, meanAccum, mean2Accum, outf);
  }

  outf.close();  // Close the file after writing uniform sampling results
  //Reset global accumulators
  meanAccum = 0;
  mean2Accum = 0;

  /* Test code for generating random numbers using custom distribution
  std::ofstream test("test.dat");
  for (int i = 0; i < 10000; i++) {
    test << rng.AR2_1() << std::endl;
  }
  */

  // Open a file to write the results for importance sampling
  outf.open("lsn2.1_mine.dat");
  if (!outf.is_open()) {
    std::cerr << "Error trying to open lsn2.1_mine.dat. Terminating.";
    return -1;
  }

  // Perform Monte Carlo integration using importance sampling
  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < throwsPerBlock; j++) {
      val = rng.AR2_1();  // Generate a random number according to the custom PDF
      blockAccumulator += integrand(val) / myPDF(val);  // Correct the integrand value by the PDF
    }

    blockAccumulator /= throwsPerBlock;  // Average over the block
    // Update the running mean and error estimates, and write to the file
    computeUpdateMeanAndError(i, blockAccumulator, meanAccum, mean2Accum, outf);
  }

  return 0;
}
