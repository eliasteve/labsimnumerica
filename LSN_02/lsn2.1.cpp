#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"
#include "../lib_NSL/lib_NSL.h"

double integrand(double x) {
  return M_PI/2*cos(M_PI*x/2);
}

double myPDF(double x) {
  return 3./2*(1-pow(x, 2));
}

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  int nBlocks = 100, throwsPerBlock = 1000;
  double blockAccumulator, val = 0, meanAccum = 0, mean2Accum = 0;

  std::ofstream outf("lsn2.1_uniform.dat");
  if (!outf.is_open()) {
    std::cerr << "Error trying to open lsn2.1_uniform.dat. Terminating.";
    return -1;
  }

  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < throwsPerBlock; j++) {
      val = rng.Rannyu();
      blockAccumulator += integrand(val);
    }
    blockAccumulator /= throwsPerBlock;
    computeUpdateMeanAndError(i, blockAccumulator, meanAccum, mean2Accum, outf);
  }

  outf.close();
  meanAccum=0;
  mean2Accum=0;

  /*std::ofstream test("test.dat");

  for (int i = 0; i < 10000; i++) {
    test << rng.AR2_1() << std::endl;
  }*/

  outf.open("lsn2.1_mine.dat");
  if (!outf.is_open()) {
    std::cerr << "Error trying to open lsn2.1_mine.dat. Terminating.";
    return -1;
  }

  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < throwsPerBlock; j++) {
      val = rng.AR2_1();
      blockAccumulator += integrand(val)/myPDF(val);
    }

    blockAccumulator /= throwsPerBlock;
    computeUpdateMeanAndError(i, blockAccumulator, meanAccum, mean2Accum, outf);
  }

  return 0;
}
