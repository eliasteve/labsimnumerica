#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"

double integrand(double x) {
  return M_PI/2*cos(M_PI*x/2);
}

double myPDF(double x) {
  return 3./2*(1-pow(x, 2));
}

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  int nBlocks = 100, throwsPerBlock = 1000;
  double blockAccumulator, val = 0, meanAccumUnif = 0, mean2AccumUnif = 0, meanUnif = 0;

  std::ofstream uniform("lsn2.1_uniform.dat");
  if (!uniform.is_open()) {
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
    meanAccumUnif += blockAccumulator;
    mean2AccumUnif += pow(blockAccumulator, 2);
    meanUnif = meanAccumUnif/(i+1);
    uniform << meanUnif << " ";
    uniform << pow((mean2AccumUnif/(i+1) - pow(meanUnif, 2))/i, 0.5) << std::endl; //Error nans if i = 0, which indicates an undefined sample mean
  }

  uniform.close();

  /*std::ofstream test("test.dat");

  for (int i = 0; i < 10000; i++) {
    test << rng.AR2_1() << std::endl;
  }*/

  std::ofstream mine("lsn2.1_mine.dat");
  if (!mine.is_open()) {
    std::cerr << "Error trying to open lsn2.1_mine.dat. Terminating.";
    return -1;
  }

  double meanAccumMine = 0, mean2AccumMine = 0, meanMine = 0;
  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < throwsPerBlock; j++) {
      val = rng.AR2_1();
      blockAccumulator += integrand(val)/myPDF(val);
    }

    blockAccumulator /= throwsPerBlock;
    meanAccumMine += blockAccumulator;
    mean2AccumMine += pow(blockAccumulator, 2);
    meanMine = meanAccumMine/(i+1);
    mine << meanMine << " ";
    mine << pow((mean2AccumMine/(i+1) - pow(meanMine, 2))/i, 0.5) << std::endl; //Error nans if i = 0, which indicates an undefined sample mean
  }

  return 0;
}
