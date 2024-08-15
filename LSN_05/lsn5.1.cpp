#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"
#include "../lib_NSL/misc.h"
#include "../lib_NSL/metropolis.h"
#include "../lib_NSL/point.h"

/* double prob_distr(double);
double gen_next_point(double, std::function<double(double)>, std::function<double(double, double, Random)>, Random);
double propose_unif(double, double, Random); */

double prob_distr(double x) {
  if (x < 0) return 0;
  else return exp(-x);
}

double psiSqGS(point p) {
  double r = sqrt(norm2(p));
  return pow(exp(-r)/sqrt(M_PI), 2);
}

double psiSq2p(point p) {
  double r = sqrt(norm2(p));
  return pow(exp(-r/2.)*p.z*sqrt(2./M_PI)/8., 2);
}

double sampleAndCalculate(
  Random &rng,
  int burn_in_steps,
  int nBlocks,
  int ptsPerBlock,
  pointAndAccept ptA,
  std::string filenamePts,
  std::string filenameBlkVals,
  std::function<double(point)> prob_distr,
  std::function<point(point, double, Random&)> propose,
  double cubeHalfSide
) {

  int accepted = 0;
  double blockAccumulator = 0, meanAccumulator = 0, mean2Accumulator = 0;

  for (int i=0; i<burn_in_steps; i++) {
    ptA = gen_next_point(ptA.p, prob_distr, propose, rng, cubeHalfSide);
  }
  
  std::ofstream out(filenamePts);
  if (!out.is_open()) {
    std::cerr << "I/O error opening file" + filenamePts + ". Program terminates." << std::endl;
    return -1;
  }

  std::ofstream outBV(filenameBlkVals);
  if (!out.is_open()) {
    std::cerr << "I/O error opening file" + filenameBlkVals + ". Program terminates." << std::endl;
    return -1;
  }

  for (int i=0; i<nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < ptsPerBlock; j++) {
      ptA = gen_next_point(ptA.p, prob_distr, propose, rng, cubeHalfSide);
      accepted += ptA.accepted;
      out << ptA.p.x << " " << ptA.p.y << " " << ptA.p.z << std::endl;

      blockAccumulator += sqrt(norm2(ptA.p));
    }
    computeUpdateMeanAndError(i, blockAccumulator/ptsPerBlock, meanAccumulator, mean2Accumulator, outBV);
  }

  out.close();
  outBV.close();

  return double(accepted)/(nBlocks*ptsPerBlock);
}

int main() {
  int burn_in_steps = 1000, N = 100000, accepted = 0;
  int nBlocks = 100, ptsPerBlock = 1000;

  Random rng("../random/seed.in", "../random/primes32001.in", 2);

  valAndAccept valA = {4, 0};

  for (int i = 0; i<burn_in_steps; i++) {
    valA = gen_next_point(valA.val, &prob_distr, static_cast<double(*)(double, double, Random&)> (&propose_unif), rng, 1.5);
  }

  std::ofstream out("out_test.dat");
  if (!out.is_open()) {
    std::cerr << "I/O error opening file out_test.dat. Program terminates." << std::endl;
    return -1;
  }

  for (int i = 0; i<N; i++) {
    valA = gen_next_point(valA.val, &prob_distr, static_cast<double(*)(double, double, Random&)> (&propose_unif), rng, 1.5);
    accepted += valA.accepted;
    out << valA.val << std::endl;
  }
  std::cout << "Acceptance ratio for test " << double(accepted)/N << std::endl;

  out.close();

  pointAndAccept ptA = {{1, 0, 0}, 0};
  double ar = 0;

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_GS.dat",
    "out_r_GS.dat",
    &psiSqGS,
    static_cast<point(*)(point, double, Random&)> (&propose_unif),
    1.25
  );
  std::cout << "Acceptance ratio for GS: " << ar << std::endl;

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_2p.dat",
    "out_r_2p.dat",
    &psiSq2p,
    static_cast<point(*)(point, double, Random&)> (&propose_unif),
    2.7
  );
  std::cout << "Acceptance ratio for 2p: " << ar << std::endl;

  return 0;
}


