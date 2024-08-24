#include <iostream>
#include "../random/random.h"
#include "../lib_NSL/metropolis.h"
#include <stdlib.h> //EXIT_FAILURE
#include "mcmc.h"


int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  double mu = 1, sigma = 0.5;
  meanErrorAccept val = computeEnergy(
    rng,
    1000,
    1500,
    100,
    0,
    squareAbsWf,
    static_cast<double(*)(double, double, Random&)> (&propose_unif),
    localEnergy,
    1,
    mu,
    sigma,
    ""
  );

  std::cout << "Energy for mu = " << mu << ", sigma = " << sigma << " is " << val.mean << " +- " << val.error << " with acceptance ratio " << val.acceptanceRatio << std::endl;

  val = computeEnergy(
    rng,
    1000,
    1500,
    100,
    0,
    squareAbsWf,
    static_cast<double(*)(double, double, Random&)> (&propose_unif),
    stableLocalEnergy,
    1,
    mu,
    sigma,
    "stable_energy_test.dat"
  );

  std::cout << "STABLE energy for mu = " << mu << ", sigma = " << sigma << " is " << val.mean << " +- " << val.error << " with acceptance ratio " << val.acceptanceRatio << std::endl;

  double xs[6] = {0, 2, 1, -0.5, 0.234, -3};
  for (int i = 0; i < 6; i++) {
    std::cout << localEnergy(xs[i], mu, sigma) << " " << stableLocalEnergy(xs[i], mu, sigma) << std::endl;
  }

  return 0;
}
