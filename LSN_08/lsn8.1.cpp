#include <iostream>
#include "../random/random.h"
#include "../lib_NSL/metropolis.h"
#include <stdlib.h> //EXIT_FAILURE
#include "mcmc.h"


int main() {
  //Computes the energy for particular values of mu, sigma as a test

  //Random number generator
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  //output file to store the data relative to the test of the energy
  std::ofstream out("energy_test.dat");
  if(!out.is_open()) {
    std::cerr << "I/O error!" << std::endl;
    exit(EXIT_FAILURE);
  }
  out << "mu sigma energy_naive error_naive energy_stable error_stable" << std::endl;

  double mu = 1, sigma = 0.5; //Wavefunction and local energy parameters
  out << mu << " " << sigma << " ";

  meanErrorAccept val = computeEnergy( //Result of energy computation
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
    "",
    ""
  );

  std::cout << "NAIVE energy for mu = " << mu << ", sigma = " << sigma << " is " << val.mean << " +- " << val.error << " with acceptance ratio " << val.acceptanceRatio << std::endl;
  out << val.mean << " " << val.error << " ";

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
    "",
    ""
  );

  std::cout << "STABLE energy for mu = " << mu << ", sigma = " << sigma << " is " << val.mean << " +- " << val.error << " with acceptance ratio " << val.acceptanceRatio << std::endl;
  out << val.mean << " " << val.error << " ";

  out.close();
  return 0;
}
