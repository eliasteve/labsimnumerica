#ifndef __MCMC__
#define __MCMC__

#include <iostream>
#include "../random/random.h"
#include "../lib_NSL/metropolis.h"
#include <stdlib.h> //EXIT_FAILURE
#include <fstream>

//Struct to store the results of the energy computation
struct meanErrorAccept {
  double mean; //Energy
  double error; //Statistical error on the energy
  double acceptanceRatio; //Acceptance ratio of the Metropolis algorithm
};

//Adjusts the proposal distribution width based on the acceptance ratio
double adjustProposalWidth(double acceptanceRatio, double width);

//Computes the energy (mean value of the Hamiltonian) of a wavefunction via
//Monte Carlo using Metropolis.
meanErrorAccept computeEnergy(
  Random &rng,
  int burn_in_steps,
  int nBlocks,
  int ptsPerBlock,
  double initial_point,
  std::function<double(double, double, double)> prob_distr,
  std::function<double(double, double, Random&)> propose,
  std::function<double(double, double, double)> localEnergy,
  double width,
  double mu,
  double sigma,
  std::string filename,
  std::string filenameEnergy
);

//Computes the squared absolute value of the wavefunction
double squareAbsWf(double x, double mu, double sigma);

//Computes the local energy (\psi^-1*H\psi) for the given wavefunction
double localEnergy(double x, double mu, double sigma);

//A more numerically stable version of the local energy computation,
//uses the logsumexp trick (see Jupyter notebook)
double stableLocalEnergy(double x, double mu, double sigma);

#endif //__MCMC__
