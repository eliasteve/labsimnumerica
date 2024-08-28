#include "mcmc.h"

//Adjusts the proposal distribution width based on the acceptance ratio
//The adjustment law was chosen empirically
double adjustProposalWidth(double acceptanceRatio, double width) {
  if (acceptanceRatio < 0.5) {
    width *= 0.8;
  } else if (acceptanceRatio > 0.55) {
    width *= 1.1;
  }
  return width;
}


//Computes the energy (mean value of the Hamiltonian) of a wavefunction via
//Monte Carlo using Metropolis.
meanErrorAccept computeEnergy(
  Random &rng, //Random number generator
  int burn_in_steps, //Number of burn-in (equilibration) steps
  int nBlocks, //Number of blocks for the blocking method
  int ptsPerBlock, //Number of sampled points per block
  double initial_point, //Initial point for the Metropolis algorithm
  //Probability distribution function (\psi^2)
  std::function<double(double, double, double)> prob_distr,
  //Distribution for move proposal
  std::function<double(double, double, Random&)> propose,
  //Local energy function (\psi^-1*H\psi)
  std::function<double(double, double, double)> localEnergy,
  //Initial width of the proposal distribution (will be dynamically
  //adjusted during the equilibration phase)
  double width,
  double mu, //Parameter mu (e.g., center of the wavefunction)
  double sigma, //Parameter sigma (e.g., spread of the wavefunction)
  //Filename for saving points sampled from the distribution, empty string
  //disables saving
  std::string filenameDistro,
  //Filename for saving the progressive values of the energy as a function
  //of the block: empty string disables saving
  std::string filenameEnergy
) {

  int accepted = 0; //Counter for accepted moves
  int flipCount = 0; //Counter for moves before a flip move
  int movesBeforeFlip = 30; //How many "standard" moves before a flip
  int periodForWidthAdj = 20; //Steps between width adjustments
  double blockAccumulator = 0; //Accumulator for block value
  double meanAccumulator = 0; //Global accumulator for mean
  double mean2Accumulator = 0; //Global accumulator for energy
  double blockMean = 0; //Block mean
  double progressiveMean = 0; //Progressive mean
  //Stores value and acceptance status, initialize it to the initial value
  valAndAccept v = {initial_point, 0};
  meanErrorAccept res = {}; //Stores result (mean, error, acceptance ratio)

  std::ofstream outfd;
  if(filenameDistro != "") {
    outfd.open(filenameDistro);
    if(!outfd.is_open()) {
      std::cerr << "I/O error. Program terminates" << std::endl;
      exit(EXIT_FAILURE);
    }
    //File header
    outfd << "mu " << mu << std::endl << "sigma " << sigma << std::endl;
    outfd << "nBlocks " << nBlocks << std::endl;
    outfd << "ptsPerblock " << ptsPerBlock << std::endl;
  }

  std::ofstream outfe;
  if(filenameEnergy != "") {
    outfe.open(filenameEnergy);
    if(!outfe.is_open()) {
      std::cerr << "I/O error. Program terminates" << std::endl;
      exit(EXIT_FAILURE);
    }
    outfe << "Energy Error" << std::endl; //File header
  }
  //Burn-in phase: also adjust width to achieve desired acceptance ratio
  for (int i = 0; i < burn_in_steps / periodForWidthAdj; i++) {
    accepted = 0;
    for (int j = 0; j < periodForWidthAdj; j++) {
      v = gen_next_point(v.val, prob_distr, propose, rng, width, mu, sigma, flipCount, movesBeforeFlip);
      accepted += v.accepted;
    }
    //Adjust proposal width based on acceptance ratio
    width = adjustProposalWidth(double(accepted) / periodForWidthAdj, width);
  }

  accepted = 0; //Reset acceptance counter for main simulation
  //Main Monte Carlo loop: compute energy using the blocking method
  for (int i = 0; i < nBlocks; i++) {
    blockAccumulator = 0;
    for (int j = 0; j < ptsPerBlock; j++) {
      v = gen_next_point(v.val, prob_distr, propose, rng, width, mu, sigma, flipCount, movesBeforeFlip);
      accepted += v.accepted;
      if(filenameDistro != "") {
        outfd << v.val << std::endl;
      }
      blockAccumulator += localEnergy(v.val, mu, sigma);
    }
    blockMean = blockAccumulator / ptsPerBlock;
    meanAccumulator += blockMean;
    mean2Accumulator += pow(blockMean, 2);
    progressiveMean = meanAccumulator/(i+1);
    if(filenameEnergy != "") 
      outfe << progressiveMean << " " << sqrt((mean2Accumulator/(i+1) - pow(progressiveMean, 2))/i) << std::endl;
  }

  //Final calculations: mean energy, statistical error, and acceptance ratio
  res.mean = meanAccumulator / nBlocks;
  res.error = sqrt((mean2Accumulator / nBlocks - pow(res.mean, 2)) / (nBlocks - 1));
  res.acceptanceRatio = double(accepted) / (nBlocks * ptsPerBlock);

  if (filenameDistro != "") outfd.close();
  if (filenameEnergy!= "") outfe.close();

  return res; //Return result (mean energy, error, and acceptance ratio)
}

//Computes the squared absolute value of the wavefunction
double squareAbsWf(
  double x, //Point to calculate the wavefunction in
  //Wavefunction parameters
  double mu,
  double sigma
){
  return pow(exp(-0.5*pow((x-mu)/sigma, 2)) + exp(-0.5*pow((x+mu)/sigma, 2)), 2);
}

//Computes the local energy (\psi^-1*H\psi) for the given wavefunction
double localEnergy(
  double x, //Point to calculate the local energy in
  //Wavefunction parameters
  double mu,
  double sigma
){
  double argminus = pow((x-mu)/sigma, 2), argplus = pow((x+mu)/sigma, 2);
  double laplpsi = ( exp(-0.5*argminus)*(argminus-1) + exp(-0.5*argplus)*(argplus-1) ) / pow(sigma, 2);
  double psi = exp(-0.5*argminus) + exp(-0.5*argplus);
  return - laplpsi / (2*psi) + pow(x, 4) - 5./2.*pow(x, 2);
}

//A more numerically stable version of the local energy computation,
//uses the logsumexp trick (see Jupyter notebook)
double stableLocalEnergy(
  double x, //Point to calculate the wavefunction in
  //Wavefunction parameters
  double mu,
  double sigma
){
  double argminus = pow((x-mu)/sigma, 2), argplus = pow((x+mu)/sigma, 2);
  double argmax = std::max(-0.5*argplus, -0.5*argminus);
  double sumexp = exp(-0.5*argplus - argmax) + exp(-0.5*argminus - argmax);
  double logprobplus = -0.5*argplus  - argmax - log(sumexp);
  double logprobminus = -0.5*argminus - argmax - log(sumexp);
  return -0.5/pow(sigma, 2) * ( exp(logprobplus) * (argplus-1) + exp(logprobminus) * (argminus-1) ) + pow(x, 4) - 5./2.*pow(x, 2);
}


