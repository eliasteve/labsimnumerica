#include <iostream>
#include "../random/random.h"
#include "../lib_NSL/metropolis.h"
#include "mcmc.h"
#include <numeric> // std::accumulate

//Struct to hold information about the simulated annealing step
struct SAStep {
  double beta; //"Inverse temperature"
  int n; //Step number
  int nBlocks; //How many blocks we use to compute the energy
};

//Implements the Metropolis step to generate a new pair of values for the
//parameters
void generateNewParameters(
  double oldMu,
  double oldSigma,
  meanErrorAccept oldEnergy,
  double& mu,
  double& sigma,
  meanErrorAccept& energy,
  int& accepted,
  int& accepted2,
  double parametersStepWidth,
  double energyStepWidth,
  double beta,
  int nBlocks,
  Random& rng
);

int main() {
  //Random number generator
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  meanErrorAccept val; //Stores the result of an MCMC energy computation

  //Initial width for the proposal distribution for the MCMC energy 
  //calculation (will be adjusted during equilibration to get an about 50%
  //acceptance ratio)
  double energyStepWidth = 2;
  //Initial width for the proposal distribution for parameter calculation.
  //Is adjusted dynamically during sampling, since as we're sampling from 
  //different distribution as beta changes we do not care as much about
  //detailed balanca.
  double parametersStepWidth = 0.5;
  //Acceptance fraction for parameters calculation
  double acceptanceFraction = 0;
  //Minimum value for beta
  double beta_min = 5;
  //maximum vslue for beta
  double beta_max = 110;
  //Exponent in the scaling law for beta,
  //beta_n = K(n-1)^exponent + beta_min
  double exponent = 1.2;
  //Steps in the annealing process
  int nSteps = 100;
  //Generated values of the parameters per annealing step
  int pointsPerStep = 100;
  //Factor in the scaling law for beta. Fixed so beta_nSteps = beta_max
  double K = (beta_max - beta_min)/pow(nSteps-1, exponent);
  //Acceoted points in every simulated annealing step
  int accepted = 0;
  //Accepted points between width adjustments
  int accepted2 = 0;
  //Counts the values generated between width adjustments
  int counterForWidthAdjustment = 0;
  //Perform a width adjustment after this many moves
  int widthAdjustmentPeriod = 20;
  //Output files for the values of the parameters, the energies and beta
  std::ofstream outfParameters, outfEnergies, outfBetas;

  //Construct the the annealing schedule

  //One per steps plus a fistitious 1st step which fixes the initial
  //values for the parameters
  SAStep SASchedule[nSteps+1];
  //Vectors which will hold the generated parameters and energies in
  //each step.
  std::vector<double> mus = {0.}, sigmas = {1.};
  std::vector<meanErrorAccept> energy = {val};
  //Useful for resizing the vectors in the 1st (actual) step
  SASchedule[0].n = 1;
  SASchedule[0].beta = beta_min; //Could also leave this uninitialized

  //Initialize the other steps. The numbers of blocks per steps changes
  //dynamically in order to compute the energy with increasingly smaller
  //uncertainties as the variability in the parameters decreases. The
  //actual change schedule was chosen empirically
  for(int i = 1; i < nSteps+1; i++) {
    SASchedule[i].n = pointsPerStep;
    SASchedule[i].beta = K*pow(i-1, exponent) + beta_min;
    if(i < nSteps-75) SASchedule[i].nBlocks = 30;
    else if(i < nSteps-60) SASchedule[i].nBlocks = 60;
    else if(i < nSteps-35) SASchedule[i].nBlocks = 120;
    else SASchedule[i].nBlocks = 150;
  }

  outfBetas.open("optimizationData/betas.dat");
  outfBetas << "Step Beta" << std::endl;

  //Annealing process

  for(int step = 1; step < nSteps+1; step++) {
    //Initiliaze first point in the step
    mus[0] = mus[SASchedule[step-1].n-1];
    sigmas[0] = sigmas[SASchedule[step-1].n-1];
    energy[0] = energy[SASchedule[step-1].n-1];

    //If the number of steps has increased since the previous iteration, resize the vectors containing the optimization data to match.
    if (SASchedule[step].n > SASchedule[step-1].n) {
      mus.resize(SASchedule[step].n, 0.);
      sigmas.resize(SASchedule[step].n, 0.);
      energy.resize(SASchedule[step].n);
    }

    accepted2 = 0;
    std::cerr << "Beginning step " << step << " (beta = " << SASchedule[step].beta << ", nBlocks " << SASchedule[step].nBlocks << ")" << std::endl;

    for(int i = 1; i < SASchedule[step].n; i++) {
      generateNewParameters(mus[i-1], sigmas[i-1], energy[i-1], mus[i], sigmas[i], energy[i], accepted, accepted2, parametersStepWidth, energyStepWidth, SASchedule[step].beta, SASchedule[step].nBlocks, rng);

      //Adjust the width of the distribution for parameter move proposal
      //dynamically in order to obtain an acceptance ratio close to 50%
      counterForWidthAdjustment++;
      if (counterForWidthAdjustment%widthAdjustmentPeriod == 0) {
        acceptanceFraction = double(accepted)/widthAdjustmentPeriod;
        parametersStepWidth = adjustProposalWidth(
          acceptanceFraction,
          parametersStepWidth
        );
        accepted = 0;
        counterForWidthAdjustment = 0;
      }
    }

    //Print values to console
    std::cout << "Acceptance: " << double(accepted2)/(SASchedule[step].n-1) << std::endl << "Last value of energy: " << energy[SASchedule[step].n-1].mean << " ± " << energy[SASchedule[step].n-1].error << std::endl << "Parameters: mu = " << mus[SASchedule[step].n-1] << ", sigma = " << sigmas[SASchedule[step].n-1] << std::endl;

    //Print values to file
    outfParameters.open("optimizationData/parameters_step" + std::to_string(step) + ".dat");
    outfEnergies.open("optimizationData/energies_step" + std::to_string(step) + ".dat");
    outfEnergies << "Energy Error Acceptance" << std::endl;
    outfParameters << "mu sigma" << std::endl;
    for (int i = 0; i < SASchedule[step].n; i++) {
      outfEnergies << energy[i].mean << " " << energy[i].error << " " << energy[i].acceptanceRatio << std::endl;
      outfParameters << mus[i] << " " << sigmas[i] << std::endl;
    }
    outfBetas << step << " " << SASchedule[step].beta << std::endl;
    outfEnergies.close();
    outfParameters.close();
  }
  outfBetas.close();

  //Final calculation of the energy, to write the sampled point of the
  //wavefunction to file
  val = computeEnergy(
    rng,
    300,
    150,
    100,
    *(mus.rbegin()),
    squareAbsWf,
    static_cast<double(*)(double, double, Random&)> (&propose_unif),
    stableLocalEnergy,
    10*(*(sigmas.rbegin())),
    *(mus.rbegin()),
    *(sigmas.rbegin()),
    "final_energy_distro.dat"
  );
  std::cout << "Parameters: mu " << *(mus.rbegin()) << ", sigma " << *(sigmas.rbegin()) << std::endl;
  std::cout << "energy = " << val.mean << " ± " << val.error << ", acc = " << val.acceptanceRatio << std::endl << std::endl;

  return 0;
}

//Implements the Metropolis step to generate a new pair of values for the
//parameters
void generateNewParameters(
  //Old values for mu, sigma
  double oldMu,
  double oldSigma,
  //Energy for the olv values of the parameters
  meanErrorAccept oldEnergy,
  //Variables in which to write the new values of the parameters
  double& mu,
  double& sigma,
  //Variable in which to write the new value of the energy
  meanErrorAccept& energy,
  //counts how many values are accepted in the period between width 
  //adjustments
  int& accepted,
  //Counts how many values are accepted in the current annealing step
  int& accepted2,
  //(Initial) width for the proposal distribution for parameter moves
  double parametersStepWidth,
  //(Initial) width for the proposal distribution for moves in the
  //calculation of energy
  double energyStepWidth,
  double beta, //"Inverse temperature"
  int nBlocks, //How many blocks to use for the calculation of the energy
  Random& rng //Random number generator
) {
  //Endpoints of interval from which to draw proposal moves with a uniform
  //distribution
  double lower = 0, upper = 0;
  //Probability to determine wether to accept a proposed value
  double prob = 0;
  
  //Propose next point from a uniform distribution. Since both mu and sigma
  //should be positive, if the lower endpoint of the interval is negative,
  //shift the interval so it contains only nonnegative values.
  lower = oldMu-parametersStepWidth;
  upper = oldMu+parametersStepWidth;
  if(lower < 0) {
    upper -= lower;
    lower = 0;
  }
  mu = rng.Rannyu(lower, upper);
  lower = oldSigma-parametersStepWidth;
  upper = oldSigma+parametersStepWidth;
  if(lower < 0) {
    upper -= lower;
    lower = 0;
  }
  sigma = rng.Rannyu(lower, upper);

  //Compute the energy of the point (use nBlocks blocks of 100 steps each
  //and 300 equilibration steps, in which se also determine the width of
  //the proposal distribution)
  energy = computeEnergy(
    rng,
    300,
    nBlocks,
    100,
    mu,
    squareAbsWf,
    static_cast<double(*)(double, double, Random&)> (&propose_unif),
    stableLocalEnergy,
    sigma,
    mu,
    sigma,
    ""
  );

  //Accept new point?
  prob = exp(-beta * (energy.mean - oldEnergy.mean));
  if (rng.Rannyu() < prob) {
    accepted++;
    accepted2++;
  }
  else {
    mu = oldMu;
    sigma = oldSigma;
    energy = oldEnergy;
  }
}
