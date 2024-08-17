#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"
#include "../lib_NSL/misc.h"
#include "../lib_NSL/metropolis.h"
#include "../lib_NSL/point.h"


//Exponential prob distribution for testing
/*
double prob_distr(double x) {
  if (x < 0) return 0;
  else return exp(-x);
}
*/

//Probability density at a point p for the ground state
double psiSqGS(point p);

//Probability density at a point p for the orbital n=2, l=1, m=0
double psiSq2p(point p);

//Funtion to sample from the probability distribution of a given orbital
//and calculate the mean value of the radius.
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
);

int main() {
  int burn_in_steps = 1000; //Burn-in (equilibration steps for Metropolis)
  int nBlocks = 100; //Number of blocks (for data blocking)
  int ptsPerBlock = 1000; //Points per block

  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  //Testing the Metropolis algorithm by sampling from an exponential distribution
  /*
  int N = 100000; //Points to draw
  int accepted = 0; //Accepted points in the Metropolis sampling
  valAndAccept valA = {4, 0}; //Initial value

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
  */

  // Sampling the orbital distributions

  double ar = 0; //Acceptance rate for Metropolis

  //What happens when we start the sampling very far away from the region
  //where most of the probability resides? (And we also don't equilibrate?)

  ar = sampleAndCalculate(
    rng,
    0, //No burn-in or we won't be able to observe the effects of
       //setting the initial value
    nBlocks,
    ptsPerBlock,
    {{100, 0, 0}, 0},
    "out_GS_faraway.dat",
    "out_r_GS_faraway.dat",
    &psiSqGS,
    //Need to static cast this, otherwise compiler can't tell the
    //difference between overloaded versions of this finction.
    static_cast<point(*)(point, double, Random&)> (&propose_unif),
    1.25 //Empirically chosen so acceptance ratio is about half
  );

  std::cout << "Acceptance ratio for GS when starting far away: ";
  std::cout << ar << std::endl;

  //Draw moves from a uniform distribution

  pointAndAccept ptA = {{1, 0, 0}, 0}; //Initial point for all calculations

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_GS.dat",
    "out_r_GS.dat",
    &psiSqGS,
    //Need to static cast this, otherwise compiler can't tell the
    //difference between overloaded versions of this finction.
    static_cast<point(*)(point, double, Random&)> (&propose_unif),
    1.25 //Empirically chosen so acceptance ratio is about half
  );
  std::cout << "Acceptance ratio for GS (uniform): " << ar << std::endl;

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_2p.dat",
    "out_r_2p.dat",
    &psiSq2p,
    //Need to static cast this, otherwise compiler can't tell the
    //difference between this function and an overloaded version of it.
    static_cast<point(*)(point, double, Random&)> (&propose_unif),
    2.7 //Empirically chosen so acceptance ratio is about half
  );
  std::cout << "Acceptance ratio for 2p (uniform): " << ar << std::endl;

  //Draw moves from a multivariate normal

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_GS_gauss.dat",
    "out_r_GS_gauss.dat",
    &psiSqGS,
    propose_gauss,
    0.75 //Empirically chosen so acceptance ratio is about half
  );
  std::cout << "Acceptance ratio for GS (normal): " << ar << std::endl;

  ar = sampleAndCalculate(
    rng,
    burn_in_steps,
    nBlocks,
    ptsPerBlock,
    ptA,
    "out_2p_gauss.dat",
    "out_r_2p_gauss.dat",
    &psiSq2p,
    propose_gauss,
    1.7 //Empirically chosen so acceptance ratio is about half
  );
  std::cout << "Acceptance ratio for 2p (normal): " << ar << std::endl;

  return 0;
}

//Funtion to sample from the probability distribution of a given orbital
//and calculate the mean value of the radius.
double sampleAndCalculate(
  Random &rng, //Random number generator
  int burn_in_steps, //Burn-in (equilibration) steps
  int nBlocks, //Number of blocks (for data blocking)
  int ptsPerBlock, //Point per block
  pointAndAccept ptA, //Initial point for Metropolis algorithm
  std::string filenamePts, //Name of the file containing the
                           //coordinates of the sampled points
  std::string filenameBlkVals, //Name of the file comtaining the
                               //progressive values of the radius
  std::function<double(point)> prob_distr, //Probability distribution to be
                                           //sampled via Metropolis
  std::function<point(point, double, Random&)> propose, //Prob distribution
                                                        //for move
                                                        //proposal
  double cubeHalfSide //Width for the move proposal distribution
) {

  int accepted = 0; //Number of accepted points
  double blockAccumulator = 0; //Accumulator for block values of radius
  double meanAccumulator = 0; //Global accumulator for the mean of radius
  double mean2Accumulator = 0; //Global accumulator for mean^2 of radius

  for (int i=0; i<burn_in_steps; i++) {
    ptA = gen_next_point(ptA.p, prob_distr, propose, rng, cubeHalfSide);
  }
  
  //Output file for the coordinates of the sampled points
  std::ofstream out(filenamePts);
  if (!out.is_open()) {
    std::cerr << "I/O error opening file" + filenamePts + ". Program terminates." << std::endl;
    return -1;
  }

  //Output file for the progressive values of the radius
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

  return double(accepted)/(nBlocks*ptsPerBlock); //Return the acceptance
                                                 //ratio
}

//Probability density at a point p for the ground state
double psiSqGS(point p) {
  double r = sqrt(norm2(p));
  return pow(exp(-r)/sqrt(M_PI), 2);
}

//Probability density at a point p for the orbital n=2, l=1, m=0
double psiSq2p(point p) {
  double r = sqrt(norm2(p));
  return pow(exp(-r/2.)*p.z*sqrt(2./M_PI)/8., 2);
}
