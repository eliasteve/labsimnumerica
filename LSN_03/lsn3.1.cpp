#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"
#include "../lib_NSL/misc.h"

//Computes the spot price by sampling a Geometric Brownian Motion
double getS(double deltaT, double sqrtDeltaT, Random &rng, double S0, double mu, double sigma);

//Function to accumulate block values for call, put using only one
//comperison (instead of calling the max function twice). Not used because
//upon testing the two versions run in about the same time, and this is
//less legible.
void accumulateBlockValues(double T, double r, double S, double K, double& CBlockAccumulator, double& PBlockAccumulator);

int main() {
  double S0 = 100; //Asset price at time 0
  double T = 1; //Delivery time
  double K = 100; //Strike price
  double r = 0.1; //Risk-free interest rate
  double sigma = 0.25; //Volatility

  double sqrtT = sqrt(T); //sqrt of T. Cache it so I don't have to calculate it every time I call getS

  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  int nBlocks = 100; //number of blocks
  int trialsPerBlock = 1000; //computations of option price per block
  int nTimeIntervals = 100; //No. of time intervals for discretized
                            //computation
  double deltaT = T/nTimeIntervals; //Duration of each time interval
  double sqrtDeltaT = sqrt(deltaT); //Cached sqrt of deltaT
  //Block and global accumulators for call value
  double CBlockAccumulator = 0, CMeanAccumulator = 0, CMean2Accumulator = 0;
  //Block and global accumulators for put value
  double PBlockAccumulator = 0, PMeanAccumulator = 0, PMean2Accumulator = 0;
  double S = 0; //Spot price

  //Direct calculation

  std::ofstream callOut("call_direct.dat");
  if (!callOut.is_open()) {
    std::cerr << "Unable to open output file call_direct.dat. The program terminates." << std::endl;
    return -1;
  }

  std::ofstream putOut("put_direct.dat");
  if (!putOut.is_open()) {
    std::cerr << "Unable to open output file put_direct.dat. The program terminates." << std::endl;
    return -1;
  }

  for (int i = 0; i < nBlocks; i++) {
    CBlockAccumulator = 0;
    PBlockAccumulator = 0;
    for (int j = 0; j < trialsPerBlock; j++) {
      S = getS(T, sqrtT, rng, S0, r, sigma);
      //accumulateBlockValues(T, r, S, K, CBlockAccumulator, PBlockAccumulator);
      CBlockAccumulator += exp(-r*T)*std::max(0., S-K);
      PBlockAccumulator += exp(-r*T)*std::max(0., K-S);
    }

    CBlockAccumulator /= trialsPerBlock;
    PBlockAccumulator /= trialsPerBlock;

    computeUpdateMeanAndError(i, CBlockAccumulator, CMeanAccumulator, CMean2Accumulator, callOut);
    computeUpdateMeanAndError(i, PBlockAccumulator, PMeanAccumulator, PMean2Accumulator, putOut);

  }

  callOut.close();
  putOut.close();

  //Discretized calculation

  callOut.open("call_discretized.dat");
  if (!callOut.is_open()) {
    std::cerr << "Unable to open output file call_discretized.dat. The program terminates." << std::endl;
    return -1;
  }

  putOut.open("put_discretized.dat");
  if (!putOut.is_open()) {
    std::cerr << "Unable to open output file put_discretized.dat. The program terminates." << std::endl;
    return -1;
  }

  //Reset global accumulators
  CMeanAccumulator = 0;
  CMean2Accumulator = 0;
  PMeanAccumulator = 0;
  PMean2Accumulator = 0;

  for (int i = 0; i < nBlocks; i++) {
    CBlockAccumulator = 0;
    PBlockAccumulator = 0;
    for (int j = 0; j < trialsPerBlock; j++) {
      S = S0;
      for (int k = 0; k < nTimeIntervals; k++) {
        S = getS(deltaT, sqrtDeltaT, rng, S, r, sigma);
      }
      //accumulateBlockValues(T, r, S, K, CBlockAccumulator, PBlockAccumulator);
      CBlockAccumulator += exp(-r*T)*std::max(0., S-K);
      PBlockAccumulator += exp(-r*T)*std::max(0., K-S);
    }
    
    CBlockAccumulator /= trialsPerBlock;
    PBlockAccumulator /= trialsPerBlock;

    computeUpdateMeanAndError(i, CBlockAccumulator, CMeanAccumulator, CMean2Accumulator, callOut);
    computeUpdateMeanAndError(i, PBlockAccumulator, PMeanAccumulator, PMean2Accumulator, putOut);
  }

  callOut.close();
  putOut.close();

  return 0;
}

//Computes the spot price by sampling a Geometric Brownian Motion
double getS (
  double deltaT, //Time interval
  double sqrtDeltaT, //sqrt of time interval: deltaT is not going to change across calls of the function, cache its sqrt for efficiency.
  Random &rng, //Random number generator
  double S0, //Value of S at initial time
  double mu, //Drift
  double sigma //Diffision constant
) {
  return S0 * exp( (mu-pow(sigma, 2)/2.) * deltaT + sigma*rng.Gauss(0, sqrtDeltaT));
}

//Function to accumulate block values for call, put using only one
//comperison (instead of calling the max function twice). Not used because
//upon testing the two versions run in about the same time, and this is
//less legible.
void accumulateBlockValues(double T, double r, double S, double K, double& CBlockAccumulator, double& PBlockAccumulator) {
  double difference = S-K;
  if (difference > 0) {
    CBlockAccumulator += exp(-r*T)*difference;
  } else {
    PBlockAccumulator -= exp(-r*T)*difference;
  }
  return;
}
