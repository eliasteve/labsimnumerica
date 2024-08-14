#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"

double getS(double deltaT, double sqrtDeltaT, Random &rng, double S0, double mu, double sigma);

void computeUpdateMeanAndError(int trialsPerBlock, int progressiveIndex, double blockAccumulator, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out);

int main() {
  double S0 = 100; //Asset price at time 0
  double T = 1; //Delivery time
  double K = 100; //Strike price
  double r = 0.1; //Risk-free interest rate
  double sigma = 0.25; //Volatility

  double sqrtT = sqrt(T); //Cache the sqrt so I don't have to calculate it every time I call S

  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  std::ofstream call_direct("lsn3.1_call_direct.dat");
  if (!call_direct.is_open()) {
    std::cerr << "Unable to open output file lsn3.1_call_direct.dat. The program terminates." << std::endl;
    return -1;
  }

  std::ofstream put_direct("lsn3.1_put_direct.dat");
  if (!put_direct.is_open()) {
    std::cerr << "Unable to open output file lsn3.1_put_direct.dat. The program terminates." << std::endl;
    return -1;
  }

  std::ofstream call_discretized("lsn3.1_call_discretized.dat");
  if (!call_discretized.is_open()) {
    std::cerr << "Unable to open output file lsn3.1_call_discretized.dat. The program terminates." << std::endl;
    return -1;
  }

  std::ofstream put_discretized("lsn3.1_put_discretized.dat");
  if (!put_discretized.is_open()) {
    std::cerr << "Unable to open output file lsn3.1_put_discretized.dat. The program terminates." << std::endl;
    return -1;
  }


  int nBlocks = 100, trialsPerBlock = 1000, nTimeIntervals = 100;
  double deltaT = T/nTimeIntervals, sqrtDeltaT = sqrt(deltaT); //Cache the sqrt
  double CBlockAccumulator = 0, CMeanAccumulator = 0, CMean2Accumulator = 0;
  double PBlockAccumulator = 0, PMeanAccumulator = 0, PMean2Accumulator = 0;
  double S = 0;

  for (int i = 0; i < nBlocks; i++) {
    CBlockAccumulator = 0;
    PBlockAccumulator = 0;
    for (int j = 0; j < trialsPerBlock; j++) {
      S = getS(T, sqrtT, rng, S0, r, sigma);
      CBlockAccumulator += exp(-r*T)*std::max(0., S-K);
      PBlockAccumulator += exp(-r*T)*std::max(0., K-S);
    }

    computeUpdateMeanAndError(trialsPerBlock, i, CBlockAccumulator, CMeanAccumulator, CMean2Accumulator, call_direct);
    computeUpdateMeanAndError(trialsPerBlock, i, PBlockAccumulator, PMeanAccumulator, PMean2Accumulator, put_direct);

  }

  call_direct.close();
  put_direct.close();


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
      CBlockAccumulator += exp(-r*T)*std::max(0., S-K);
      PBlockAccumulator += exp(-r*T)*std::max(0., K-S);
    }

    computeUpdateMeanAndError(trialsPerBlock, i, CBlockAccumulator, CMeanAccumulator, CMean2Accumulator, call_discretized);
    computeUpdateMeanAndError(trialsPerBlock, i, PBlockAccumulator, PMeanAccumulator, PMean2Accumulator, put_discretized);
  }

  call_discretized.close();
  put_discretized.close();

  return 0;
}

double getS (
  double deltaT,
  double sqrtDeltaT, //deltaT is not going to change across calls of the function, cache its sqrt for efficiency.
  Random &rng,
  double S0,
  double mu,
  double sigma
) {
  return S0 * exp( (mu-pow(sigma, 2)/2.) * deltaT + sigma*rng.Gauss(0, sqrtDeltaT));
}

void computeUpdateMeanAndError(int trialsPerBlock, int progressiveIndex, double blockAccumulator, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out) {
  double blockMean, mean;

  blockMean = blockAccumulator/trialsPerBlock;
  meanAccumulator += blockMean;
  mean2Accumulator += pow(blockMean, 2);
  mean = meanAccumulator/(progressiveIndex+1);
  out << mean << " " << sqrt((mean2Accumulator/(progressiveIndex+1) - pow(mean, 2))/progressiveIndex) << std::endl;
}
