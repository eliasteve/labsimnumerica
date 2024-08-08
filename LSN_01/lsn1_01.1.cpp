#include <iostream>
#include <fstream>
#include <cmath>
#include "../random/random.h"

void fillWithZeros(double* array, int length){
  for (int i = 0; i < length; i++) {
    array[i] = 0;
  }
}

int main() {
  // Random generator initialization
  Random gen;

  std::ifstream seed("random/seed.in");
  int s[4];

  if (seed.is_open()) {
    seed.ignore(10); // First word will be RANDOMSEED
    for (int i=0; i<4; i++) {
      seed >> s[i]; 
    }
    seed.close();
  }
  else {
    std::cout << "IO ERROR" << std::endl;
    return -1;
  }

  std::ifstream primes("random/primes32001.in");
  int p1=0, p2=0;
  if (primes.is_open()) {
    primes >> p1 >> p2;
    primes.close();
  }
  else {
    std::cout << "IO ERROR" << std::endl;
    return -1;
  }

  std::cout << "Seed digits: " << s[0] << " " << " " << s[1] << " " << " " << s[2] << " " << " " << s[3] << std::endl;
  std::cout << "Primes digits: " << p1 << " " << p2 << std::endl;

  gen.SetRandom(s, p1, p2);
  
  //Calculation of mean, variance

  int nBlocks = 100, numbersPerBlock = 1000, totalNumbers = nBlocks*numbersPerBlock;

  std::ifstream rng_galli("rng_galli.dat");

  std::ofstream outputMean("lsn1_01.1_means.dat");
  if (!outputMean.is_open()) {
    std::cerr << "Failed to open output file lsn1_01.1_means.dat. Terminating." << std::endl;
    return -1;
  }

  std::ofstream outputVariance("lsn1_01.1_variances.dat");
  if (!outputVariance.is_open()) {
    std::cerr << "Failed to open output file lsn1_01.1_variances.dat. Terminating." << std::endl;
    return -1;
  }

  double mean = 0, errorMean = 0;
  double meanAccumulator = 0, blockMeanAccumulator = 0, mean2Accumulator = 0;
  
  double variance = 0, errorVariance = 0;
  double blockVarianceAccumulator = 0, varianceAccumulator = 0, variance2Accumulator = 0;

  double val = 0;

  for (int i=0; i<nBlocks; i++) {
    // Extract numbers for a single block
    blockMeanAccumulator = 0;
    blockVarianceAccumulator = 0;

    for (int j=0; j<numbersPerBlock; j++) {
      //rng_galli >> val;
      val = gen.Rannyu();
      blockMeanAccumulator += val;
      blockVarianceAccumulator += pow(val-0.5, 2);
    }

    // Compute quantities of interest for a single block
    blockMeanAccumulator /= numbersPerBlock;
    blockVarianceAccumulator /= numbersPerBlock;

    meanAccumulator += blockMeanAccumulator;
    mean2Accumulator += pow(blockMeanAccumulator, 2);
    errorMean = mean2Accumulator/(i+1) - pow(meanAccumulator/(i+1), 2);
    errorMean = pow(errorMean/i, 0.5); //Will nan if i=0, which I'm fine with (it indicates an undefined sample variance).
    mean = meanAccumulator/(i+1);
    outputMean << mean << " " << errorMean << std::endl;

    varianceAccumulator += blockVarianceAccumulator;
    variance2Accumulator += pow(blockVarianceAccumulator, 2);
    errorVariance = variance2Accumulator/(i+1) - pow(varianceAccumulator/(i+1), 2);
    errorVariance = pow(errorVariance/i, 0.5); //Will nan if i=0, which I'm fine with (it indicates an undefined sample variance).
    variance = varianceAccumulator/(i+1);
    outputVariance << variance << " " << errorVariance << std::endl;
    }

  // Not striclty necessary as the destructor should take care of this when the stream objects come out of scope, but good practice nevertheless.
  outputMean.close();
  outputVariance.close();

  //Chi^2 calculation

  std::ofstream outputChi2("lsn1_01.1_chi2s.dat");
  if (!outputChi2.is_open()) {
    std::cerr << "Failed to open output file lsn1_01.1_chi2.dat. Terminating." << std::endl;
    return -1;
  }
  int intervalSubdivisions = 100, chi2Trials = 100, valuesPerTrial = 10000;

  double successes[intervalSubdivisions], chi2 = 0;
  double expectedSuccesses = valuesPerTrial/intervalSubdivisions;
  int integer_index = 0;

  std::cout << "Chi-squared:" << std::endl;
  for (int i=0; i < chi2Trials; i++) {
    fillWithZeros(successes, intervalSubdivisions);
    chi2 = 0;

    // Run single trial
    for (int j=0; j < valuesPerTrial; j++) {
      val = gen.Rannyu();
      val *= intervalSubdivisions; // 0 <= val < intervalSubdivisions
      integer_index = int(val); // Drops the decimal part, == i if val was in the i-th interval.
      successes[integer_index]++;
    }

    // Accumulate values for chi2_i
    for (int j=0; j < intervalSubdivisions; j++) {
      chi2 += pow(successes[j] - expectedSuccesses, 2)/expectedSuccesses;
    }
    std::cout << i << " " << chi2 << std::endl;
    outputChi2 << chi2 << std::endl;
  }

  // Not striclty necessary as the destructor should take care of this when the stream objects come out of scope, but good practice nevertheless.
  outputChi2.close();

  return 0;
}
