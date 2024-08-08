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
  Random gen("../random/seed.in", "../random/primes32001.in", 1);

  //Calculation of mean, variance

  int nBlocks = 100; //Number of blocks for blocking method
  int numbersPerBlock = 1000;

  //Output file for mean, variance
  std::ofstream outputMean("lsn1.1_means.dat");
  if (!outputMean.is_open()) {
    std::cerr << "Failed to open output file lsn1.1_means.dat. Terminating." << std::endl;
    return -1;
  }

  std::ofstream outputVariance("lsn1.1_variances.dat");
  if (!outputVariance.is_open()) {
    std::cerr << "Failed to open output file lsn1.1_variances.dat. Terminating." << std::endl;
    return -1;
  }

  double blockMeanAccumulator = 0;
  double meanAccumulator = 0; //Global mean accumulator
  //Global accumulator for block values of mean squared
  double mean2Accumulator = 0; 
  
  double blockVarianceAccumulator = 0;
  double varianceAccumulator = 0; //Global variance accumulator
  //Global accumulator for block values of variance squared
  double variance2Accumulator = 0;

  double val = 0; //To store the generated random numbers

  //Loop over blocks
  for (int i=0; i<nBlocks; i++) {
    // Reset values of block accumulators
    blockMeanAccumulator = 0;
    blockVarianceAccumulator = 0;

    //Generate random numbers within each block
    for (int j=0; j<numbersPerBlock; j++) {
      val = gen.Rannyu();
      blockMeanAccumulator += val;
      blockVarianceAccumulator += pow(val-0.5, 2);
    }

    // Compute mean and error for the current block, write it to file, update global accumulators
    blockMeanAccumulator /= numbersPerBlock;

    meanAccumulator += blockMeanAccumulator;
    mean2Accumulator += pow(blockMeanAccumulator, 2);
    outputMean << meanAccumulator/(i+1) - 0.5 << " ";
    //Error will nan if i=0, which I'm fine with (it indicates an undefined sample variance and is handled wothout problems by C++ and Python for data analysis).
    outputMean << sqrt((mean2Accumulator/(i+1) - 0.5 - pow(meanAccumulator/(i+1) - 0.5, 2))/i);
    outputMean << std::endl;

    // Compute variance and error for the current block, write it to file, update global accumulators
    blockVarianceAccumulator /= numbersPerBlock;

    varianceAccumulator += blockVarianceAccumulator;
    variance2Accumulator += pow(blockVarianceAccumulator, 2);
    outputVariance << varianceAccumulator/(i+1) << " ";
    //Error will nan if i=0, which I'm fine with (it indicates an undefined sample variance and is handled wothout problems by C++ and Python for data analysis).
    outputVariance << sqrt((variance2Accumulator/(i+1) - pow(varianceAccumulator/(i+1), 2))/i);
    outputVariance << std::endl;
  }
  // Not striclty necessary as the destructor should take care of this when the stream objects come out of scope, but good practice nevertheless.
  outputMean.close();
  outputVariance.close();

  //Chi^2 calculation

  std::ofstream outputChi2("lsn1.1_chi2s.dat");
  if (!outputChi2.is_open()) {
    std::cerr << "Failed to open output file lsn1.1_chi2.dat. Terminating." << std::endl;
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
