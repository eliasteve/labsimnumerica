#include <iostream>
#include <fstream>
#include "../random/random.h"

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);
  
  // Testing the new methods
  int nSamples = 100000;
  std::ofstream CLOut("test_CL.dat"), expOut("test_exp.dat");
  if (CLOut.is_open() && expOut.is_open()){
    for (int i = 0; i < nSamples; i++) {
      CLOut << rng.CauchyLorentz(0.5, 1.7) << std::endl;
      expOut << rng.Exp(0.3) << std::endl;
    }
    CLOut.close();
    expOut.close();
  } else {
    std::cerr << "I/O error when testing, program continues." << std::endl;
  }

  // Die throwing!
  const int valuesOfN = 4;
  int N[valuesOfN] = {1, 2, 10, 100}, realizations = 10000;
  double accumulatorStd, accumulatorExp, accumulatorCauchy;
  std::ofstream stdDie("std_die.dat"), expDie("exp_die.dat"), cauchyDie("cauchy_die.dat");
  if(stdDie.is_open() && expDie.is_open() && cauchyDie.is_open()) {
    for (int i = 0; i < realizations; i++) {
      for (int j = 0; j < valuesOfN; j++) {
        accumulatorStd = 0, accumulatorExp = 0, accumulatorCauchy = 0;

        for (int k = 0; k < N[j]; k++) {
          accumulatorStd += rng.Rannyu();
          accumulatorExp += rng.Exp(1);
          accumulatorCauchy += rng.CauchyLorentz(0, 1);
        }

        accumulatorStd /= N[j];
        accumulatorExp /= N[j];
        accumulatorCauchy /= N[j];

        stdDie << accumulatorStd << " ";
        expDie << accumulatorExp << " ";
        cauchyDie << accumulatorCauchy << " ";
      }
      stdDie << std::endl;
      expDie << std::endl;
      cauchyDie << std::endl;
    }
  } else {
    std::cerr << "I/O error when throwing dice, program terminates." << std::endl;
  }


  return 0;
}
