#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"
#include "../lib_NSL/point.h"
#include "../lib_NSL/misc.h"

//
// Function to perform a discrete step in the random walk
point discreteStep(Random& rng, double latticeConst);

// Function to perform a continuous step in the random walk
point continuousStep(Random& rng, double stepLength);

// Function to perform a random walk simulation
void performRW(
  int nBlocks,
  int walkersPerBlock,
  int nSteps,
  int stepLength,
  std::string outFileName,
  std::function<point(Random&, double)> stepFunction,
  Random& rng
);

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  int nBlocks = 100;           // Number of blocks for the Monte Carlo simulation
  int walkersPerBlock = 100;   // Number of walkers per block
  int nSteps = 100;            // Number of steps for each walker
  double a = 2;                // Lattice constant for the discrete random walk (and
                               // length of step for the continuous one)

  // Perform random walk with discrete steps
  performRW(nBlocks, walkersPerBlock, nSteps, a, "discreteRW.dat", discreteStep, rng);

  // Perform random walk with continuous steps
  performRW(nBlocks, walkersPerBlock, nSteps, a, "continuousRW.dat", continuousStep, rng);

  return 0;
}

// Function to perform a discrete step in the random walk
point discreteStep(Random& rng, double latticeConst) {
  //Random numbers to decide what the step will be:
  //0<=y<1: step in the x direction, 1<=y<2: step in the y direction,
  //2<=y<3: step in the z direction
  //x<0.5: step in the positive direction, x>=0.5: step in the negative
  //direction
  double x = rng.Rannyu(), y = rng.Rannyu(0, 3);
  point p = {0, 0, 0}; //Will store the step
  double toAdd = 0; //Value which will be added to one of the coordinates
                    //of p in order to make the step  
  if (x < 0.5) {
    toAdd = -latticeConst;
  } else {
    toAdd = latticeConst;
  }

  if (y < 1) {
    p.x = toAdd;
  } else if (y < 2) {
    p.y = toAdd;
  } else {
    p.z = toAdd;
  }

  return p;
}

// Function to perform a continuous step in the random walk
point continuousStep(Random& rng, double stepLength) {
  point step = {0, 0, 0};
  double norm;

  //Using a simple accept-reject method: sample a point uniformly in the
  //"unit cube", reject it if it's not inside the unit sphere; when you
  //accept a point then rescale il so it's on the surface of the sphere
  do {
    step.x = rng.Rannyu(-1, 1);
    step.y = rng.Rannyu(-1, 1);
    step.z = rng.Rannyu(-1, 1);
    norm = norm2(step);
  } while (norm > 1);

  return step/(sqrt(norm)/stepLength);
}

// Function to perform a random walk simulation
void performRW(
  int nBlocks, //Number of blocks in qhich to divide the simulation
  int walkersPerBlock, //Walks per block
  int nSteps, //Steps for each walk
  int stepLength, //Length of the step
  std::string outFileName, //output file
  std::function<point(Random&, double)> stepFunction, //Function to gene-
                                                      //rate the step
  Random& rng //Random number generator
) {
  //Output file for average of norm of posizions
  std::ofstream outf(outFileName);
  //Output file for vector posizions, for debugging
  //std::ofstream posOut("positions.dat");
  if(!outf.is_open()) {
    std::cerr << "I/O error opening " << outFileName << ", program continues." << std::endl;
    return;
  }

  //Will store the positions of the walkers.
  //Allocate this dynamically as a best practice, as there is no guarantee
  //that the arguments of a function will be known at compile time (even
  //though in practice I'm the only one who will ever call this function,
  //and I _know_ that its arguments are known at compile time)
  point** positions = new point*[nBlocks];
  for (int i = 0; i < nBlocks; i++) {
    positions[i] = new point[walkersPerBlock];
    for (int j = 0; j < walkersPerBlock; j++) {
      positions[i][j] = {0, 0, 0};
    }
  }

  double blockNorm2Accumulator = 0;        // Accumulator for the norm in a block
  double meanAccumulator = 0;              // Accumulator for the mean value
  double mean2Accumulator = 0;             // Accumulator for the square of the mean value
  double blockMean = 0;                    // Mean value of the norm per block
  double stepMean = 0;                     // Mean value of the norm at each step

  for (int i = 0; i < nSteps; i++) {
    meanAccumulator = 0;
    mean2Accumulator = 0;
    for (int k = 0; k < nBlocks; k++) {
      blockNorm2Accumulator = 0;
      for (int j = 0; j < walkersPerBlock; j++) {
        positions[k][j] = positions[k][j] + stepFunction(rng, stepLength);
        //posOut << positions[k][j] << " "; //debug
        blockNorm2Accumulator += norm2(positions[k][j]);
      }
      blockMean = sqrt(blockNorm2Accumulator/walkersPerBlock);
      meanAccumulator += blockMean;
      mean2Accumulator += pow(blockMean, 2);
    }
    //posOut << std::endl; //debug
    stepMean = meanAccumulator/nBlocks;
    outf << stepMean << " ";
    outf << sqrt((mean2Accumulator/nBlocks - pow(stepMean, 2))/(nBlocks-1)) << std::endl;
  }

  outf.close();
  //posOut.close(); //debug
  for (int i = 0; i < nBlocks; i++) {
    delete[] positions[i];
  }
  delete[] positions;
}
