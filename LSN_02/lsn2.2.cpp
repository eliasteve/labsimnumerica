#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"
#include "../lib_NSL/lib_NSL.h"

// Function to initialize an array of points with zeros
void fillWithZeros(point *array, int len);

// Function to calculate the squared norm of a point
double normSquared(point p);

// Function to perform a discrete step in the random walk
point discreteStep(Random& rng, double latticeConst);

// Function to perform a continuous step in the random walk
point continuousStep(Random& rng, double stepLength);

// Function to print a point to an output file
void printPoint(point p, std::ofstream& out);

// Function to print a point to the console
void printPoint(point p);

// Function to add two points (vector addition)
point add(point p1, point p2);

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
  double a = 1;                // Lattice constant for the discrete random walk (and
                               // length of step for the continuous one)

  // Perform random walk with discrete steps
  performRW(nBlocks, walkersPerBlock, nSteps, a, "discreteRW.dat", discreteStep, rng);

  // Perform random walk with continuous steps
  performRW(nBlocks, walkersPerBlock, nSteps, a, "continuousRW.dat", continuousStep, rng);

  return 0;
}

// Function to initialize an array of points with zeros
void fillWithZeros(point *array, int len) {
  for (int i = 0; i < len; i++) {
    array[i].x = 0;
    array[i].y = 0;
    array[i].z = 0;
  }
}

// Function to calculate the squared norm of a point
double normSquared(point p) {
  return pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
}

// Function to perform a discrete step in the random walk
point discreteStep(Random& rng, double latticeConst) {
  double x = rng.Rannyu(), y = rng.Rannyu(0, 3);
  double toAdd = 0;
  point p = {0, 0, 0};
  
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
  double phi = rng.Rannyu(0, 2*M_PI), x = rng.Rannyu();
  double theta = acos(1-2*x);
  point step = {stepLength*sin(theta)*cos(phi), stepLength*sin(theta)*sin(phi), stepLength*cos(theta)};
  return step;
}

// Function to print a point to an output file
void printPoint(point p, std::ofstream& out) {
  out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

// Function to print a point to the console
void printPoint(point p) {
  std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

// Function to add two points (vector addition)
point add(point p1, point p2) {
  point result;
  result.x = p1.x + p2.x;
  result.y = p1.y + p2.y;
  result.z = p1.z + p2.z;
  return result;
}

// Function to perform a random walk simulation
void performRW(
  int nBlocks,
  int walkersPerBlock,
  int nSteps,
  int stepLength,
  std::string outFileName,
  std::function<point(Random&, double)> stepFunction,
  Random& rng
) {
  std::ofstream outf(outFileName);
  std::ofstream posOut("positions.dat");
  if(!outf.is_open()) {
    std::cerr << "I/O error opening " << outFileName << ", program continues." << std::endl;
    return;
  }

  point positions[nBlocks][walkersPerBlock];
  for (int i = 0; i < nBlocks; i++) {
    for (int j = 0; j < walkersPerBlock; j++) {
      positions[i][j] = {0, 0, 0};
    }
  }

  double blockNorm2Accumulator = 0;        // Accumulator for the norm squared in a block
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
        positions[k][j] = add(positions[k][j], stepFunction(rng, stepLength));
        printPoint(positions[k][j], posOut);
        posOut << " ";
        blockNorm2Accumulator += normSquared(positions[k][j]);
      }
      blockMean = sqrt(blockNorm2Accumulator/walkersPerBlock);
      meanAccumulator += blockMean;
      mean2Accumulator += pow(blockMean, 2);
    }
    posOut << std::endl;
    stepMean = meanAccumulator/nBlocks;
    outf << stepMean << " ";
    outf << sqrt((mean2Accumulator/nBlocks - pow(stepMean, 2))/(nBlocks-1)) << std::endl;
  }

  outf.close();
  posOut.close();
}
