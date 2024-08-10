#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"
#include "../lib_NSL/lib_NSL.h"

/*struct point{
  double x;
  double y;
  double z;
};*/

void fillWithZeros(point *array, int len) {
  for (int i = 0; i < len; i++) {
    array[i].x = 0;
    array[i].y = 0;
    array[i].z = 0;
  }
}

double normSquared(point p) {
  return pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
}

point discreteStep(Random& rng, double latticeConst) {
  // y: which direction, x: up or down
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

point continuousStep(Random& rng, double stepLength) {
  double phi = rng.Rannyu(0, 2*M_PI), x = rng.Rannyu();
  double theta = acos(1-2*x);
  point step = {stepLength*sin(theta)*cos(phi), stepLength*sin(theta)*sin(phi), stepLength*cos(theta)};
  return step;
}

void printPoint(point p, std::ofstream& out) {
  out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

void printPoint(point p) {
  std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

point add(point p1, point p2) {
  //printPoint(p1); std::cout << " ";
  //printPoint(p2); std::cout << " ";
  point result;
  //std::cout << p1.x << " " << p1.y << " " << p1.z << " " << p2.x << " " << p2.y << " " << p2.z << std::endl;
  result.x = p1.x + p2.x;
  result.y = p1.y + p2.y;
  result.z = p1.z + p2.z;
  //printPoint(result); std::cout << std::endl;
  return result;
}

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

  double blockNorm2Accumulator = 0, meanAccumulator = 0, mean2Accumulator = 0, blockMean = 0, stepMean = 0;

  for (int i = 0; i < nSteps; i++) {
    meanAccumulator = 0;
    mean2Accumulator = 0;
    for (int k = 0; k < nBlocks; k++) {
      blockNorm2Accumulator = 0;
      for (int j = 0; j < walkersPerBlock; j++) {
        //std::cout << positions[k][j] << " -> ";
        positions[k][j] = add(positions[k][j], stepFunction(rng, stepLength));
        //std::cout << positions[k][j] << std::endl;
        printPoint(positions[k][j], posOut);
        posOut << " ";
        blockNorm2Accumulator += normSquared(positions[k][j]);
      }
      blockMean = sqrt(blockNorm2Accumulator/walkersPerBlock);
      //std::cout << blockMean << std::endl;
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

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 1);

  //std::ofstream discrete("lsn2.2_discreteRW.dat");
  //std::ofstream posOut("positions.dat");

  /*if (!discrete.is_open()) {
    std::cerr << "Error trying to open output file lsn2.2_discreteRW.dat. The program terminates." << std::endl;
    return -1;
  }*/

  int nBlocks = 100, walkersPerBlock = 100, nSteps = 100;
  double a = 1; // Lattice constant for the discrete case
  //double blockNorm2Accumulator = 0, meanAccumulator = 0, mean2Accumulator = 0, blockMean = 0, mean = 0;

  performRW(nBlocks, walkersPerBlock, nSteps, a, "discreteRW.dat", discreteStep, rng);
  performRW(nBlocks, walkersPerBlock, nSteps, a, "continuousRW.dat", continuousStep, rng);
  /*point positions[nBlocks*walkersPerBlock];
  fillWithZeros(positions, nBlocks*walkersPerBlock);

  for (int i = 0; i < nSteps; i++) {
    meanAccumulator = 0;
    mean2Accumulator = 0;
    for (int k = 0; k < nBlocks; k++) {
      blockNorm2Accumulator = 0;
      for (int j = 0; j < walkersPerBlock; j++) {
        positionIndex = walkersPerBlock*k + j;
        positions[positionIndex] = add(positions[positionIndex], discreteStep(&rng, a));
        printPoint(positions[positionIndex], &posOut);
        posOut << " ";
        blockNorm2Accumulator += normSquared(positions[positionIndex]);
      }
      blockMean = sqrt(blockNorm2Accumulator/walkersPerBlock);
      meanAccumulator += blockMean;
      mean2Accumulator += pow(blockMean, 2);
    }
    posOut << std::endl;
    mean = meanAccumulator/nBlocks;
    discrete << mean << " ";
    discrete << sqrt((mean2Accumulator/nBlocks - pow(mean, 2))/(nBlocks-1)) << std::endl;
  }

  discrete.close();
  posOut.close();
  */

  /*
  point test = {3, 1, -2}, res, step;
  printPoint(test); std::cout << std::endl;
  for (int i = 0; i < 100; i++) {
    step = discreteStep(&rng, a);
    res = add(step, test);
    printPoint(step);
    std::cout << " ";
    printPoint(res);
    std::cout << std::endl;
  } */

  return 0;
}
