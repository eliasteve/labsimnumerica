#include "individual.h"
#include "../random/random.h"
#include <cmath>
#include <cstdlib> // exit()
#include <iomanip> // setprecision()


void randomSearchIteration(std::vector<individual>&, Random&, arma::mat&);
void mutate(individual&, Random&, double, double, double, double, arma::mat&);
std::vector<individual> geneticAlgoIteration(std::vector<individual>&, Random&, arma::mat&);
double sortPopulation(std::vector<individual>&);
void optimizeGeneticAlgo(std::vector<individual>&, Random&, arma::mat&, int, int);
void optimizeRandomSearch(std::vector<individual>&, Random&, arma::mat&, int, int);
void writeCities(arma::mat&);
void writeIndividual(individual&, std::string);

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 2);

  //std::cout << std::setprecision(17); 
  int nCities = 34, populationSize = 100;
  double circleRadius = 1, theta = 0;
  arma::mat cities(nCities, 2), distanceLUT(nCities, nCities);
  std::vector<individual> population(populationSize), populationCopy;

  //Initialize city position
  for (int i = 0; i < nCities; i++) {
    theta = rng.Rannyu(0, 2*M_PI);
    cities(i, 0) = circleRadius*cos(theta); // x-coordinate of city
    cities(i, 1) = circleRadius*sin(theta); // y-coordinate of city
  }

  /*
  for (int i = 0; i < nCities; i++) {
    cities(i, 0) = rng.Rannyu(); // x-coordinate of city
    cities(i, 1) = rng.Rannyu(); // y-coordinate of city
  }
  */

  //Compute lookup table of distances between cities
  for (int i = 0; i < nCities; i++) {
    for (int j = i+1; j < nCities; j++) {
      distanceLUT(i, j) = vecnorm(cities.row(i) - cities.row(j));
    }
  }
  distanceLUT = distanceLUT + trans(distanceLUT);

  /*
  //Testing of the check function
  individual test{arma::conv_to<arma::uvec>::from(arma::linspace(1, nCities, nCities)), 0};
  for (int i=0; i<100; i++) {
    test.path = arma::randperm(100, nCities);
    std::cout << checkConstraint(test) << std::endl; //Should fail
    test.path = arma::randperm(nCities);
    std::cout << checkConstraint(test) << std::endl; //Should succeed
  } */

  //Initialize population at random
  for (int i = 0; i < populationSize; i++) {
    population[i] = {join_cols(arma::uvec(1).zeros(), arma::randperm(nCities-1) + 1), 0}; //Ensure the 1st city is 0, randomly permute the others.
    checkConstraintAndTerminate(population[i]);
    population[i].fitness = computeFitness(population[i], distanceLUT);
    //std::cout << population[i].path << std::endl << population[i].fitness << std::endl;
  }
  //Copy the initial state in order to compare random search to GAs
  populationCopy = population;

  //Test of sorting by fitness
  /*std::sort(population.begin(), population.end());
  for (int i = 0; i < populationSize; i++) {
    std::cout << population[i].fitness << std::endl;
  } */

  //Test of pair permutation
  /*
  int index = 0;
  for (int i = 0; i < 100; i++) {
    index = rng.Rannyu(0, populationSize);
    //std::cout << population[index].path << std::endl;
    pairPermutation(population[index]);
    //std::cout << population[index].path << std::endl;
    //std::cout << "\n\n\n";
  } */

  //Test of shift
  /*
  for (int i = 0; i < populationSize; i++) {
    //std::cout << population[i].path << std::endl;
    shift(population[i], rng, distanceLUT);
    //std::cout << population[i].path << std::endl;
    //std::cout << "\n\n\n";
  } */

  //Test of permutation of ranges
  /*
  for (int i = 0; i < populationSize; i++) {
    //std::cout << population[i].path.t() << std::endl;
    permuteRanges(population[i], rng, distanceLUT);
    //std::cout << population[i].path.t() << std::endl;
    //std::cout << "\n\n\n";
  }*/

  //Test of inversion
  /*
  for (int i = 0; i < populationSize; i++) {
    std::cout << population[i].path.t() << std::endl;
    inversion(population[i], rng, distanceLUT);
    std::cout << population[i].path.t() << std::endl;
    std::cout << "\n\n\n";
  } */

  //Random search using only mutations

  sortPopulation(population);
  writeIndividual(population[0], "initial_best_path.dat");
  optimizeRandomSearch(population, rng, distanceLUT, 100000, 3000);
  /*int timesUnchanged = 0, i;
  double bestFitness = 0, bestFitnessPrev = 0;
  std::sort(population.begin(), population.end());
  writeIndividual(population[0], "initial_best_path.dat");
  bestFitnessPrev = population[0].fitness;

  for (i = 0; (i < 100000) && (timesUnchanged < 3000); i++) {
    bestFitness = randomSearchIteration(population, rng, distanceLUT);
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      std::cout << "Improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }
  std::cout << "Exited after " << i+1 << " iterations" << std::endl;*/


  writeIndividual(population[0], "path_random_search.dat");
  writeCities(cities);

  //Testing crossover
  /*
  std::array<individual, 2> children;
  for (int i = 0; i < populationSize; i+=2) {
    crossover(population[i], population[i+1], children[0], children[1], rng, distanceLUT);
    std::cout << population[i].path.t() << population[i+1].path.t() << children[0].path.t() << children[1].path.t() << std::endl;
  }
  */

  population = populationCopy;
  optimizeGeneticAlgo(population, rng, distanceLUT, 100000, 3000);
  writeIndividual(population[0], "path_genetic_algo.dat");
  return 0;
}



void randomSearchIteration(std::vector<individual>& population, Random& rng, arma::mat& distLUT) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to randomSearchIteration. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  //Keep the fittest individual unchanged, suppress the least fit
  population[population.size()-1] = population[0];
  //Mutate all remaining individuals
  for (int i = 0; i < population.size()-1; i++) { 
    mutate(population[i], rng, 0.25, 0.25, 0.25, 0.25, distLUT);
  }
}



void mutate(
  individual& ind,
  Random& rng, 
  double pPairPerm, 
  double pShift, 
  double pPermRanges,
  double pInv,
  arma::mat& distLUT
) {
  double randomVal = rng.Rannyu();
  double final_prob = pPairPerm + pShift + pPermRanges + pInv;
  if(final_prob > 1) {
    std::cerr << "Sum of probabilities in mutate exceeds 1. I'm pretty sure you didn't want this." << std::endl;
  }

  if (randomVal < pPairPerm) pairPermutation(ind, distLUT);
  if (randomVal < pPairPerm + pShift) shift(ind, rng, distLUT);
  if (randomVal < pPairPerm + pShift + pPermRanges) permuteRanges(ind, rng, distLUT);
  if (randomVal < final_prob) inversion(ind, rng, distLUT);
}



std::vector<individual> geneticAlgoIteration(std::vector<individual>& population, Random& rng, arma::mat& distLUT) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to geneticAlgoIteration. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  double pCrossover = 0.90, pPairPerm = 0.025, pShift = 0.025, pPermRanges = 0.025, pInv = 0.025;

  std::vector<individual> newPopulation(population.size());
  std::array<int, 2> selectedIndices;

  for (int i = 0; i < population.size(); i+=2) { 
    //std::cout << "Generating places " << i << ", " << i+1 << std::endl;
    selectedIndices = select(population.size(), rng);
    //std::cout << "Selected " << selectedIndices[0] << " " << selectedIndices[1] << std::endl;
    if (rng.Rannyu() < pCrossover) {
      //std::cout << "a";
      crossover(population[selectedIndices[0]], population[selectedIndices[1]], newPopulation[i], newPopulation[i+1], rng, distLUT);
    } else {
      newPopulation[i] = population[selectedIndices[0]];
      newPopulation[i+1] = population[selectedIndices[1]];
    }
    //std::cout << "b";
    mutate(newPopulation[i], rng, pPairPerm, pShift, pPermRanges, pInv, distLUT);
    //std::cout << "c" << std::endl;
    mutate(newPopulation[i+1], rng, pPairPerm, pShift, pPermRanges, pInv, distLUT);
  }

  return newPopulation;
}



double sortPopulation(std::vector<individual>& population) {
  //Sorts the population and returns the fitness of the best individual
  std::sort(population.begin(), population.end());
  return population[0].fitness;
}



void optimizeGeneticAlgo(std::vector<individual>& population, Random& rng, arma::mat& distanceLUT, int nSteps, int maxTimesUnchanged) {
  double bestFitness = 0, bestFitnessPrev = sortPopulation(population);
  int i, timesUnchanged = 0;
  for (i = 0; (i < nSteps) && (timesUnchanged < maxTimesUnchanged); i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);
    //std::cout << "best " << bestFitness << ", prev " << bestFitnessPrev << std::endl;
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      std::cout << "Improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }
  std::cout << "Exited after " << i+1 << " iterations" << std::endl;
}



void optimizeRandomSearch(std::vector<individual>& population, Random& rng, arma::mat& distanceLUT, int nSteps, int maxTimesUnchanged) {
  int timesUnchanged = 0, i;
  double bestFitness = 0, bestFitnessPrev = sortPopulation(population);

  for (i = 0; (i < 100000) && (timesUnchanged < 3000); i++) {
    randomSearchIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      std::cout << "Improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }
  std::cout << "Exited after " << i+1 << " iterations" << std::endl;
}



void writeCities(arma::mat& cities) {
  std::ofstream outf("cities.dat");
  if(!outf.is_open()) {
    std::cout << "IO error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }
  outf << cities;
  outf.close();
}



void writeIndividual(individual& ind, std::string filename) {
  std::ofstream outf(filename);
  if(!outf.is_open()) {
    std::cout << "IO error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }
  outf << "Fitness: " << ind.fitness << std::endl << ind.path;
  outf.close();
}
