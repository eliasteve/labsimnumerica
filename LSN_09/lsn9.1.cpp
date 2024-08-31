#include <iostream>
#include <fstream>
#include <string>
#include "individual.h"
#include "../random/random.h"
#include <cmath>
#include <cstdlib> // exit()
#include <iomanip> // setprecision()


//Run one iteration of random search
void randomSearchIteration(std::vector<individual>&, Random&, arma::mat&);
//Mutate one individual applying one of: pair permutation, shift,
//permutation of ranges, inversion, or no mutation at all according to the
//probabilities passed as argument
void mutate(individual&, Random&, double, double, double, double, arma::mat&);
//Run one iteration of the genetic algorithm
std::vector<individual> geneticAlgoIteration(std::vector<individual>&, Random&, arma::mat&);
//Sort the population and returns the fitness of the best individual
double sortPopulation(std::vector<individual>&);
//Run a genetic algorithm in order to find the solution to the travelling
//salesman problem
void optimizeGeneticAlgo(std::vector<individual>&, Random&, arma::mat&, int, int, std::string, std::string);
//Run a random search (mutations only, no crossover) in order to find the
//solution to the travelling salesman problem
void optimizeRandomSearch(std::vector<individual>&, Random&, arma::mat&, int, int, std::string, std::string);
//Write the matrix of cities to tile
void writeCities(arma::mat&, std::string);
//Write an individual path and its fitness to file
void writeIndividual(individual&, std::string);
//Compute mean and standard deviation of the best half of a population of
//paths
void computeStatisticsOfBestHalf(std::vector<individual>&, double&, double&);
//Compute the lookup table of the distance between pairs of cities
void computeDistanceLUT(int, arma::mat&, arma::mat&);

int main() {
  Random rng("../random/seed.in", "../random/primes32001.in", 2);

  //std::cout << std::setprecision(17); 
  int nCities = 34; //number of cities
  //Number of individual in the population of solutions
  int populationSize = 100;
  double theta = 0; //For generating cities in a circle
  arma::mat cities(nCities, 2); //(nCities, 2) matrix to hold the cities
  //Lookup table for distances between pairs of cities. Technically
  //something like memoization with a sparse matrix would have been more
  //efficient, but for a small number of cities we can eat the overhead
  //of computing, storing the whole lookup table (armadillo documentation
  //discourages the use of sparse matrices for size <= (100, 100) anyways)
  arma::mat distanceLUT(nCities, nCities);
  //Population of solutions
  std::vector<individual> population(populationSize);
  //copy of the population, to start both RS and GA from the same state
  std::vector<individual> populationCopy;

  // CITIES UNIFRORMY DISTRIBUTED IN A CIRCLE

  //Initialize city position
  for (int i = 0; i < nCities; i++) {
    theta = rng.Rannyu(0, 2*M_PI);
    cities(i, 0) = cos(theta); // x-coordinate of city
    cities(i, 1) = sin(theta); // y-coordinate of city
  }

  //Compute lookup table of distances between cities
  computeDistanceLUT(nCities, cities, distanceLUT);

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

  std::cout << "Beginning random search (circle)..." << std::endl;
  optimizeRandomSearch(population, rng, distanceLUT, 100000, 3000, "statistics_random_search_circle.dat", "prog_best_random_search_circle.dat");

  //Write best individual (if population is sorted it's the first one)
  writeIndividual(population[0], "path_random_search_circle.dat");
  writeCities(cities, "cities_circle.dat");

  //Testing crossover
  /*
  std::array<individual, 2> children;
  for (int i = 0; i < populationSize; i+=2) {
    crossover(population[i], population[i+1], children[0], children[1], rng, distanceLUT);
    std::cout << population[i].path.t() << population[i+1].path.t() << children[0].path.t() << children[1].path.t() << std::endl;
  }
  */

  //Optimization via genetic algorithm
  population = populationCopy; //Start from the same initial state as RS
  std::cout << "Beginning optimization via genetic algorithm (circle)..." << std::endl;
  optimizeGeneticAlgo(population, rng, distanceLUT, 100000, 3000, "statistics_genetic_algo_circle.dat", "prog_best_genetic_algo_circle.dat");
  //Write best individual (if population is sorted it's the first one)
  writeIndividual(population[0], "path_genetic_algo_circle.dat");


  // CITIES UNIFORMLY DISTRIBUTED ON A SQUARE
  
  //Initialize city position
  for (int i = 0; i < nCities; i++) {
    cities(i, 0) = rng.Rannyu(); // x-coordinate of city
    cities(i, 1) = rng.Rannyu(); // y-coordinate of city
  }

  //Recompute lookup table of distance bewteen cities
  computeDistanceLUT(nCities, cities, distanceLUT);

  //Use the same initial state as the circle case to avoid generating
  //another initial state, recompute the fitness using the new lookup table.
  for (int i = 0; i < populationSize; i++) {
    populationCopy[i].fitness = computeFitness(populationCopy[i], distanceLUT);
  }

  //Random search

  population = populationCopy;
  std::cout << "Beginning random search (square)..." << std::endl;
  optimizeRandomSearch(population, rng, distanceLUT, 100000, 3000, "statistics_random_search_square.dat", "prog_best_random_search_square.dat");

  //Write best individual (if population is sorted it's the first one)
  writeIndividual(population[0], "path_random_search_square.dat");
  writeCities(cities, "cities_square.dat");

  //Optimization via genetic algorithm
  population = populationCopy;
  std::cout << "Beginning optimization via genetic algorithm (square)..." << std::endl;
  optimizeGeneticAlgo(population, rng, distanceLUT, 100000, 3000, "statistics_genetic_algo_square.dat", "prog_best_genetic_algo_square.dat");
  //Write best individual (if population is sorted it's the first one)
  writeIndividual(population[0], "path_genetic_algo_square.dat");

  return 0;
}



//Run one iteration of random search
void randomSearchIteration(
    std::vector<individual>& population, //Population of paths
    Random& rng, //Random number generator
    arma::mat& distLUT //Lookup table of distances between pairs of cities
) {
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



//Mutate one individual applying one of: pair permutation, shift,
//permutation of ranges, inversion, or no mutation at all according to the
//probabilities passed as argument
void mutate(
  individual& ind, //Individual to be mutated
  Random& rng, //Random number generator
  double pPairPerm, //Probability of pair permutation
  double pShift, //Probability of shift
  double pPermRanges, //Probability of permutation of ranges
  double pInv, //Probability of inversion
  arma::mat& distLUT //Lookup table of distances between pairs of cities
) {
  double randomVal = rng.Rannyu(); //To select which mutation to apply
  //Compute this separately to check the probabilities sum to something
  //less than 1
  double final_prob = pPairPerm + pShift + pPermRanges + pInv;
  if(final_prob > 1) {
    std::cerr << "Sum of probabilities in mutate exceeds 1. I'm pretty sure you didn't want this." << std::endl;
  }

  if (randomVal < pPairPerm) pairPermutation(ind, distLUT);
  if (randomVal < pPairPerm + pShift) shift(ind, rng, distLUT);
  if (randomVal < pPairPerm + pShift + pPermRanges) permuteRanges(ind, rng, distLUT);
  if (randomVal < final_prob) inversion(ind, rng, distLUT);
}



//Run one iteration of the genetic algorithm
std::vector<individual> geneticAlgoIteration(
    std::vector<individual>& population, //Population of paths
    Random& rng, //Random number generator
    arma::mat& distLUT //Lookup table of distances between pairs of cities
) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to geneticAlgoIteration. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  double pCrossover = 0.90; //Probability of crossover
  double pPairPerm = 0.10; //Probability of pair permutation
  double pShift = 0.10; //Probability of shift
  double pPermRanges = 0.10; //Probability of permutation of ranges
  double pInv = 0.10; //Probability of inversion

  //New population of paths
  std::vector<individual> newPopulation(population.size());
  //To store indices for individual selected to crossover
  std::array<int, 2> selectedIndices;

  for (int i = 0; i < population.size(); i+=2) { 
    selectedIndices = select(population.size(), rng);
    
    //Either we apply crossover of we copy the parents into the new
    //population
    if (rng.Rannyu() < pCrossover) {
      crossover(population[selectedIndices[0]], population[selectedIndices[1]], newPopulation[i], newPopulation[i+1], rng, distLUT);
    } else {
      newPopulation[i] = population[selectedIndices[0]];
      newPopulation[i+1] = population[selectedIndices[1]];
    }
    mutate(newPopulation[i], rng, pPairPerm, pShift, pPermRanges, pInv, distLUT);
    mutate(newPopulation[i+1], rng, pPairPerm, pShift, pPermRanges, pInv, distLUT);
  }

  return newPopulation;
}



//Sort the population and returns the fitness of the best individual
double sortPopulation(
  std::vector<individual>& population //Population to be sorted
) {
  std::sort(population.begin(), population.end());
  return population[0].fitness;
}



//Run a genetic algorithm in order to find the solution to the travelling
//salesman problem
void optimizeGeneticAlgo(
    std::vector<individual>& population, //Initlal population of paths
    Random& rng, //Random number generator
    arma::mat& distanceLUT, //Lookup table for distance of pairs of cities
    int nSteps, //Max number of steps to perform
    //After this many steps without an improvement in fitness, stop the
    //algorithm
    int maxTimesUnchanged,
    //Filename for saving the mean fitness of the best half of paths. Empty
    //string disables saving
    std::string filenameMean, 
    //Filename for saving the best fitness per iteration. Empty string 
    //disables saving
    std::string filenameBest
) {
  //Will store the best fitness of the current iteration
  double bestFitness = 0; 
  double bestFitnessPrev = sortPopulation(population); //Best fitness of the
                                                       //previous iteration
  //Counter for iterations. Need to access it from outside of its for loop
  //to print the number of iteration once we're done with the algorithm
  int i;
  //Counts how many iterations have passed without an improvement in fitness
  int timesUnchanged = 0;
  //Mean standard deviation of the fitness of the best half of paths
  double mean = 0, std = 0;
  //File for saving the mean, standard deviation of the fitness of the best
  //half of paths
  std::ofstream outfMean;
  //File for saving the best fitness at each iteration
  std::ofstream outfBest;

  if (filenameMean != "") {
    outfMean.open(filenameMean);
    if (!outfMean.is_open()) {
      std::cerr << "Unable to open file " + filenameMean + ", statistics of the best half will NOT be written." << std::endl;
    }
    outfMean << "Iteration Mean StdDeviation" << std::endl;
  }

  if (filenameBest != "") {
    outfBest.open(filenameBest);
    if (!outfBest.is_open()) {
      std::cerr << "Unable to open file " + filenameBest + ", fitness of best individual per iteration will NOT be written." << std::endl;
    }
    outfBest << "Iteration Best" << std::endl;
  }

  for (i = 0; (i < nSteps) && (timesUnchanged < maxTimesUnchanged); i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);
    //Saving of best fitness and statistics of the best half
    if(filenameMean != "") {
      computeStatisticsOfBestHalf(population, mean, std);
      outfMean << i+1 << " " << mean << " " << std << std::endl;
    }
    if(filenameBest != "") {
      outfBest << i+1 << " " << bestFitness << std::endl;
    }

    //Check if there have been any improvements
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      std::cout << "Improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }

  if(filenameMean != "") outfMean.close();
  if(filenameBest != "") outfBest.close();

  std::cout << "Exited after " << i << " iterations" << std::endl;
}



//Run a random search (mutations only, no crossover) in order to find the
//solution to the travelling salesman problem
void optimizeRandomSearch(
    std::vector<individual>& population, //Initlal population of paths
    Random& rng, //Random number generator
    arma::mat& distanceLUT, //Lookup table for distance of pairs of cities
    int nSteps, //Max number of steps to perform
    //After this many steps without an improvement in fitness, stop the
    //algorithm
    int maxTimesUnchanged,
    //Filename for saving the mean fitness of the best half of paths. Empty
    //string disables saving
    std::string filenameMean, 
    //Filename for saving the best fitness per iteration. Empty string 
    //disables saving
    std::string filenameBest
) {
  //Will store the best fitness of the current iteration
  double bestFitness = 0; 
  double bestFitnessPrev = sortPopulation(population); //Best fitness of the
                                                       //previous iteration
  //Counter for iterations. Need to access it from outside of its for loop
  //to print the number of iteration once we're done with the algorithm
  int i;
  //Counts how many iterations have passed without an improvement in fitness
  int timesUnchanged = 0;
  //Mean standard deviation of the fitness of the best half of paths
  double mean = 0, std = 0;
  //File for saving the mean, standard deviation of the fitness of the best
  //half of paths
  std::ofstream outfMean;
  //File for saving the best fitness at each iteration
  std::ofstream outfBest;

  if (filenameMean != "") {
    outfMean.open(filenameMean);
    if (!outfMean.is_open()) {
      std::cerr << "Unable to open file " + filenameMean + ", statistics of the best half will NOT be written." << std::endl;
    }
    outfMean << "Iteration Mean StdDeviation" << std::endl;
  }
  if (filenameBest != "") {
    outfBest.open(filenameBest);
    if (!outfBest.is_open()) {
      std::cerr << "Unable to open file " + filenameBest + ", fitness of best individual per iteration will NOT be written." << std::endl;
    }
    outfBest << "Iteration Best" << std::endl;
  }

  for (i = 0; (i < 100000) && (timesUnchanged < 3000); i++) {
    randomSearchIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);
    
    //Saving of best fitness and statistics of the best half
    if(filenameMean != "") {
      computeStatisticsOfBestHalf(population, mean, std);
      outfMean << i+1 << " " << mean << " " << std << std::endl;
    }
    if(filenameBest != "") {
      outfBest << i+1 << " " << bestFitness << std::endl;
    }
    
    //Check if there have been any improvements
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      std::cout << "Improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }

  if(filenameMean != "") outfMean.close();
  if(filenameBest != "") outfBest.close();

  std::cout << "Exited after " << i << " iterations" << std::endl;
}



//Write the matrix of cities to file
void writeCities(
  arma::mat& cities, //Cities to write
  std::string filename //Name of output file
) {
  std::ofstream outf(filename);
  if(!outf.is_open()) {
    std::cout << "IO error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }
  outf << cities;
  outf.close();
}



//Write an individual path and its fitness to file
void writeIndividual(
  individual& ind, //Individual to write
  std::string filename //Name of output file
) {
  std::ofstream outf(filename);
  if(!outf.is_open()) {
    std::cout << "IO error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }
  outf << "Fitness: " << ind.fitness << std::endl << ind.path;
  outf.close();
}



//Compute mean and standard deviation of the best half of a population of
//paths
void computeStatisticsOfBestHalf(
    std::vector<individual> &population, //Population of file
    double& mean, //Output: mean of best half
    double& std //Output: standard deviation of best half
) {
  int N = ceil(population.size()/2.);
  //std::cout << N << std::endl;
  double meanAccum = 0;
  double mean2Accum = 0;

  for (int i = 0; i < N; i++) {
    meanAccum += population[i].fitness;
    mean2Accum += pow(population[i].fitness, 2);
  }

  mean = meanAccum/N;
  //Take the absolute value in case the variance is zero and a small
  //negative number is the numeric result of the computation.
  std = sqrt(abs(mean2Accum/N - pow(mean, 2)));
}



//Compute the lookup table of the distance between pairs of cities
void computeDistanceLUT(
  int nCities, //Number of cities
  arma::mat& cities, //(nCities, 2) matrix of cities
  arma::mat& distanceLUT //Output: Lookup table
) {
  std::cout << "Computing lookup table of distances between cities... ";
  distanceLUT.zeros();
  for (int i = 0; i < nCities; i++) {
    for (int j = i+1; j < nCities; j++) {
      distanceLUT(i, j) = vecnorm(cities.row(i) - cities.row(j));
      //std::cout << cities.row(i) << std::endl << cities.row(j) << std::endl;
      //std::cout << vecnorm(cities.row(i) - cities.row(j)) << std::endl;
    }
  }
  //Gives me the full matrix because at this point the bottom triangle 
  //contains only zeroes
  distanceLUT = distanceLUT + trans(distanceLUT);
  std::cout << "done." << std::endl;
}
