#include "../lib_NSL/individual.h"
#include "../lib_NSL/optimization.h"
#include "../random/random.h"
#include <cmath>
#include <cstdlib> // exit()
#include <iomanip> // setprecision()
#include "mpi.h"


void optimizeGeneticAlgo(std::vector<individual>&, Random&, arma::mat&, int, int, int, int, int);
void optimizeGeneticAlgoNoMigr(std::vector<individual>&, Random&, arma::mat&, int, int, int, int);
arma::uvec generateExchangeList(int, int);
void migrate(std::vector<individual>&, int, int);
void sendIndividual(individual&, int, MPI_Request&, int, int);
void receiveIndividual(individual&, int, MPI_Status&, int, int);
individual bestOverallIndividual(std::vector<individual>&, int, int);
//void createIndividualMPIDatatype(MPI_Datatype&, int);

int main(int argc, char* argv[]) {
  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Need to use different seeds on each node, otherwise we just run the same
  //calculation on all nodes.
  Random rng("../random/seed.in", "../random/primes32001.in", 1+2*size+rank);
  arma::arma_rng::set_seed(2*size+rank); //For population generation and some mutation

  int nCities = 110, populationSize = 100;
  arma::mat cities(nCities, 2), distanceLUT(nCities, nCities);
  std::vector<individual> population(populationSize);
  std::string title = "best_path_GA_node" + std::to_string(rank) + ".dat";
  individual bestIndividual;

  //Initialize city positions
  if (rank == 0) {
    std::ifstream inf("cap_prov_ita.dat");
    if (!inf.is_open()) {
      std::cerr << "I/O ERROR: PROGRAM TERMINATES." << std::endl;
      exit(EXIT_FAILURE);
    }
    //Number of rows in the file is guaranteed to be 110, i.e. equal to nCities.
    for (int i = 0; i < nCities; i++) {
      inf >> cities(i, 0) >> cities(i, 1);
    }
    //writeCities(cities, "cities.dat");
  }
  MPI_Bcast(&cities(0, 0), 2*nCities, MPI_REAL8, 0, MPI_COMM_WORLD);
  //writeCities(cities, "cities_node" + std::to_string(rank));


  //Compute lookup table of distances between cities
  if (rank == 0) computeDistanceLUT(nCities, cities, distanceLUT);
  MPI_Bcast(&distanceLUT(0, 0), nCities*nCities, MPI_REAL8, 0, MPI_COMM_WORLD);

  //Initialize population at random
  for (int i = 0; i < populationSize; i++) {
    population[i] = {join_cols(arma::uvec(1).zeros(), arma::randperm(nCities-1) + 1), 0}; //Ensure the 1st city is 0, randomly permute the others.
    checkConstraintAndTerminate(population[i]);
    population[i].fitness = computeFitness(population[i], distanceLUT);
  }

  //Copy initial population so we can compare the results of the search with and without migration
  std::vector<individual> populationCopy(population);


  //Search with migration
  optimizeGeneticAlgo(populationCopy, rng, distanceLUT, 200000, 10000, 50, size, rank);
  bestIndividual = bestOverallIndividual(populationCopy, size, rank);
  if (rank==0) writeIndividual(bestIndividual, "best_individual.dat");


  //Search without migration
  optimizeGeneticAlgoNoMigr(population, rng, distanceLUT, 200000, 10000, size, rank);
  bestIndividual = bestOverallIndividual(population, size, rank);
  if (rank==0) writeIndividual(bestIndividual, "best_individual_nomigr.dat");
  
  MPI_Finalize();
  return 0;
}



void optimizeGeneticAlgo(std::vector<individual>& population, Random& rng, arma::mat& distanceLUT, int nSteps, int maxTimesUnchanged, int migrationPeriod, int size, int rank) {
  double bestFitness = 0, bestFitnessPrev = sortPopulation(population);
  bool wantToExit = false, exit = false;
  int i, timesUnchanged = 0;
  std::ofstream outf("node" + std::to_string(rank) + "_output.txt");
  //std::ofstream debugExit("exitlog_node" + std::to_string(rank));
  if (rank == 0) writeIndividual(population[0], "best_initial_individual.dat");
  for (i = 0; (i < nSteps) && !exit; i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);
    if(i % migrationPeriod == migrationPeriod-1) {
      outf << "Iteration " << i << ": exchanging individuals. Best fitness before the exchange: ";
      outf << sortPopulation(population) << std::endl;
      migrate(population, size, rank);
      outf << "Best individual received has fitness " << population[population.size()-1].fitness << std::endl;
    }
    bestFitness = sortPopulation(population);
    //outf << "best " << bestFitness << ", prev " << bestFitnessPrev << std::endl;
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      outf << "Iteration " << i << ": improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
    wantToExit = (timesUnchanged > maxTimesUnchanged);
    MPI_Allreduce(&wantToExit, &exit, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    //debugExit << "Iteration " << i << ": wantToExit is " << wantToExit << ", exit is " << exit << std::endl;
  }
  outf << "Exited after " << i << " iterations" << std::endl;
  outf.close();
  //debugExit.close();
}



void optimizeGeneticAlgoNoMigr(std::vector<individual>& population, Random& rng, arma::mat& distanceLUT, int nSteps, int maxTimesUnchanged, int size, int rank) {
  double bestFitness = 0, bestFitnessPrev = sortPopulation(population);
  int i, timesUnchanged = 0;
  std::ofstream outf("node" + std::to_string(rank) + "_output_NM.txt");
  if (rank == 0) writeIndividual(population[0], "best_initial_individual_NM.dat");
  for (i = 0; (i < nSteps) && (timesUnchanged < maxTimesUnchanged); i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);
    //outf << "best " << bestFitness << ", prev " << bestFitnessPrev << std::endl;
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      outf << "Iteration " << i << ": improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
  }
  outf << "Exited after " << i << " iterations" << std::endl;
  outf.close();
}



arma::uvec generateExchangeList(int size, int rank) {
  arma::uvec exchangeList(size, arma::fill::zeros);
  //Generate the list of exchanges as a permitation of numbers from 0 to size-1. Each node sends to the one after it in the array, and receives from the one before it.
  if (rank == 0) exchangeList = arma::randperm(size);
  MPI_Bcast(&exchangeList[0], size, MPI_INTEGER8, 0, MPI_COMM_WORLD);
  return exchangeList;
}



void migrate(std::vector<individual>& population, int size, int rank) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to migrate. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  int sendTo, receiveFrom, nIndividuals = 5;
  MPI_Request req;
  MPI_Status stat;
  arma::uvec exchangeList = generateExchangeList(size, rank);
  for (int i = 0; i < size; i++) {
    if(rank == exchangeList[i]) {
      sendTo = exchangeList[(i+1)%size];
      receiveFrom = exchangeList[(size+i-1)%size]; //Add size to prevent index before modulus from being negative.
    }
  }
  
  /*
  std::ofstream outf("out_" + std::to_string(rank) + ".txt");
  outf << "Exchange list: " << exchangeList.t() << std::endl;
  outf << "Node " << rank << " sends to " << sendTo << ", receives from " << receiveFrom << std::endl;
  */
  /*
  for (int i = 0; i < population.size(); i++) {
    outf << "Individual number " << i+1 << ":" << std::endl;
    outf << population[i].path << std::endl;
    outf << "Fitness " << population[i].fitness << std::endl;
  }
  */

  for(int i = 0; i < nIndividuals; i++) {
    //Comm ID is rank of sender for path, size+rank of sender for fitness
    /*
    outf << "Sending individual with path" << std::endl;
    outf << population[i].path << std::endl;
    outf << "and fitness " << population[i].fitness << std::endl;
    */
    sendIndividual(population[i], sendTo, req, size, rank);
    receiveIndividual(population[population.size()-1-i], receiveFrom, stat, size, rank);
    /*
    outf << "Received individual with path" << std::endl;
    outf << population[population.size()-1-i].path << std::endl;
    outf << "and fitness " << population[population.size()-1-i].fitness << std::endl;
    */
  }

  //outf.close();
}



void sendIndividual(individual& ind, int toNode, MPI_Request& req, int size, int rank) {
  MPI_Isend(&ind.path[0], ind.path.n_elem, MPI_INTEGER8, toNode, rank, MPI_COMM_WORLD, &req);
  MPI_Isend(&ind.fitness, 1, MPI_REAL8, toNode, size+rank, MPI_COMM_WORLD, &req);
}



void receiveIndividual(individual& ind, int receiveFrom, MPI_Status& status, int size, int rank) {
  MPI_Recv(&ind.path[0], ind.path.n_elem, MPI_INTEGER8, receiveFrom, receiveFrom, MPI_COMM_WORLD, &status);
  MPI_Recv(&ind.fitness, 1, MPI_REAL8, receiveFrom, size+receiveFrom, MPI_COMM_WORLD, &status);
}



individual bestOverallIndividual(std::vector<individual>& population, int size, int rank) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to bestOverallIndividual. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  MPI_Request req;
  MPI_Status status;
  if (rank != 0) {
    sendIndividual(population[0], 0, req, size, rank);
  }

  if (rank == 0) {
    std::vector<individual> candidates(size);
    individual tmp = population[0]; //Initialize this to some value so we receive on memory that's allocated
    candidates[0] = population[0];
    for (int i = 1; i < size; i++) {
      receiveIndividual(tmp, i, status, size, rank);
      candidates[i] = tmp;
    }
    std::sort(candidates.begin(), candidates.end());
    return candidates[0];
  }

  return population[0];
}
