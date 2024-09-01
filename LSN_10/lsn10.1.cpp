#include "../lib_NSL/individual.h"
#include "../lib_NSL/optimization.h"
#include "../random/random.h"
#include <cmath>
#include <cstdlib> // exit()
#include <iomanip> // setprecision()
#include "mpi.h"


//Produce a solution to the travelling salesman problem using genetic search
//(parallelized genetic algorithm + migration)
void optimizeGeneticAlgo(std::vector<individual>&, Random&, arma::mat&, int, int, int, int, int);
//Run a genetic algorithm in order to find the solution to the travelling
//salesman problem
void optimizeGeneticAlgoNoMigr(std::vector<individual>&, Random&, arma::mat&, int, int, int, int);
//Generate an exchange list for migration of individuals in genetic search
arma::uvec generateExchangeList(int, int);
//Perform migration while running a genetic search
void migrate(std::vector<individual>&, int, int);
//Sends an individual to another node
void sendIndividual(individual&, int, MPI_Request&, int, int);
//Receives an individual from another node
void receiveIndividual(individual&, int, MPI_Status&, int, int);
//Find the best individual among the fittest in each node
individual bestOverallIndividual(std::vector<individual>&, int, int);

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



//Produce a solution to the travelling salesman problem using genetic search
//(parallelized genetic algorithm + migration)
void optimizeGeneticAlgo(
  std::vector<individual>& population, //Initial population of paths
  Random& rng, //Random number generator
  arma::mat& distanceLUT, //Lookup table for distance of pairs of cities
  int nSteps, //Max number of steps to perform
  //After this many steps without improvement, the node sets wantToExit to
  //true (see below)
  int maxTimesUnchanged,
  //ater this many steps, execute a migration
  int migrationPeriod,
  //MPI parameters
  int size, //Number of nodes
  int rank //Label of node
) {
  //Will store the best fitness of the current iteration
  double bestFitness = 0; 
  double bestFitnessPrev = sortPopulation(population); //Best fitness of the
                                                       //previous iteration

  //With migrations, all nodes need to exit simultaneously (otherwise a
  //node might send an individual to or receive it from a node which has 
  //exited already, and the program hangs). The exit policy is this: after
  //maxTimesUnchanged iteration without improvement, a node sets wantToExit
  //to true. When all nodes have done so, exit gets set to true, and all 
  //the nodes exit simultaneously.
  bool wantToExit = false, exit = false;
  //Counter for iterations. Need to access it from outside of its for loop
  //to print the number of iteration once we're done with the algorithm
  int i;
  //Counts how many iterations have passed without an improvement in fitness
  int timesUnchanged = 0;

  //Output file for the node (for legibility, so the nodes don't all try to
  //write to stdout).
  std::ofstream outf("node" + std::to_string(rank) + "_output.txt");
  //Output file for debigging the cose for synchronized exit
  //std::ofstream debugExit("exitlog_node" + std::to_string(rank));

  if (rank == 0) writeIndividual(population[0], "best_initial_individual.dat");
  for (i = 0; (i < nSteps) && !exit; i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);

    //Migration
    if(i % migrationPeriod == migrationPeriod-1) {
      outf << "Iteration " << i << ": exchanging individuals. Best fitness before the exchange: ";
      outf << sortPopulation(population) << std::endl;
      migrate(population, size, rank);
      outf << "Best individual received has fitness " << population[population.size()-1].fitness << std::endl;
    }

    bestFitness = sortPopulation(population);
    //outf << "best " << bestFitness << ", prev " << bestFitnessPrev << std::endl;

    //Synchronized exit code
    if(!(bestFitness < bestFitnessPrev)) {
      timesUnchanged++;
    }
    else {
      timesUnchanged = 0;
      outf << "Iteration " << i << ": improved fitness: " << bestFitnessPrev << "->" << bestFitness << std::endl;
      bestFitnessPrev = bestFitness;
    }
    wantToExit = (timesUnchanged > maxTimesUnchanged);
    //After this call, exit is true on all nodes if wantToExit is true on
    //all nodes, false otherwise
    MPI_Allreduce(&wantToExit, &exit, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    //debugExit << "Iteration " << i << ": wantToExit is " << wantToExit << ", exit is " << exit << std::endl;
  }
  outf << "Exited after " << i << " iterations" << std::endl;
  outf.close();
  //debugExit.close();
}



//Run a genetic algorithm in order to find the solution to the travelling
//salesman problem
void optimizeGeneticAlgoNoMigr(
  std::vector<individual>& population, //Initlal population of paths
  Random& rng, //Random number generator
  arma::mat& distanceLUT, //Lookup table for distance of pairs of cities
  int nSteps, //Max number of steps to perform
  //After this many steps without an improvement in fitness, stop the
  //algorithm
  int maxTimesUnchanged,
  //Filename for saving the mean fitness of the best half of paths. Empty
  //string disables saving
  int size, //Number of nodes
  int rank //Label of node
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
  //Output file for the node (for legibility, so the nodes don't all try to
  //write to stdout).
  std::ofstream outf("node" + std::to_string(rank) + "_output_NM.txt");

  if (rank == 0) writeIndividual(population[0], "best_initial_individual_NM.dat");
  for (i = 0; (i < nSteps) && (timesUnchanged < maxTimesUnchanged); i++) {
    population = geneticAlgoIteration(population, rng, distanceLUT);
    bestFitness = sortPopulation(population);

    //Check if there have been improvements
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



//Generate an exchange list for migration of individuals in genetic search
arma::uvec generateExchangeList(
  int size, //Number of nodes
  int rank //Label of node
) {
  arma::uvec exchangeList(size, arma::fill::zeros); //Exchange list
  //Generate the list of exchanges as a permitation of numbers from 0 to size-1. Each node sends to the one after it in the array, and receives from the one before it.
  if (rank == 0) exchangeList = arma::randperm(size);
  MPI_Bcast(&exchangeList[0], size, MPI_INTEGER8, 0, MPI_COMM_WORLD);
  return exchangeList;
}



//Perform migration while running a genetic search
void migrate(
  std::vector<individual>& population, //Population of individuals in the
                                       //node
  int size, //Number of nodes
  int rank //Label of node
) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to migrate. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  int sendTo; //Rank of the node to send to
  int receiveFrom; //Rank of the node to receive from
  int nIndividuals = 5; //Number of individuals to be exchanged
  //Variables needed for communication with MPI
  MPI_Request req;
  MPI_Status stat; 
  //Exchange list, determines who each node sends to and receives from
  arma::uvec exchangeList = generateExchangeList(size, rank);

  for (int i = 0; i < size; i++) {
    if(rank == exchangeList[i]) {
      sendTo = exchangeList[(i+1)%size];
      receiveFrom = exchangeList[(size+i-1)%size]; //Add size to prevent index before modulus from being negative.
    }
  }
 
  //Replace the nIndividual least fit individuals in the population of the
  //receiving node with the fittest in the sending node.
  for(int i = 0; i < nIndividuals; i++) {
    sendIndividual(population[i], sendTo, req, size, rank);
    receiveIndividual(population[population.size()-1-i], receiveFrom, stat, size, rank);
  }
}



//Sends an individual to another node
void sendIndividual(
  individual& ind, //Individual to be sent
  int toNode, //Destination node
  MPI_Request& req, //Request, for sending via MPI
  int size, //Number of nodes
  int rank //Label of node
) {
  //Comm ID is rank of sender for path, size+rank of sender for fitness
  MPI_Isend(&ind.path[0], ind.path.n_elem, MPI_INTEGER8, toNode, rank, MPI_COMM_WORLD, &req);
  MPI_Isend(&ind.fitness, 1, MPI_REAL8, toNode, size+rank, MPI_COMM_WORLD, &req);
}



//Receives an individual from another node
void receiveIndividual(
  individual& ind, //Output: received individual
  int receiveFrom, //Sender node
  MPI_Status& status, //Status, for receiving via MPI
  int size, //Number of nodes
  int rank //Label of node
) {
  //Comm ID is rank of sender for path, size+rank of sender for fitness
  MPI_Recv(&ind.path[0], ind.path.n_elem, MPI_INTEGER8, receiveFrom, receiveFrom, MPI_COMM_WORLD, &status);
  MPI_Recv(&ind.fitness, 1, MPI_REAL8, receiveFrom, size+receiveFrom, MPI_COMM_WORLD, &status);
}



//Find the best individual among the fittest in each node
individual bestOverallIndividual(
  std::vector<individual>& population, //Population of individuals in the
                                       //current node
  int size, //Number of nodes
  int rank //Label of node
) {
  //THE FUNCTION ASSUMES THE POPULATION HAS BEEN ALREADY SORTED
  //Probably comment this check out if your program has been tested and found to work in order to gain some speed.
  if(!std::is_sorted(population.begin(), population.end())) {
    std::cout << "Error: unsorted population passed to bestOverallIndividual. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }

  //MPI variables, for sending/receiving
  MPI_Request req;
  MPI_Status status;

  //Node 0 receives the fittest individuals for all nodes and compares them
  //to find the fittest overall
  if (rank != 0) {
    sendIndividual(population[0], 0, req, size, rank);
  }

  if (rank == 0) {
    //Fittest individuals for every node
    std::vector<individual> candidates(size);
    //Variable to receive individuals via MPI. Initialize it to some value
    //so we receive on memory that's allocated
    individual tmp = population[0];
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
