#ifndef __OPTIMIZATION__
#define __OPTIMIZATION__

#include <iostream>
#include <fstream>
#include <string>
#include "individual.h"
#include "../random/random.h"
#include <cmath>
#include <cstdlib> // exit()

//Run one iteration of random search
void randomSearchIteration(std::vector<individual>&, Random&, arma::mat&);

//Mutate one individual applying one of: pair permutation, shift,
//permutation of ranges, inversion, or no mutation at all according to the
//probabilities passed as argument
void mutate(individual&, Random&, double, double, double, double, arma::mat&);

//Run one iteration of the genetic algorithm
std::vector<individual> geneticAlgoIteration(std::vector<individual>&, Random&, arma::mat&, double, double, double, double, double);

//Sort the population and returns the fitness of the best individual
double sortPopulation(std::vector<individual>&);

//Run a random search (mutations only, no crossover) in order to find the
//solution to the travelling salesman problem
void optimizeRandomSearch(std::vector<individual>&, Random&, arma::mat&, int, int, std::string, std::string);

//Write the matrix of cities to file
void writeCities(arma::mat&, std::string);


//Compute the lookup table of the distance between pairs of cities
void computeDistanceLUT(int, arma::mat&, arma::mat&);

//Compute mean and standard deviation of the best half of a population of
//paths
void computeStatisticsOfBestHalf(std::vector<individual>&, double&, double&);
#endif //__OPTIMIZATION__
