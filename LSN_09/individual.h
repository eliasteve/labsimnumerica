#ifndef __INDIVIDUAL__
#define __INDIVIDUAL__

#include <armadillo>
#include "../random/random.h"
#include <cstdlib> //exit()

//Struct to store a path and its length (fitness)
struct individual{
    arma::uvec path;
    double fitness;
};

//Compute the fitness of a individual (length of the path, uses the L^1
//norm)
double computeFitness(individual&, arma::mat&);
//Check that every individual vsits each city exactly once, and that the
//first visited city is the one designated with 0
bool checkConstraint(individual&);
//Checks the constraint (that each city is visited exactly once and that
//the first city to be visited is number 0): if the check fails print an
//error message and terminate
void checkConstraintAndTerminate(individual&);
//Select two individuals for crossover
std::array<int, 2> select(int, Random&);
//Function used to select an individual for crossover
int selectionFunction(int, double);
//Crossover operation. Given two parents, prodices two children by splitting
//the path of each parent in two, copying the first part of the path and
//then rearranging the cities in the other part to match the order in which 
//they appear in the other parent
void crossover(individual&, individual&, individual&, individual&, Random&, arma::mat&);
//Function which takes as input a child of the form [....|0000] and fills
//in the missing cities in the order in which they appear in another
//individual
void copyWithOtherParentOrder(individual&, individual&, int);

//Overloading of the operator <, to sort individuals in order of ascending
//length
bool operator<(const individual&, const individual&);

//Pair permutation: mutate an individual by selecting a pair of cities in
//the path with uniform probability and swapping them
void pairPermutation(individual&, arma::mat&); 
//Shift: mutate an individual by selecting an interval in the path with
//uniform probability (both the initial position and the length) and
//shifting it by an amount selected with uniform probability
void shift(individual&, Random&, arma::mat&);
//Range permutation: mutate an individual by selecting two sub-intervals of
//equal length in the path and swapping them. All relevant quantities are 
//selected with uniform probability
void permuteRanges(individual&, Random&, arma::mat&);
//Inversion: mutate an individual by selecting a sub-interval and reversing
//the order of the elements in it. All relevant quantities are chosen with
//uniform probability
void inversion(individual&, Random&, arma::mat&);
#endif //__INDIVIDUAL__
