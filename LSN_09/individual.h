#ifndef __INDIVIDUAL__
#define __INDIVIDUAL__

#include <armadillo>
#include "../random/random.h"
#include <cstdlib> //exit()

struct individual{
    arma::uvec path;
    double fitness;
};

double computeFitness(individual&, arma::mat&);
bool checkConstraint(individual&);
void checkConstraintAndTerminate(individual&);
std::array<int, 2> select(int, Random&);
int selectionFunction(int, double);
void crossover(individual&, individual&, individual&, individual&, Random&, arma::mat&);
void copyWithOtherParentOrder(individual&, individual&, int);

bool operator<(const individual&, const individual&);
//bool operator==(const individual, const individual);
//bool operator>(const individual, const individual);

void pairPermutation(individual&, arma::mat&); 
void shift(individual&, Random&, arma::mat&);
void permuteRanges(individual&, Random&, arma::mat&);
void inversion(individual&, Random&, arma::mat&);
#endif //__INDIVIDUAL__
