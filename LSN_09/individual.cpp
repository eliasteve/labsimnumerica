#include "individual.h"

individual initIndividual(arma::uvec path) {
  individual res{path, 0};
  return res;
}

double computeFitness(individual& target, arma::mat& distLUT) {
  double fitness = 0;
  for (int i = 0; i < target.path.n_elem; i++) {
    fitness += distLUT(target.path(i), target.path((i+1)%target.path.n_elem));
  }
  return fitness;
}

bool checkConstraint(individual& target) {
  arma::uvec sorted_path = sort(target.path);
  arma::uvec check = arma::conv_to<arma::uvec>::from(arma::linspace(0, target.path.n_elem-1, target.path.n_elem));
  return all(sorted_path == check) and target.path[0]==0;
}

void checkConstraintAndTerminate(individual& target) {
  bool success = checkConstraint(target);
  if (!success) {
      std::cerr << "INVALID INDIVIDUAL GENERATED. PROGRAM TERMINATES." << std::endl;
      exit(EXIT_FAILURE);
  }
}

bool operator<(const individual &ind1, const individual &ind2) {
  return ind1.fitness < ind2.fitness;
}

/*bool operator==(const individual ind1, const individual ind2) {
  return ind1.fitness == ind2.fitness;
}*/

/*bool operator>(const individual ind1, const individual ind2) {
  return ind1.fitness > ind2.fitness;
}*/

void pairPermutation(individual& ind, arma::mat& distLUT) {
  arma::uvec indices = arma::randperm(ind.path.n_elem-1, 2)+1;
  ind.path.swap_rows(indices[0], indices[1]);
  //std::cout << "Swapped indices " << indices[0] << " and " << indices[1] << std::endl;
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

void shift(individual& ind, Random& rng, arma::mat& distLUT) {
  int shift, length, offset;
  length = int(floor(rng.Rannyu(2, ind.path.n_elem))); //ranges from 2 to N-1
  shift = int(floor(rng.Rannyu(1, length)));
  offset = int(floor(rng.Rannyu(0, ind.path.n_elem - length)));
  //std::cout << "shift is " << shift << ", length is " << length << ", offset is " << offset << " and offset + length = " << offset+length << std::endl;
  std::rotate(ind.path.begin()+offset+1, ind.path.begin()+offset+1+(length - shift), ind.path.begin()+offset+1+length);
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

void permuteRanges(individual& ind, Random& rng, arma::mat& distLUT){
  int L = (ind.path.n_elem-1)/2; //Max length of the range. Truncation is the expected behaviour here
  int length = int(floor(rng.Rannyu(1, L)));
  int absOffset = int(floor(rng.Rannyu(0, ind.path.n_elem-1-2*length)));
  int relOffset = int(floor(rng.Rannyu(1, ind.path.n_elem-1-(2*length+absOffset)-1)));
  //std::cout << "[" << 1+absOffset << ", " << 1+absOffset+length << "], [" << 1+absOffset+length+relOffset << ", " << 1+absOffset+2*length+relOffset << "]" << std::endl;
  std::swap_ranges(ind.path.begin()+1+absOffset, ind.path.begin()+1+absOffset+length+1, ind.path.begin()+1+absOffset+length+relOffset); //The 1st range is defined as [first, last)
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

void inversion(individual& ind, Random& rng, arma::mat& distLUT) {
  int length = int(floor(rng.Rannyu(1, ind.path.n_elem)));
  int offset = int(floor(rng.Rannyu(0, ind.path.n_elem - length)));
  //std::cout << "[" << offset+1 << ", " << offset+length+1 << "]" << std::endl;
  std::reverse(ind.path.begin()+offset+1, ind.path.begin()+offset+length+1);
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

int selectionFunction(int M, double r) {
  return floor(M*pow(r, 5));
}

std::array<int, 2> select(int numberOfIndividuals, Random& rng) {
  std::array<int, 2> res;

  res[0] = selectionFunction(numberOfIndividuals, rng.Rannyu());
  res[1] = selectionFunction(numberOfIndividuals-1, rng.Rannyu());
  if(res[1] >= res[0]) res[1]++; //Remove the first selectee from the 2nd selection
  return res;
}

void copyWithOtherParentOrder(individual& child, individual& otherParent, int j) {
  for(int i = 0; i < otherParent.path.n_elem; i++) {
    if(arma::any(child.path == otherParent.path[i])) continue;
    child.path[j] = otherParent.path[i];
    j++;
  }
}

void crossover(individual& parent1, individual& parent2, individual& child1, individual& child2, Random& rng, arma::mat& distLUT) {

  child1.path = arma::uvec(parent1.path.n_elem, arma::fill::zeros);
  child2.path = arma::uvec(parent1.path.n_elem, arma::fill::zeros);
  
  int splitPoint = rng.Rannyu(1, parent1.path.n_elem);
  //std::cout << "Split point is " << splitPoint << std::endl;
  
  //Copy path before the split point
  std::copy(parent1.path.begin(), parent1.path.begin()+splitPoint, child1.path.begin());
  std::copy(parent2.path.begin(), parent2.path.begin()+splitPoint, child2.path.begin());

  //Copy elements after the split point in the order they appear in the other parent
  copyWithOtherParentOrder(child1, parent2, splitPoint);
  copyWithOtherParentOrder(child2, parent1, splitPoint);

  checkConstraintAndTerminate(child1);
  checkConstraintAndTerminate(child2);

  child1.fitness = computeFitness(child1, distLUT);
  child2.fitness = computeFitness(child2, distLUT);
}
