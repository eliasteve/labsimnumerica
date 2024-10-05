#include "individual.h"

//Compute the fitness of a individual (length of the path, uses the L^1
//norm)
double computeFitness(
  individual& target, //Individual whose fitness is to be computed
  arma::mat& distLUT //Lookup table of distance between pairs of cities
) {
  double fitness = 0; //Fitness of the individual
  for (int i = 0; i < target.path.n_elem; i++) {
    fitness += distLUT(target.path(i), target.path((i+1)%target.path.n_elem));
  }
  return fitness;
}

//Check that every individual vsits each city exactly once, and that the
//first visited city is the one designated with 0
bool checkConstraint(
    individual& target //Individual to be checked
) {
  arma::uvec sortedPath = sort(target.path);
  arma::uvec check = arma::conv_to<arma::uvec>::from(arma::linspace(0, target.path.n_elem-1, target.path.n_elem)); //Vector with all of the cities in
                                                                                                                   //order, if every city is visited
                                                                                                                   //exacty once it should be equal to
                                                                                                                   //the sorted path of the individual
  return all(sortedPath == check) and target.path[0]==0;
}

//Checks the constraint (that each city is visited exactly once and that
//the first city to be visited is number 0): if the check fails print an
//error message and terminate
void checkConstraintAndTerminate(
    individual& target //Individual to be checked
) {
  if (!checkConstraint(target)) {
      std::cerr << "INVALID INDIVIDUAL GENERATED. PROGRAM TERMINATES." << std::endl;
      exit(EXIT_FAILURE);
  }
}

//Overloading of the operator <, to sort individuals in order of ascending
//length
bool operator<(const individual &ind1, const individual &ind2) {
  return ind1.fitness < ind2.fitness;
}

//Pair permutation: mutate an individual by selecting a pait of cities in
//the path with uniform probability and swapping them
void pairPermutation(
  individual& ind, //Individual to be mutated
  arma::mat& distLUT //Lookup table of distances between pairs of cities
) {
  //Indices of the selected cities
  arma::uvec indices = arma::randperm(ind.path.n_elem-1, 2)+1;
  ind.path.swap_rows(indices[0], indices[1]);
  //std::cout << "Swapped indices " << indices[0] << " and " << indices[1] << std::endl;
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

//Shift: mutate an individual by selecting an interval in the path with
//uniform probability (both the initial position and the length) and
//shifting it by an amount selected with uniform probability
void shift(
  individual& ind, //Individual to be mutated
  Random& rng, //Random number generator
  arma::mat& distLUT //Lookup table of distances between pairs of cities
) {
  //Example of individual before the shift:
  //
  //offset     shift
  //  <->      <----->
  // [...|****|@@@@@@@|..]
  //      <---------->
  //      length
  //
  //
  // After the shift:
  // [...|@@@@@@@|****|..]

  int shift, length, offset;
  length = int(floor(rng.Rannyu(2, ind.path.n_elem))); //ranges from 2 to N-1
  shift = int(floor(rng.Rannyu(1, length)));
  offset = int(floor(rng.Rannyu(0, ind.path.n_elem - length)));
  //std::cout << "shift is " << shift << ", length is " << length << ", offset is " << offset << " and offset + length = " << offset+length << std::endl;
  std::rotate(ind.path.begin()+offset+1, ind.path.begin()+offset+1+(length - shift), ind.path.begin()+offset+1+length);
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

//Range permutation: mutate an individual by selecting two sub-intervals of
//equal length in the path and swapping them. All relevant quantities are 
//selected with uniform probability
void permuteRanges(individual& ind, Random& rng, arma::mat& distLUT){
  //Example of individual before swapping:
  //
  // absOffset         relOffset
  // <-------->        <------->
  //[..........|******|.........|@@@@@@|..]
  //            <-//->           <-//->
  //            length
  //
  //After swapping:
  //[..........|@@@@@@|.........|******|..]

  int L = (ind.path.n_elem-1)/2; //Max length of the range. Truncation is the expected behaviour here
  int length = int(floor(rng.Rannyu(1, L)));
  int absOffset = int(floor(rng.Rannyu(0, ind.path.n_elem-1-2*length)));
  int relOffset = int(floor(rng.Rannyu(1, ind.path.n_elem-1-(2*length+absOffset)-1)));
  //std::cout << "[" << 1+absOffset << ", " << 1+absOffset+length << "], [" << 1+absOffset+length+relOffset << ", " << 1+absOffset+2*length+relOffset << "]" << std::endl;
  std::swap_ranges(ind.path.begin()+1+absOffset, ind.path.begin()+1+absOffset+length+1, ind.path.begin()+1+absOffset+length+relOffset); //The 1st range is defined as [first, last)
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

//Inversion: mutate an individual by selecting a sub-interval and reversing
//the order of the elements in it. All relevant quantities are chosen with
//uniform probability
void inversion(individual& ind, Random& rng, arma::mat& distLUT) {
  //Example of individual before inversion:
  //
  //  offset  length
  // <------> <---->
  //[........|234567|..]
  //
  //After inversion:
  //[........|765432|..]

  int length = int(floor(rng.Rannyu(1, ind.path.n_elem)));
  int offset = int(floor(rng.Rannyu(0, ind.path.n_elem - length)));
  //std::cout << "[" << offset+1 << ", " << offset+length+1 << "]" << std::endl;
  std::reverse(ind.path.begin()+offset+1, ind.path.begin()+offset+length+1);
  checkConstraintAndTerminate(ind);
  ind.fitness = computeFitness(ind, distLUT);
}

//Function used to select an individual for crossover
int selectionFunction(int M, double r) {
  return floor(M*pow(r, 5));
}

//Select two individuals for crossover
std::array<int, 2> select(int numberOfIndividuals, Random& rng) {
  std::array<int, 2> res; //Indices of the selected individuals

  res[0] = selectionFunction(numberOfIndividuals, rng.Rannyu());
  res[1] = selectionFunction(numberOfIndividuals-1, rng.Rannyu());
  if(res[1] >= res[0]) res[1]++; //Remove the first selectee from the 2nd selection
  return res;
}

//Function which takes as input a child of the form [....|0000] and fills
//in the missing cities in the order in which they appear in another
//individual
void copyWithOtherParentOrder(
  //Child. Must be in the form [....|0000]
  individual& child,
  //Individual who determines in which order the missing cities will be
  //added to the child
  individual& otherParent,
  int j //Index of the first missing city in the child
) {
  for(int i = 0; i < otherParent.path.n_elem; i++) {
    if(arma::any(child.path == otherParent.path[i])) continue;
    child.path[j] = otherParent.path[i];
    j++;
  }
}

//Crossover operation. Given two parents, prodices two children by splitting
//the path of each parent in two, copying the first part of the path and
//then rearranging the cities in the other part to match the order in which 
//they appear in the other parent
void crossover(
  //Parents
  individual& parent1,
  individual& parent2,
  //First output child. Will have the first part of the path of parent 1
  //and the other cities in the order in which they appear in parent 2
  individual& child1,
  //Second output child. Will have the first part of the path of parent 2
  //and the other cities in the order in which they appear in parent 1
  individual& child2,
  Random& rng, //Random number generator
  arma::mat& distLUT //Lookup table of distances between pairs of cities
) {

  child1.path = arma::uvec(parent1.path.n_elem, arma::fill::zeros);
  child2.path = arma::uvec(parent1.path.n_elem, arma::fill::zeros);
  
  //Index of the first element of the second part of the path (see above)
  int splitPoint = rng.Rannyu(1, parent1.path.n_elem);
  //std::cerr << "Split point is " << splitPoint << std::endl;
  
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

//Write an individual path and its fitness to file
void writeIndividual(
  individual& ind, //Individual to write
  std::string filename //Name of output file
) {
  std::ofstream outf(filename);
  if(!outf.is_open()) {
    std::cerr << "IO error. Program terminates." << std::endl;
    exit(EXIT_FAILURE);
  }
  outf << "Fitness: " << ind.fitness << std::endl << ind.path;
  outf.close();
}




