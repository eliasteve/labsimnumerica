#include "metropolis.h"

//Code to implement the Metropolis-Hastings algorithm

//Propose move from a uniform distribution for scalar Metropolis
double propose_unif(
    double val, //Current value
    double radius, //Width of the distribution
    Random &rng //Random number generator
) {
  return rng.Rannyu(val-radius, val+radius);
}

//Generates the next point in the Metropolis algorithm (for functions of
//a scalar variables)
valAndAccept gen_next_point(
  double prev, //Previous value
  std::function<double(double)> distro, //Distribution to be sampled
  //Distribution for move proposal
  std::function<double(double, double, Random&)> propose,
  Random &gen, //Random number generator
  double move_width //Width for the move proposal distribution
) {
  valAndAccept res; //Struct for next value and wether the proposal
                    //was accepted
  double next = propose(prev, move_width, gen); //Next value
  double r = distro(next)/distro(prev); //Ratio for point acceptance
  double random_val = gen.Rannyu(); //Random value for determining wether
                                    //to accept the point

  if (random_val < r) {
    res.val = next;
    res.accepted = 1;
  }
  else {
    res.val = prev;
    res.accepted = 0;
  }

  return res;
}


//Generates the next point in the Metropolis algorithm for sampling an even
//function with the same signature as the wavefunction of ex. 8.
//Since the function we are going to sample is made of the superposition
//of two bell curves, to mitigate the risk of only sampling one half of
//the function when the bells get narrow, we add a flip (x -> -x) once every
//set number of moves.
valAndAccept gen_next_point(
  double prev, //Previous point in the Markov chain
  //Probability distribution function to sample from
  std::function<double(double, double, double)> distro,
  //Distribution function for move proposal
  std::function<double(double, double, Random&)> propose,
  Random &gen, //Random number generator
  double move_width, //Width parameter for proposal distribution
  double mu, //Parameter mu (e.g., center of the wavefunction)
  double sigma, //Parameter sigma (e.g., spread of the wavefunction)
  int &flipCount, //Counter for moves before a flip
  int movesBeforeFlip //Number of moves before a flip is performed
) {
  //Struct to store the resulting point and acceptance status
  valAndAccept res;

  //Incorporate flips to sample both sides when sigma is small
  if (flipCount == movesBeforeFlip) {
    flipCount = 0;
    res.val = -prev;
    res.accepted = 1;
  }
  else {
    double next = propose(prev, move_width, gen); //Proposed next point
    //Acceptance ratio
    double r = distro(next, mu, sigma)/distro(prev, mu, sigma);
    double random_val = gen.Rannyu(); //Generated random number

    if (random_val < r) {
      res.val = next;
      res.accepted = 1;
    }
    else {
      res.val = prev;
      res.accepted = 0;
    }
  }

  return res;
}



//Propose move from a uniform distribution for point (3D vector) Metropolis
point propose_unif(
  point val, //Current point
  double cubeHalfSide, //Width of the uniform distribution
  Random &rng //Random number generator
) {
  point proposal = {0, 0, 0}; //Proposed point
  proposal.x = rng.Rannyu(val.x-cubeHalfSide, val.x+cubeHalfSide);
  proposal.y = rng.Rannyu(val.y-cubeHalfSide, val.y+cubeHalfSide);
  proposal.z = rng.Rannyu(val.z-cubeHalfSide, val.z+cubeHalfSide);
  return proposal;
}

//Propose move from a multivariate gaussian distribution for point
//(3D vector) Metropolis
point propose_gauss(
  point val, //Current point
  double width, //Standard deviation of the distribution
  Random &rng //Random number generator
) {
  point proposal = {0, 0, 0};
  proposal.x = rng.Gauss(val.x, width);
  proposal.y = rng.Gauss(val.y, width);
  proposal.z = rng.Gauss(val.z, width);
  return proposal;
}

//Generates the next point in the Metropolis algorithm (for functions
//of a point (3D vector) variable)
pointAndAccept gen_next_point(
  point prev, //Current point in the Metropolis algorithm
  std::function<double(point)> distro, //Distribution to be sampled
  std::function<point(point, double, Random&)> propose, //Distribution for
                                                        //move proposal
  Random &gen, //Random Number Generator
  double cubeHalfSide //Width of the proposal distribution
) {
  pointAndAccept res; //Struct to hold next point and wether proposal was
                      //accepted
  point next = propose(prev, cubeHalfSide, gen); //Proposed point
  double r = distro(next)/distro(prev); //Ratio for point acceptance
  double random_val = gen.Rannyu(); //Random value to decide wether to
                                    //accept point

  if (random_val < r) {
    res.p = next;
    res.accepted = 1;
  }
  else {
    res.p = prev;
    res.accepted = 0;
  }

  return res;
}




