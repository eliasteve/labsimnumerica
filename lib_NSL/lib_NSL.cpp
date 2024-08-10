#include "lib_NSL.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"

std::ostream& operator<< (std::ostream& out, point p) {
  return out << "(" << p.x << ", " << p.y << ", " << p.z << ")"; 
}

point operator+ (const point& p1, const point& p2) {
  point result;
  result.x = p1.x + p2.x;
  result.y = p1.y + p2.y;
  result.z = p1.z + p2.z;
  return result;
}

double propose_unif(double val, double radius, Random &rng) {
  return rng.Rannyu(val-radius, val+radius);
}

valAndAccept gen_next_point(
  double prev,
  std::function<double(double)> distro,
  std::function<double(double, double, Random&)> propose,
  Random &gen,
  double move_width
) {
  valAndAccept res;
  double next = propose(prev, move_width, gen);
  double r = distro(next)/distro(prev);
  double random_val = gen.Rannyu();

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

valAndAccept gen_next_point(
  double prev,
  std::function<double(double, double, double)> distro,
  std::function<double(double, double, Random&)> propose,
  Random &gen,
  double move_width,
  double mu,
  double sigma
) {
  valAndAccept res;
  double next = propose(prev, move_width, gen);
  double r = distro(next, mu, sigma)/distro(prev, mu, sigma);
  double random_val = gen.Rannyu();

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


double norm(point p) {return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));}

point propose_unif(point val, double cubeHalfSide, Random &rng) {
  point proposal = {0, 0, 0};
  proposal.x = rng.Rannyu(val.x-cubeHalfSide, val.x+cubeHalfSide);
  proposal.y = rng.Rannyu(val.y-cubeHalfSide, val.y+cubeHalfSide);
  proposal.z = rng.Rannyu(val.z-cubeHalfSide, val.z+cubeHalfSide);
  return proposal;
}

pointAndAccept gen_next_point(
  point prev,
  std::function<double(point)> distro,
  std::function<point(point, double, Random&)> propose,
  Random &gen,
  double cubeHalfSide
) {
  pointAndAccept res;
  point next = propose(prev, cubeHalfSide, gen);
  double r = distro(next)/distro(prev);
  double random_val = gen.Rannyu();

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

//Function to handle computing the mean and error for a certain progressive
//number of blocks, wrting the values to file and updating the relevant
//accumulators
void computeUpdateMeanAndError(int progressiveIndex, double blockMean, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out) {
  // Parameters:
  //
  // progressiveIndex: index of the current block
  // blockMean: mean for the current block
  // meanAccumulator: accumulator for the block means. When the function is
  // called contains the means of all the previous blocks
  // mean2Accumulator: accumulator for the squared of the block means
  // out: output file

  double mean; //Progressive mean

  meanAccumulator += blockMean;
  mean2Accumulator += pow(blockMean, 2);
  mean = meanAccumulator/(progressiveIndex+1);
  out << mean << " " << sqrt((mean2Accumulator/(progressiveIndex+1) - pow(mean, 2))/progressiveIndex) << std::endl;
}


