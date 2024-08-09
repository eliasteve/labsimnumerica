#include "lib_NSL.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"

std::ostream& operator<< (std::ostream& out, point p) {
  return out << "(" << p.x << ", " << p.y << ", " << p.z << ")"; 
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

void computeUpdateMeanAndError(int progressiveIndex, double blockMean, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out) {
  double mean;

  meanAccumulator += blockMean;
  mean2Accumulator += pow(blockMean, 2);
  mean = meanAccumulator/(progressiveIndex+1);
  out << mean << " " << sqrt((mean2Accumulator/(progressiveIndex+1) - pow(mean, 2))/progressiveIndex) << std::endl;
}


