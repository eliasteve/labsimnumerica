#ifndef __LIB_NSL__
#define __LIB_NSL__

#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include "../random/random.h"

struct valAndAccept {
  double val;
  int accepted;
};

double propose_unif(double val, double radius, Random &rng);

valAndAccept gen_next_point(
  double prev,
  std::function<double(double)> distro,
  std::function<double(double, double, Random&)> propose,
  Random &gen,
  double move_width
);

struct point {
  double x;
  double y;
  double z;
};

std::ostream& operator<< (std::ostream& out, point p);

double norm(point p);

struct pointAndAccept {
  point p;
  int accepted;
};

point propose_unif(point val, double cubeHalfSide, Random &rng);

pointAndAccept gen_next_point(
  point prev,
  std::function<double(point)> distro,
  std::function<point(point, double, Random&)> propose,
  Random &gen,
  double cubeHalfSide
);

//Function to handle computing the mean and error for a certain progressive
//number of blocks, wrting the values to file and updating the relevant
//accumulators
void computeUpdateMeanAndError(int progressiveIndex, double blockValue, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out);

valAndAccept gen_next_point(
  double prev,
  std::function<double(double, double, double)> distro,
  std::function<double(double, double, Random&)> propose,
  Random &gen,
  double move_width,
  double mu,
  double sigma
);

#endif //__LIB_NSL__
