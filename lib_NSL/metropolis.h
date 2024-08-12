#ifndef __LIB_NSL__
#define __LIB_NSL__

#include <functional>
#include "../random/random.h"
#include "point.h"

//Code to implement the Metropolis-Hastings algorithm

struct pointAndAccept {
  point p;
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


point propose_unif(point val, double cubeHalfSide, Random &rng);

pointAndAccept gen_next_point(
  point prev,
  std::function<double(point)> distro,
  std::function<point(point, double, Random&)> propose,
  Random &gen,
  double cubeHalfSide
);

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
