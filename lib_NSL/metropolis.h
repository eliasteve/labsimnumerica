#ifndef __LIB_NSL__
#define __LIB_NSL__

#include <functional>
#include "../random/random.h"
#include "point.h"

//Code to implement the Metropolis-Hastings algorithm

//Struct to hold the generated scalar value for a Metropolis step
//and wether the proposed value was accepted or not
struct valAndAccept {
  double val;
  int accepted;
};

//Struct to hold the generated point (3d vector) value for a Metropolis step
//and wether the proposed value was accepted or not
struct pointAndAccept {
  point p;
  int accepted;
};

//Propose move from a uniform distribution for scalar Metropolis
double propose_unif(double val, double radius, Random &rng);

//Generates the next point in the Metropolis algorithm (for functions
//of a scalar variable)
valAndAccept gen_next_point(
  double prev,
  std::function<double(double)> distro,
  std::function<double(double, double, Random&)> propose,
  Random &gen,
  double move_width
);


//Propose move from a uniform distribution for point (3D vector) Metropolis
point propose_unif(point val, double cubeHalfSide, Random &rng);

//Propose move from a multivariate gaussian distribution for point
//(3D vector) Metropolis
point propose_gauss(point val, double width, Random &rng);

//Generates the next point in the Metropolis algorithm (for functions
//of a point (3D vector) variable)
pointAndAccept gen_next_point(
  point prev,
  std::function<double(point)> distro,
  std::function<point(point, double, Random&)> propose,
  Random &gen,
  double cubeHalfSide
);

//lezione 8
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
