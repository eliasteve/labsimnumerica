#ifndef __POINT__
#define __POINT__

#include <ostream>

//Struct representing a point in 3D and related functions

struct point {
  double x;
  double y;
  double z;
};

std::ostream& operator<< (std::ostream& out, point p);

point operator+ (const point& p1, const point& p2);

point operator/ (const point& p, const double& d);

//Computes the norm^2 of a point
double norm2(point p);


#endif
