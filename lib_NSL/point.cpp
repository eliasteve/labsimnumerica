#include "point.h"

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

point operator/ (const point& p, const double& d) {
  point result;
  result.x = p.x / d;
  result.y = p.y / d;
  result.z = p.z / d;
  return result;
}

//Computes the norm^2 of a point
double norm2(point p) {return pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);}
