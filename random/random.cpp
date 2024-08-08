/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"

using namespace std;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: Random(std::string seedFile, std::string primesFile, int primesLineToRead) {
  // Initializes the generator reading seed and primes from supplied files.
  std::ifstream seed(seedFile);

  if (!seed.is_open()) {
    std::cout << "Unable to open file " << seedFile << ". The constructor terminates without performing any action." << std::endl;
    return;
  }

  std::ifstream primes(primesFile);

  if (!primes.is_open()) {
    std::cout << "Unable to open file " << primesFile << ". The constructor terminates without performing any action." << std::endl;
    return;
  }

  int s[4];
  string checkWord; 

  getline(seed, checkWord, '\t');
  if (checkWord != "RANDOMSEED") { //The first word of a seed file will be RANDOMSEED, check if it is true.
    std::cerr << "Error while reading seed from " << seedFile << ": ";
    std::cerr << "Expected the first word to be RANDOMSEED, found " << checkWord << " instead. The constructor terminates without performing any action." << std::endl;
    return;
  }
  else {
    for (int i=0; i<4; i++) {
      seed >> s[i]; 
    }
  }
  seed.close();

  int p1=0, p2=0;
  for (int lineNumber = 0; lineNumber < primesLineToRead; lineNumber++) {
    primes >> p1 >> p2;
  }
  primes.close();

  SetRandom(s, p1, p2);

  std::cout << "Generator initialized with ";
  std::cout << "seed digits " << s[0] << " " << " " << s[1] << " " << " " << s[2] << " " << " " << s[3] << ", ";
  std::cout << "primes digits " << p1 << " " << p2 << std::endl;
  
  
}

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max)
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random :: Exp(double lambda) {
  // Extracts a value from an exponential distribution, lambda*exp(-lambda*x), x >= 0
  double x = Rannyu();
  return -log(1-x)/lambda;
}

double Random :: CauchyLorentz(double mu, double gamma) {
  // Generates a random number with a Cauchy-Lorentz distribution (pdf is gamma/pi*1/((x-mu)^2+gamma^2))
  double x = Rannyu();
  return mu + gamma*tan(M_PI*(x-0.5));
}

double Random :: AR2_1() {
  // Extracts a value from the distribution 3/2*(1-x^2) with the accept-reject method, see exercise 2.1
  double x, y;
  while(1) {
    x = Rannyu();
    y = Rannyu(0, 1.6);

    if( y <= 3*(1-pow(x, 2))/2 ) {
      return x;
    }
  }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
