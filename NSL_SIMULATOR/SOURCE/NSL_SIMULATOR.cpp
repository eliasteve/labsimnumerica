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
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  ofstream outf("../OUTPUT/output.dat", ios::app);
  outf << "Beginning equilibration... ";
  int i;
  for (i = 0; i < SYS.get_neqsteps(); i++) {
    SYS.step();
  }
  outf << "done (" << i << " steps)." << std::endl;
  outf.close();

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    //std::cout << "Beginning block " << i+1 << std::endl;
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      //std::cout << "Beginning step " << j+1 << std::endl;
      SYS.step();
      //std::cout << "Stepped... ";
      SYS.measure();
      //std::cout << "and measured!" << std::endl;
      if(j%10 == 0){
//        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

  return 0;
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
