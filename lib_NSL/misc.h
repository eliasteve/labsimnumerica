#ifndef __MISC__
#define __MISC__

#include <cmath>
#include <fstream>

//Function to handle computing the mean and error for a certain progressive
//number of blocks, wrting the values to file and updating the relevant
//accumulators
void computeUpdateMeanAndError(int progressiveIndex, double blockValue, double &meanAccumulator, double &mean2Accumulator, std::ofstream &out);

#endif
