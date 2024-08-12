#include "misc.h"

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
