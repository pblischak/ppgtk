// System headers
#include <iostream>
#include <cstdlib>

// Boost headers

// OpenMP
#ifdef _OPENMP
  #include <omp.h>
#endif

// Local headers
#include "main.hpp"
#include "CmdLineParser.hpp"
#include "DataParser.hpp"

MbRandom *r = new MbRandom;

int main(int argc, char **argv){

    std::cout << "\n\nThis is PPGtk version " << VERSION << " " << VERSION_DATE << "\n\n";

    CmdLineParser cmd(argc, argv);

    // Reset the random number seed if it is provided as an argument.
    // It is set to -999 in the constructor and we only change it here
    // if the user specifies a new seed.
    if(cmd.seed != -999){
      r->setSeed(cmd.seed);
    }

    DataParser data;

    data.getReadData(cmd.totFile, cmd.refFile, cmd.nInd, cmd.nLoci);

    data.printMat(cmd.nInd, cmd.nLoci);

    double s = 5.0;
    double t = 2.0;

    # pragma omp parallel for ordered
    for(unsigned i = 0; i < 10; i++){
      #pragma omp ordered
      std::cout << i << "\t" << r->gammaRv(s, t) << "\n";
    }

    delete r;

    return 0;

}
