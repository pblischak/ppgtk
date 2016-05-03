// System headers
#include <iostream>
#include <iomanip>
#include <cstdlib>

// Boost headers

// OpenMP
#ifdef _OPENMP
  #include <omp.h>
#endif

// Local headers
#include "MbRandom.hpp"
#include "CmdLineParser.hpp"
#include "DataParser.hpp"
#include "Frequency.hpp"
#include "Genotype.hpp"
#include "main.hpp"

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

    int zero = 0;

    DataParser data;

    data.getReadData(cmd.totFile, cmd.refFile, cmd.errFile, cmd.nInd, cmd.nLoci);

    data.printMat(cmd.nInd, cmd.nLoci);


    std::cout << "\n";
    for(int l = 0; l < cmd.nLoci; l++){
      std::cout << std::setw(10) << std::setprecision(10) << data.err[l] << "\t";
    }
    std::cout << "\n";


    Genotype G(cmd.nInd, cmd.nLoci, cmd.ploidy, data.totMat, data.refMat, data.err);

    Frequency P(cmd.nLoci);
    P.writeFrequency(zero);
    P.getLogLiks(G.liks, cmd.nInd, cmd.nLoci, cmd.ploidy);

    for(int m = 1; m <= cmd.mcmc_gen; m++){
      P.mhUpdate(G.liks, cmd.nInd, cmd.nLoci, cmd.ploidy);
      P.writeFrequency(m);
      if(m % cmd.thin == 0){
        P.printMeanAcceptRatio();
      }
    }

    delete r;

    return 0;

}
