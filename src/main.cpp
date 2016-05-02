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

    DataParser data;

    data.getReadData(cmd.totFile, cmd.refFile, cmd.errFile, cmd.nInd, cmd.nLoci);

    data.printMat(cmd.nInd, cmd.nLoci);


    std::cout << "\n";
    for(int l = 0; l < cmd.nInd; l++){
      std::cout << std::setw(10) << std::setprecision(10) << data.err[l] << "\t";
    }
    std::cout << "\n";


    Genotype G(cmd.nInd, cmd.nLoci, cmd.ploidy, data.totMat, data.refMat, data.err);

    Frequency P(cmd.nLoci);
    P.writeFrequency(cmd.nInd);
    P.getLogLiks(G.logLiks, cmd.nInd, cmd.nLoci, cmd.ploidy);



    std::cout << "\n\nAllele frequencies...\n";

    for(int l = 0; l < P.vals.size(); l++){
      std::cout << P.vals[l] << "\n";
    }

    delete r;

    return 0;

}
