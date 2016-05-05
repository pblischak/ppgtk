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

    //std::cout << "\n\nThis is PPGtk version " << VERSION << " " << VERSION_DATE << "\n\n";

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

    //data.printMat(cmd.nInd, cmd.nLoci);

    Genotype G(cmd.nInd, cmd.nLoci, cmd.ploidy, data.totMat, data.refMat, data.err);

    /*std::cout << G.liks.size() << "\t" << cmd.nInd * cmd.nLoci * (cmd.ploidy + 1) << "\n";
    for(int g = 0; g < G.liks.size(); g++){
      std::cout << G.liks[g] << "\n";
    }*/

    Frequency P(cmd.nLoci);
    P.getLogLiks(G.liks, cmd.nInd, cmd.nLoci, cmd.ploidy);

    double fFreq;
    std::vector<double> res;
    for(int ff = 1; ff <= 200; ff++){
      fFreq = (double) ff / 201.0;
      res = P.calcLogLik(data.tTotMat, data.tRefMat, data.err, cmd.nInd, cmd.nLoci, cmd.ploidy, fFreq);

      std::cout << fFreq << "\t";
      for(int r = 0; r < res.size(); r++){
        std::cout << res[r] << "\t";
      }
      std::cout << "\n";
    }

    for(int m = 1; m <= cmd.mcmc_gen; m++){
      //P.mhUpdate(G.liks, cmd.nInd, cmd.nLoci, cmd.ploidy);
      if(m % cmd.thin == 0 && m > cmd.burn){
        //P.writeFrequency(m);
        //P.printMeanAcceptRatio();
      }
    }

    delete r;

    return 0;

}
