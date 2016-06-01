// Standard headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Boost headers
#include <boost/program_options.hpp>

// Local headers
#include "ModelFreqs.hpp"
#include "Frequency.hpp"
#include "Genotype.hpp"
#include "DataParser.hpp"

namespace po = boost::program_options;

ModelFreqs::ModelFreqs(std::string cFile, bool q, bool p){

  quiet = q;
  print = p;
  std::ifstream cFileStream;
  cFileStream.open(cFile, std::ios::in);

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
    ("num_ind", po::value<int>(&nInd)->required(), "the number of individuals.")
    ("num_loci", po::value<int>(&nLoci)->required(), "the number of loci.")
    ("ploidy", po::value<int>(&ploidy)->required(), "the ploidy level of the individuals.")
    ("mcmc_gen", po::value<int>(&mcmc_gen)->required(), "number of MCMC generations.")
    ("thin", po::value<int>(&thin)->required(), "how often to thin/sample the chain.")
    ("burn", po::value<int>(&burn)->required(), "the number of burn-in generations before samples are saved.")
    ("total_reads", po::value<std::string>(&totFile)->required(), "file name with total read counts.")
    ("reference_reads", po::value<std::string>(&refFile)->required(), "file name with reference read counts.")
    ("error_file", po::value<std::string>(&errFile)->required(), "file name containing per locus read error rates.")
    ("glikelihoods", po::value<std::string>(&glFile)->required(), "file name with genotype likelihoods.")
    ("freq_tune", po::value<double>(&freq_tune)->required(), "allele frequency tuning parameter for M-H algorithm.");

    po::variables_map vm;
    try{
      po::store(po::parse_config_file(cFileStream, desc), vm);
      //po::store(po::parse_command_line(argCount, argVar, desc), vm);

      po::notify(vm);

    }

    // catch any errors in boost argument parsing.
    catch(po::error& e){
      std::cout << "\nerror: " << e.what() << std::endl;
      std::cout << "\n\n" << "Usage: ./ppgtk --model <model-name> -c <config-file>" << std::endl;
      std::cout << "\nFor config-file options type: ./ppgtk --model freqs -h\n" << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  // catch any other exceptions/errors in standard argument parsing.
  catch(std::exception& e){
    std::cout << std::endl << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

}

void ModelFreqs::run(){

  bool openmp=0;

  #ifdef _OPENMP
    openmp=1;
  #endif

  //std::cout << openmp << "\n";

  DataParser data;
  data.getReadData(totFile, refFile, errFile, nInd, nLoci);

  Genotype G(nInd, nLoci, ploidy, data.totMat, data.refMat, data.err);

  Freqs::Frequency P(nLoci);
  P.getLogLiks(G.tLiks, nInd, nLoci, ploidy);

  P.setTune(freq_tune);

  for(int m = 1; m <= mcmc_gen; m++){

    if(openmp){
      P.mhUpdateParallel(G.tLiks, nInd, nLoci, ploidy);
    } else {
      P.mhUpdate(G.tLiks, nInd, nLoci, ploidy);
    }

    if(m % thin == 0 && m > burn){

      P.writeFrequency(m);

      if(quiet != 1){
        P.printMeanAcceptRatio();
      }

    }
  }

}
