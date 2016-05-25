// System headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Boost headers
#include <boost/program_options.hpp>

// Local headers
#include "ModelDiseq.hpp"
#include "MbRandom.hpp"
#include "main.hpp"
#include "Frequency.hpp"
#include "Phi.hpp"
#include "Genotype.hpp"
#include "DataParser.hpp"

namespace po = boost::program_options;

ModelDiseq::ModelDiseq(std::string cFile, bool q, bool p){

  quiet = q;
  print = p;
  std::ifstream cFileStream;
  cFileStream.open(cFile, std::ios::in);

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
    ("num_ind", po::value<int>(&nInd)->required(), "REQUIRED: the number of individuals.")
    ("num_loci", po::value<int>(&nLoci)->required(), "REQUIRED: the number of loci.")
    ("ploidy", po::value<int>(&ploidy)->required(), "REQUIRED: the ploidy level of the individuals.")
    ("mcmc_gen", po::value<int>(&mcmc_gen)->required(), "Number of MCMC generations.")
    ("thin", po::value<int>(&thin)->required(), "How often to thin/sample the chain.")
    ("burn", po::value<int>(&burn)->required(), "The number of burn-in generations before samples are saved.")
    ("total_reads", po::value<std::string>(&totFile)->required(), "REQUIRED: file name with total read counts.")
    ("reference_reads", po::value<std::string>(&refFile)->required(), "REQUIRED: file name with reference read counts.")
    ("error_file", po::value<std::string>(&errFile)->required(), "REQUIRED: file name containing per locus read error rates.")
    ("glikelihoods", po::value<std::string>(&glFile)->required(), "File name with genotype likelihoods.")
    ("freq_tune", po::value<double>(&freq_tune)->required(), "Allele frequency tuning parameter for M-H algorithm.")
    ("phi_tune", po::value<double>(&phi_tune)->required(), "Phi tuning parameter for M-H algorithm");

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
      std::cout << "\nFor additional options type: ./pgpsi --model diseq -h\n" << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  // catch any other exceptions/errors in standard argument parsing.
  catch(std::exception& e){
    std::cout << std::endl << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

}

void ModelDiseq::run(){

  DataParser data;
  data.getReadData(totFile, refFile, errFile, nInd, nLoci);

  Genotype G(nInd, nLoci, ploidy, data.totMat, data.refMat, data.err);

  Diseq::Frequency P(nLoci);
  Diseq::Phi F(nLoci);

  P.getLogLiks(G.tLiks, F.vals, nInd, nLoci, ploidy);
  F.getLogLiks(G.tLiks, P.vals, nInd, nLoci, ploidy);

  P.setTune(freq_tune);
  F.setTune(phi_tune);

  int ran;

  for(int m = 1; m <= mcmc_gen; m++){

    ran = r->sampleInteger(1,2);

    if(ran == 1){
      P.mhUpdate(G.tLiks, F.vals, nInd, nLoci, ploidy);
    } else if(ran == 2){
      F.mhUpdate(G.tLiks, P.vals, nInd, nLoci, ploidy);
    } else {
      std::cout << "Invalid integer in Diseq model M-H algorithm...\n";
      exit(1);
    }

    if(m % thin == 0 && m > burn){

      P.writeFrequency(m);
      F.writePhi(m);

      if(quiet != 1){
        P.printMeanAcceptRatio();
        F.printMeanAcceptRatio();
      }

    }
  }

}
