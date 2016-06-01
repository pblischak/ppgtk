// Standard headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

// Boost headers
#include <boost/program_options.hpp>

// Local headers
#include "ModelAlloSNP.hpp"
#include "Frequency.hpp"
#include "Genotype.hpp"
#include "Theta.hpp"
#include "AncFreq.hpp"
#include "DataParser.hpp"

namespace po = boost::program_options;

ModelAlloSNP::ModelAlloSNP(std::string cFile, bool q, bool p){

  quiet = q;
  print = p;
  std::ifstream cFileStream;
  cFileStream.open(cFile, std::ios::in);

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
    ("num_ind", po::value<int>(&nInd)->required(), "the number of individuals.")
    ("num_loci", po::value<int>(&nLoci)->required(), "the number of loci.")
    ("ploidy1", po::value<int>(&ploidy1)->required(), "the ploidy level of subgenome 1.")
    ("ploidy2", po::value<int>(&ploidy2)->required(), "the ploidy level of subgenome 2 (ploidy2 >= ploidy1).")
    ("mcmc_gen", po::value<int>(&mcmc_gen)->required(), "number of MCMC generations.")
    ("thin", po::value<int>(&thin)->required(), "how often to thin/sample the chain.")
    ("burn", po::value<int>(&burn)->required(), "the number of burn-in generations before samples are saved.")
    ("total_reads", po::value<std::string>(&totFile)->required(), "file name with total read counts.")
    ("reference_reads", po::value<std::string>(&refFile)->required(), "file name with reference read counts.")
    ("error_file", po::value<std::string>(&errFile)->required(), "file name containing per locus read error rates.")
    ("glikelihoods", po::value<std::string>(&glFile)->required(), "file name with genotype likelihoods.")
    ("freq_tune", po::value<double>(&freq_tune)->required(), "allele frequency tuning parameter for M-H algorithm.")
    ("theta_tune", po::value<double>(&theta_tune)->required(), "theta tuning parameter for M-H algorithm.")
    ("anc_tune", po::value<double>(&anc_tune)->required(), "ancestral allele frequency tuning parameter for M-H algorithm.");

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
      std::cout << "\nFor config-file options type: ./ppgtk --model alloSNP -h\n" << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  // catch any other exceptions/errors in standard argument parsing.
  catch(std::exception& e){
    std::cout << std::endl << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

}

}

void ModelAlloSNP::run(){

  int ploidy = ploidy1 + ploidy2;

  bool openmp = 0;

  #ifdef _OPENMP
    openmp = 1;
  #endif

  DataParser data;
  data.getReadData(totFile, refFile, errFile, nInd, nLoci);

  Genotype G(nInd, nLoci, ploidy, data.totMat, data.refMat, data.err);

  AlloSNP::Frequency P(nLoci);
  P.getLogLiks(G.tLiks, nInd, nLoci, ploidy1);
  P.setTune(freq_tune);

  AlloSNP::Theta T;
  T.getLogLiks();
  T.setTune(theta_tune);

  AlloSNP::AncFreq A;
  A.getLogLiks();
  A.setTune(anc_tune);

  int ran;

  if(openmp){

    ran = r->sampleInteger(1,3);

    if(ran == 1){
      P.mhUpdateParallel();
    } else if(ran == 2) {
      T.mhUpdateParallel();
    } else if(ran ==3) {
      A.mhUpdateParallel();
    } else {
      std::cout << "Invalid random integer during M-H algorithm for alloSNP model...\n";
      exit(1);
    }

  } else {

    ran = r->sampleInteger(1,3);

    if(ran == 1){
      P.mhUpdate();
    } else if(ran == 2) {
      T.mhUpdate();
    } else if(ran ==3) {
      A.mhUpdate();
    } else {
      std::cout << "Invalid random integer during M-H algorithm for alloSNP model...\n";
      exit(1);
    }

  }

}
