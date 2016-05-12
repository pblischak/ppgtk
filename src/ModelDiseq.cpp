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
#include "Frequency.hpp"
#include "Phi.hpp"
#include "Genotype.hpp"
#include "DataParser.hpp"

namespace po = boost::program_options;

ModelDiseq::ModelDiseq(int argCount, char **argVar, bool q, bool p){

  quiet = q;
  print = p;
  freq_tune = -999;
  phi_tune = -999;

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "Prints help message.")
    ("num_ind,n", po::value<int>(&nInd)->required(), "REQUIRED: the number of individuals.")
    ("num_loci,l", po::value<int>(&nLoci)->required(), "REQUIRED: the number of loci.")
    ("ploidy_level,p", po::value<int>(&ploidy)->required(), "REQUIRED: the ploidy level of the individuals.")
    ("mcmc_gen,m", po::value<int>(&mcmc_gen)->default_value(5000), "Number of MCMC generations.")
    ("thin", po::value<int>(&thin)->default_value(10), "How often to thin/sample the chain.")
    ("burn_in,b", po::value<int>(&burn)->default_value(2000), "The number of burn-in generations before samples are saved.")
    ("total_reads,t", po::value<std::string>(&totFile)->required(), "REQUIRED: file name with total read counts.")
    ("reference_reads,r", po::value<std::string>(&refFile)->required(), "REQUIRED: file name with reference read counts.")
    ("error_rates,e", po::value<std::string>(&errFile)->required(), "REQUIRED: file name containing per locus read error rates.")
    (",gl", po::value<std::string>(&glFile), "File name with genotype likelihoods.")
    ("freq_tune", po::value<double>(&freq_tune), "Allele frequency tuning parameter for M-H algorithm.")
    ("phi_tune", po::value<double>(&phi_tune), "Allele frequency tuning parameter for M-H algorithm.");

    po::variables_map vm;
    try{
      po::store(po::command_line_parser(argCount, argVar).options(desc).allow_unregistered().run(), vm);
      //po::store(po::parse_command_line(argCount, argVar, desc), vm);


      if(vm.count("help")){
        std::cout << "\n\n" << desc << "\n\n" << std::endl;
        exit(0);
      }

      po::notify(vm);

    }

    // catch any errors in boost argument parsing.
    catch(po::error& e){
      std::cout << "\nerror: " << e.what() << std::endl;
      std::cout << "\n\n" << "Usage: ./ppgtk --model mod -n #taxa -l #loci -p ploidy_level -t total_reads.txt -r ref_reads.txt -e error_file.txt" << std::endl;
      std::cout << "\nFor additional options type: ./pgpsi --model freqs -h\n" << std::endl;
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
  P.getLogLiks(G.tLiks, nInd, nLoci, ploidy);

  Diseq::Phi F(nLoci);
  F.getLogLiks(G.tLiks, nInd, nLoci, ploidy);

  if(freq_tune != -999){
    P.setTune(freq_tune);
  }

  if(phi_tune != -999){
    F.setTune(phi_tune);
  }

  double ran;

  for(int m = 1; m <= mcmc_gen + burn; m++){

    ran = r->uniformRv();

    if(ran <= 0.5){
      P.mhUpdate(G.tLiks, F.vals, nInd, nLoci, ploidy);
    } else {
      F.mhUpdate(G.tLiks, P.vals, nInd, nLoci, ploidy);
    }

    if(m % thin == 0 && m > burn){

      P.writeFrequency(m - burn);
      F.writeFrequency(m - burn);

      if(quiet != 1){
        P.printMeanAcceptRatio();
        F.printMeanAcceptRatio();
      }

    }
  }

}
