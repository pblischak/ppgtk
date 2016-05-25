// System headers
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>

// Boost headers
#include <boost/program_options.hpp>

// Local headers
#include "MbRandom.hpp"
#include "ModelFreqs.hpp"
#include "ModelDiseq.hpp"
//#include "ModelBetaMix.hpp"
//#include "ModelPopAdmix.hpp"
#include "main.hpp"

namespace po = boost::program_options;

MbRandom *r = new MbRandom;

int main(int argc, char **argv){

  std::string model = "none", configFile;
  long int seed = -999;
  bool print = 0, quiet = 0;

  try{
    po::options_description desc("Global options (for model specific options use `ppgtk --model <model-name> -h`)");
    desc.add_options()
    ("help,h", "Prints help message.")
    ("version,v","Print software version information.")
    ("model,m", po::value<std::string>(&model)->required(), "the model to be run:\n -> freqs\n -> diseq\n -> betaMix\n -> popAdmix")
    ("config,c", po::value<std::string>(&configFile)->required(), "configuration file with model options.")
    ("seed,s", po::value<long int>(&seed), "random number seed.")
    ("quiet,q", "Turn off printing run information to stdout.")
    ("print,p", "Print updates to screen.");

    po::variables_map vm;
    try{
      po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
      //po::store(po::parse_command_line(argc, argv, desc), vm);

      if(vm.count("help") && !vm.count("model")){
        std::cout << "\n\n" << desc << "\n\n" << std::endl;
        exit(0);
      }


      if(vm.count("help") && vm.count("model")) {

        if(vm["model"].as<std::string>() == "freqs") {

          std::cout << "\n\nModel arguments for site allele frequency estimation:\n\n"
                    << "An example config file (freqs.txt) can be found in the `config-files/` folder.\n\n"
                    << "  num_ind             The number of individuals.\n"
                    << "  num_loci            The number of loci.\n"
                    << "  ploidy              The ploidy level of the population.\n"
                    << "  total_reads         The name of the total reads file.\n"
                    << "  reference_reads     The name of the reference reads file.\n"
                    << "  error_file          The name of the per locus error rate file.\n\n"
                    << "  mcmc_gen            The number of MCMC generations (default=5000).\n"
                    << "  burn                The number of burn-in generations (default=1000).\n"
                    << "  thin                How often to sample the MCMC chain (default=10).\n"
                    << "  freq_tune           Allele frequency tuning parameter for the M-H algorithm (default=0.1).\n"
                    << "\n\n";
          exit(0);

        } else if(vm["model"].as<std::string>() == "diseq") {

          std::cout << "\n\nModel arguments for HW disequilibrium:\n\n"
                    << "An example config file (diseq.txt) can be found in the `config-files/` folder.\n\n"
                    << "  num_ind             The number of individuals.\n"
                    << "  num_loci            The number of loci.\n"
                    << "  ploidy              The ploidy level of the population.\n"
                    << "  total_reads         The name of the total reads file.\n"
                    << "  reference_reads     The name of the reference reads file.\n"
                    << "  error_file          The name of the per locus error rate file.\n\n"
                    << "  mcmc_gen            The number of MCMC generations (default=5000).\n"
                    << "  burn                The number of burn-in generations (default=1000).\n"
                    << "  thin                How often to sample the MCMC chain (default=10).\n"
                    << "  freq_tune           Allele frequency tuning parameter for M-H algorithm (default=0.1).\n"
                    << "  phi_tune            Phi tuning parameter for M-H algorithm (default=0.1).\n"
                    << "\n\n";
          exit(0);

        } else if(vm["model"].as<std::string>() == "betaMix") {

          std::cout << "\n\nThis model is still in the works...\n\n";
          exit(0);

        } else if(vm["model"].as<std::string>() == "popAdmix") {

          std::cout << "\n\nThis model is still in the works...\n\n";
          exit(0);

        } else {

          std::cout << "\nError: Invalid model specified (" << vm["model"].as<std::string>()
                    << "). \nPlease choose one of the following: freqs, diseq, betaMix, popAdmix.\n\n";
          exit(0);

        }

      }

      if(vm.count("version")){

        // Obligatory figlet logo...
        std::cout << "\n   ____  ____   ____ _   _    \n"
                  << "  |  _ \\|  _ \\ / ___| |_| | __\n"
                  << "  | |_) | |_) | |  _| __| |/ / \n"
                  << "  |  __/|  __/| |_| | |_|   <  \n"
                  << "  |_|   |_|    \\____|\\__|_|\\_\\\n\n";

        std::cout << "**************************************************************************************" << std::endl;
        std::cout << "**  This is PPGtk version " << VERSION << " " << VERSION_DATE << std::endl;
        std::cout << "**  For help using the sofware type: ./ppgtk -h or ./ppgtk --model <model-name> -h" << std::endl;
        std::cout << "**  Documentation can be also found online at http://pblischak.github.io/ppgtk/." << std::endl;
        std::cout << "**************************************************************************************\n\n" << std::endl;
        exit(0);
      }

      // if the quiet flag is passed via the command line
      // set "quiet"=1 so that no printing to stdout occurs.
      if(vm.count("quiet")){
        quiet = 1;
      }

      if(vm.count("print")){
        print = 1;
      }

      po::notify(vm);

    }

    // catch any errors in boost argument parsing.
    catch(po::error& e){
      std::cout << "\nerror: " << e.what() << std::endl;
      std::cout << "\n\n" << "Usage: ./ppgtk --model <model-name> -c <config-file>" << std::endl;
      std::cout << "\nFor additional options type: ./ppgtk -h\n" << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  // catch any other exceptions/errors in standard argument parsing.
  catch(std::exception& e){
    std::cout << std::endl << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  if(seed != -999){
    r->setSeed(seed);
  }

  if(model == "freqs"){
    ModelFreqs mod(configFile, quiet, print);
    mod.run();
  } else if(model == "diseq"){
    ModelDiseq mod(configFile, quiet, print);
    mod.run();
  } else if(model == "betaMix"){
    //ModelBetaMix mod(configFile, quiet, print);
    //mod.run();
  } else if(model == "popAdmix"){
    //ModelPopAdmix mod(configFile, quiet, print);
    //mod.run();
  } else {
    std::cout << "\"" << model << "\"" << " is not a valid model. Please choose one of the following:\n";
    std::cout << "  -> freqs\n  -> diseq\n  -> betaMix\n  -> popAdmix\n";
  }

  //DataParser data;

  //data.getReadData(cmd.totFile, cmd.refFile, cmd.errFile, cmd.nInd, cmd.nLoci);



  //Genotype G(cmd.nInd, cmd.nLoci, cmd.ploidy, data.totMat, data.refMat, data.err);

  /*std::cout << G.liks.size() << "\t" << cmd.nInd * cmd.nLoci * (cmd.ploidy + 1) << "\n";
  for(int g = 0; g < G.liks.size(); g++){
    std::cout << G.liks[g] << "\n";
  }*/

  //Freqs::Frequency P(cmd.nLoci);
  //Diseq::Phi F(cmd.nLoci);
  //P.getLogLiks(G.tLiks, cmd.nInd, cmd.nLoci, cmd.ploidy);

  /*double fFreq;
  std::vector<double> res1, res2;
  for(int ff = 1; ff <= 200; ff++){
    fFreq = (double) ff / 201.0;
    res1 = P.calcLogLik(data.tTotMat, data.tRefMat, data.err, cmd.nInd, cmd.nLoci, cmd.ploidy, fFreq);
    //std::cout << "\n";
    res2 = P.calcLogLik(G.tLiks, data.tTotMat, data.tRefMat, cmd.nInd, cmd.nLoci, cmd.ploidy, fFreq);
    //std::cout << "\n";


    //std::cout << fFreq << "\t";
    for(int r = 0; r < res1.size(); r++){
      //std::cout << res1[r] + 0.5 * log(fFreq) + 0.5 * log(1 - fFreq) << "\t" << res2[r] + 0.5 * log(fFreq) + 0.5 * log(1 - fFreq) << "\t";
    }
    //std::cout << "\n";
  }*/

  /*for(int m = 1; m <= cmd.mcmc_gen; m++){

    P.mhUpdate(G.tLiks, cmd.nInd, cmd.nLoci, cmd.ploidy);

    if(m % cmd.thin == 0 && m > cmd.burn){

      P.writeFrequency(m);

      if(cmd.quiet != 1){
        P.printMeanAcceptRatio();
      }

    }
  }*/

  delete r;

  return 0;

}
