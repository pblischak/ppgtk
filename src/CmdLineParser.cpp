// System headers
#include <string>
#include <iostream>
#include <stdlib.h>

// Boost headers
#include <boost/program_options.hpp>


// Local headers
#include "CmdLineParser.hpp"

namespace po = boost::program_options;

/*!

*/
CmdLineParser::CmdLineParser(int arg_count, char **arg_var){

  seed = -999;

  try{
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "Prints help message.")
    ("version,v","Print software version information.")
    ("num_ind,n", po::value<int>(&nInd)->required(), "REQUIRED: the number of individuals.")
    ("num_loci,l", po::value<int>(&nLoci)->required(), "REQUIRED: the number of loci.")
    ("ploidy_level,p", po::value<int>(&ploidy)->required(), "REQUIRED: the ploidy level of the individuals.")
    ("mcmc_gen,m", po::value<int>(&mcmc_gen)->default_value(5000), "Number of MCMC generations.")
    ("thin", po::value<int>(&thin)->default_value(10), "How often to thin/sample the chain.")
    ("total_reads,t", po::value<std::string>(&totFile)->required(), "REQUIRED: file name with total read counts.")
    ("reference_reads,r", po::value<std::string>(&refFile)->required(), "REQUIRED: file name with reference read counts.")
    ("error_rates,e", po::value<std::string>(&errFile)->required(), "REQUIRED: file name containing per locus read error rates.")
    ("seed,s", po::value<long int>(&seed), "random number seed.")
    ("quiet,q", "Turn off printing run information to stdout.")
    ("print", "Print updates to screen.");

    po::variables_map vm;
    try{
      po::store(po::parse_command_line(arg_count, arg_var, desc), vm);


      if(vm.count("help")){
        std::cout << "\n\n" << desc << "\n\n" << std::endl;
        exit(0);
      }

      if(vm.count("version")){
        std::cout << "\n\n**********************************************************************************" << std::endl;
        std::cout << "**  This is PPGtk version " << VERSION << " " << VERSION_DATE << std::endl;
        std::cout << "**  For help using the sofware type: ./ppgtk -h" << std::endl;
        std::cout << "**  Documentation can be also found online at http://pblischak.github.io/ppgtk/." << std::endl;
        std::cout << "**********************************************************************************\n\n" << std::endl;
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
      std::cout << "\n\n" << "Usage: ./ppgtk -n #taxa -l #loci -p ploidy_level -t total_reads.txt -r ref_reads.txt" << std::endl;
      std::cout << "\nFor additional options type: ./pgpsi -h\n" << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  // catch any other exceptions/errors in standard argument parsing.
  catch(std::exception& e){
    std::cout << std::endl << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

}
