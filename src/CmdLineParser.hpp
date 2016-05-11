#ifndef CmdLineParser_HPP
#define CmdLineParser_HPP

#define VERSION "0.1.0"
#define VERSION_DATE "(May 2016)"

class CmdLineParser {
public:
  CmdLineParser(int, char**);
  ~CmdLineParser(){};
  std::string model, totFile, refFile, errFile, glFile;
  int nInd, nLoci, ploidy, mcmc_gen, thin, burn;
  long int seed;
  bool quiet, print;
};

#endif
