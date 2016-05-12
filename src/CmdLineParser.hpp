#ifndef CmdLineParser_HPP
#define CmdLineParser_HPP



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
