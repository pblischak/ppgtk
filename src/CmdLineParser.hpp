#ifndef CmdLineParser_HPP
#define CmdLineParser_HPP

#define VERSION "0.1.0"
#define VERSION_DATE "(April 2016)"

class CmdLineParser {
public:
  CmdLineParser(int, char**);
  std::string totFile, refFile;
  int nInd, nLoci, ploidy;
  long int seed;
  bool quiet, print;
};

#endif
