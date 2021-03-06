#ifndef ModelFreqs_HPP
#define ModelFreqs_HPP

class ModelFreqs {

public:

  ModelFreqs(std::string cFile, bool q, bool p);
  void run();

private:

  std::string totFile, refFile, errFile, glFile;
  bool print, quiet;
  int nInd, nLoci, ploidy, mcmc_gen, thin, burn;
  double freq_tune;

};

#endif
