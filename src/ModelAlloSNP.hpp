#ifndef ModelAlloSNP_HPP
#define ModelAlloSNP_HPP

class ModelAlloSNP {

public:
  ModelAlloSNP(std::string cFile, bool q, bool p);
  run();

private:
  std::string totFile, refFile, errFile, glFile;
  bool print, quiet;
  int nInd, nLoci, ploidy1, ploidy2, mcmc_gen, thin, burn;
  double freq_tune, theta_tune, anc_tune;
  
};

#endif
