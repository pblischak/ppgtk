#ifndef ModelDiseq_HPP
#define ModelDiseq_HPP

class ModelDiseq {

public:

  ModelDiseq(std::string cFile, bool q, bool p);
  void run();

private:

  std::string totFile, refFile, errFile, glFile;
  bool print, quiet;
  int nInd, nLoci, ploidy, mcmc_gen, thin, burn;
  double freq_tune, phi_tune;

};

#endif
