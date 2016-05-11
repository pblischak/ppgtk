#ifndef Phi_HPP
#define Phi_HPP

class Phi {
public:
  Phi(int &loci);
  ~Phi(){};

  // Member functions
  void getLogLiks(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy);
  std::vector<double> calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy);
  double calcLogLik(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loc, int ploidy);
  void mhUpdate(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy);
  void writePhi(const int &iter);
  void setTune(const double &newTune){ tune = newTune; }
  void setOutFile(const std::string &newOutFile){ outFile = newOutFile; }

  // Member variables
  std::vector<double> vals;

private:
  double tune, aa, bb;
  int nRow, size;
  std::vector<int> nAccepted, nProposals;
  std::string outFile;
  std::vector<double> currLogLiks, acceptRatio;

};

#endif
