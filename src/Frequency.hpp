#ifndef Frequency_HPP
#define Frequency_HPP

class Frequency {
public:
  // Constructor/destructor
  Frequency(int &loci);
  ~Frequency(){};

  // Member functions
  void getLogLiks(std::vector<double> &gLiks, int ind, int loci, int ploidy);
  void mhUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy);
  void brentUpdate();
  void emUpdate();
  void writeFrequency(int &iter);
  void printMeanAcceptRatio();
  void setTune(double &newTune){ tune = newTune; }
  void setOutFile(std::string &newOutFile){ outFile = newOutFile; }

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
