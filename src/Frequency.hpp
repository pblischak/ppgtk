#ifndef Frequency_HPP
#define Frequency_HPP

namespace Freqs {

  class Frequency {
  public:
    // Constructor/destructor
    Frequency(int &loci);
    ~Frequency(){};

    // Member functions
    void getLogLiks(std::vector<double> &gLiks, int ind, int loci, int ploidy);
    void mhUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy);
    std::vector<double> calcLogLik(std::vector<double> &gLiks, std::vector<int> &tot, std::vector<int> &ref, int ind, int loci, int ploidy, double f);
    std::vector<double> calcLogLik(std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err, int ind, int loci, int ploidy, double f);
    //void brentUpdate();
    //std::vector<double> emUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy);
    void writeFrequency(const int &iter);
    void printMeanAcceptRatio();
    void setTune(const double &newTune){ tune = newTune; }
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

}

namespace Diseq {

  class Frequency {

  public:
    Frequency(const int &loci);
    ~Frequency(){};

    // Member functions
    void getLogLiks(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy);
    std::vector<double> calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy);
    double calcLogLik(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loc, int ploidy);
    void mhUpdate(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy);
    void writeFrequency(const int &iter);
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

}

#endif
