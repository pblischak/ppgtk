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
    void mhUpdateParallel(std::vector<double> &gLiks, int ind, int loci, int ploidy);
    std::vector<double> calcLogLik(std::vector<double> &gLiks, std::vector<int> &tot, std::vector<int> &ref, int ind, int loci, int ploidy, double f);
    std::vector<double> calcLogLik(std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err, int ind, int loci, int ploidy, double f);
    double calcLogLik(std::vector<double> &gLiks, int ind, int loc, int ploidy, double f);
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
    double calcLogLik(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loc, int ploidy, double f);
    void mhUpdate(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy);
    void mhUpdateParallel(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy);
    void writeFrequency(const int &iter);
    void printMeanAcceptRatio();
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

namespace AlloSNP {

  class Frequency {
    Frequency(const int &loci, std::vector<double> &theta, double anc);
    ~Frequency(){};

    void getLogLiks(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2);
    std::vector<double> calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2);
    double calcLogLik(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loc, int ploidy1, int ploidy2, double f1, double f2);
    void mhUpdate(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2);
    void mhUpdateParallel(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2);
    void writeFrequency(const int &iter);
    void printMeanAcceptRatio();
    void setTune(const double &newTune){ tune = newTune; }


    std::vector<double> vals1, vals2;

  private:
    double tune, aa, bb;
    int nRow, size;
    std::vector<int> nAccepted1, nProposals1, nAccepted2, nProposals2;
    std::string outFile1, outFile2;
    std::vector<double> currLogLiks, acceptRatio;

  };

}

#endif
