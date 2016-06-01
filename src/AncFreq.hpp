#ifndef AncFreq_HPP
#define AncFreq_HPP

namespace AlloSNP {

  class AncFreq {

  public:
    AncFreq();
    void setPrior(const double &alph, const double &bet){ aa = alph; bb = bet; }
    void setTune(const double &newTune){ tune = newTune; }
    mhUpdate();
    mhUpdateParallel();

  private:
    double aa, bb, tune;
    int nAccepted, nProposals;
    std::string outFile;
    double currLogLiks, acceptRatio;

  };

}

#endif
