#ifndef Theta_HPP
#define Theta_HPP

namespace AlloSNP {

  class Theta {

  public:
    Theta();
    void setPrior(const double &alph, const double &bet){ aa = alph; bb = bet; }
    void setTune(const double &newTune){ tune = newTune; }
    void mhUpdate();
    void mhUpdateParallel();

  private:
    double tune, aa, bb;
  };

}

#endif
