// System headers
#include <vector>

// Boost headers


// Local headers
#include "Genotype.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

Genotype::Genotype(int ind, int loci, int ploidy, std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err){

  double gEpsilon, lnVal;
  int pos_il;

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      pos_il = i*loci + l;
      for(int a = 0; a <= ploidy; a++){

        gEpsilon = (a / (double) ploidy) * (1 - err[l]) + (1 - (a / (double) ploidy)) * err[l];

        if(tot[pos_il] != 0){
          lnVal = r->lnBetaPdf(tot[pos_il], ref[pos_il], gEpsilon);
          logLiks.push_back(lnVal);
        } else {
          logLiks.push_back(-999);
        }

      }
    }
  }

}
