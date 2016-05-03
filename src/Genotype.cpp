// System headers
#include <vector>
//#include <iostream>
//#include <iomanip>

// Boost headers


// Local headers
#include "Genotype.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

Genotype::Genotype(int ind, int loci, int ploidy, std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err){

  double gEpsilon, val;
  int pos_il;

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      pos_il = i*loci + l;
      for(int a = 0; a <= ploidy; a++){

        gEpsilon = (a / (double) ploidy) * (1 - err[l]) + (1 - (a / (double) ploidy)) * err[l];


        if(tot[pos_il] != 0){
          val = r->binomPdf(tot[pos_il], ref[pos_il], gEpsilon);
          liks.push_back(val);
        } else {
          liks.push_back(-9999.0);
        }

      }
    }
  }

}
