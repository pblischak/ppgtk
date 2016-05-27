// Standard headers
#include <vector>
#include <iostream>
//#include <iomanip>

// Boost headers


// Local headers
#include "Genotype.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

Genotype::Genotype(int ind, int loci, int ploidy, std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err){

  double gEpsilon, val;
  int pos_il;
  tLiks.resize(ind*loci*(ploidy+1),0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        gEpsilon = (a / (double) ploidy) * (1 - err[l]) + (1 - (a / (double) ploidy)) * err[l];


        if(tot[i*loci + l] != 0){
          tLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] = r->lnBinomPdf(tot[i*loci + l], ref[i*loci + l], gEpsilon);
          if(l == 0){
            //std::cout << l*ind*(ploidy+1) + i*(ploidy+1) + a << " " << tot[i*loci + l] << "\t" << ref[i*loci + l] << "\t" << tLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] << "\n";
          }
        } else {
          tLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] = -9999.0;
        }
      }
      if(l == 0){
        //std::cout << "----\n";
      }
    }
  }

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      pos_il = i*loci + l;
      for(int a = 0; a <= ploidy; a++){

        gEpsilon = (a / (double) ploidy) * (1 - err[l]) + (1 - (a / (double) ploidy)) * err[l];


        if(tot[pos_il] != 0){
          val = r->lnBinomPdf(tot[pos_il], ref[pos_il], gEpsilon);
          //std::cout << tot[pos_il] << ", " << ref[pos_il] << ", " << val << "\n";
          liks.push_back(val);
        } else {
          liks.push_back(-9999.0);
        }

      }
      //std::cout << "----\n";
    }
  }

  //t(ind, loci, ploidy);

}

void Genotype::t(const int ii, const int ll, const int pp){

  int pos_ila;

  for(int l = 0; l < ll; l++){
    for(int i = 0; i < ii; i++){
      for(int a = 0; a <= pp; a++){

        pos_ila = i*ll*pp + l*pp + a;

        //std::cout << i + 1 << "\t" << l + 1 << "\t" << a << "\n";

        tLiks.push_back(liks[pos_ila]);
        if(l == 0){
          //std::cout << liks[pos_ila] << "\n";
        }

      }
    }
  }
}
