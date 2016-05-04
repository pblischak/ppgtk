// System headers
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// Boost headers
#include <boost/math/tools/minima.hpp>

#ifdef _OPENMP
  #include<omp.h>
#endif

// Local headers
#include "Frequency.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

//

namespace brent = boost::math::tools;

Frequency::Frequency(int &loci){

  outFile = "frequencies.txt";
  tune = 0.25;
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted.resize(loci, 0);
  nProposals.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran, aa = 0.5, bb = 0.5;

  for(int l = 0; l < loci; l++){
    ran = r->uniformRv();
    vals.push_back(ran);
  }

}

void Frequency::getLogLiks(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  double tmpVal1, tmpVal2, tmpLik;

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      tmpLik = 0.0;
      for(int a = 0; a <= ploidy; a++){

        if(gLiks[i*loci*ploidy + l*ploidy + a] != -9999.0){
          tmpVal1 = log(r->binomPdf(ploidy, a, vals[l]));
          tmpVal2 = log(gLiks[i*loci*ploidy + l*ploidy + a]);
          //std::cout << tmpVal1 << "," << tmpVal2 << "\t";
          tmpLik += exp(tmpVal1 + tmpVal2);
        } else {
          continue;
        }
        //std::cout << tmpLik << "\n";
      }
      //std::cout << "----\n";
      currLogLiks[l] += log(tmpLik);
    }
  }

}

void Frequency::mhUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  std::vector<double> newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1);
  double lnMetropRatio, lnU, tmpVal1, tmpVal2, tmpLik;

  for(int l = 0; l < vals.size(); l++){

    // propose  new values using normal RV centered at current val with s.d. equal to sqrt(tune).
    // Keep making proposals if the proposed new values is outside the range [0,1].
    while(newVals[l] < 0 || newVals[l] > 1){
      newVals[l] = r->normalRv(vals[l], tune);
    }

  }

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      tmpLik = 0.0;
      for(int a = 0; a <= ploidy; a++){
        tmpVal1 = log(r->binomPdf(ploidy, a, newVals[l]));
        tmpVal2 = log(gLiks[i*loci*ploidy + l*ploidy + a]);
        tmpLik +=  exp(tmpVal1 + tmpVal2);
      }
      newLogLiks[l] += log(tmpLik);
    }
  }

  for(int l = 0; l < loci; l++){

    lnMetropRatio = (newLogLiks[l] + (aa - 1)*log(newVals[l]) + (bb - 1)*log(1 - newVals[l]))
                    - (currLogLiks[l]  + (aa - 1)*log(vals[l]) + (bb - 1)*log(1 - vals[l]));
    lnU = log(r->uniformRv());
    //std::cout << exp(lnMetropRatio) << "," << exp(lnU) << "\n";

    if(lnU < lnMetropRatio){
      vals[l] = newVals[l];
      currLogLiks[l] = newLogLiks[l];
      nAccepted[l]++;
      nProposals[l]++;
      acceptRatio[l] = nAccepted[l] / (double) nProposals[l];
    } else {
      nProposals[l]++;
      acceptRatio[l] = nAccepted[l] / (double) nProposals[l];
    }

  }

}

void Frequency::brentUpdate(){

}

void emUpdate(){

}

void Frequency::writeFrequency(int &iter){

  std::ofstream outFileStream;
  outFileStream.open(outFile, std::ios::out | std::ios::app);

  if(outFileStream.is_open()){

    outFileStream << iter << "\t";

    for(int l = 0; l < vals.size(); l++){
      outFileStream << vals[l] << "\t";
    }

    outFileStream << "\n";

  } else {
    std::cout << "Failed to open file: " << outFile << "...\n";
    exit(1);
  }

}

void Frequency::printMeanAcceptRatio(){

  //double meanAcceptRatio = 0.0;
  std::cout << "Allele frequency acceptance ratio: ";
  for(int l = 0; l < acceptRatio.size(); l++){
     std::cout << acceptRatio[l] << "\t";
  }

  //meanAcceptRatio /= (double) acceptRatio.size();

  std::cout  << "\n";

}
