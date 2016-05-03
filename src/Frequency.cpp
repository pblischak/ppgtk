// System headers
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// Boost headers
#include <boost/math/tools/minima.hpp>

// Local headers
#include "Frequency.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

//

namespace brent = boost::math::tools;

Frequency::Frequency(int &loci){

  outFile = "frequencies.txt";
  tune = 0.008;
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

  double tmpLik;

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){

        if(gLiks[i*loci*ploidy + l*ploidy + a] != -9999.0){
          tmpLik += gLiks[i*loci*ploidy + l*ploidy + a] * r->binomPdf(ploidy, a, vals[l]);
        } else {
          continue;
        }

      }
      currLogLiks[l] += log(tmpLik);
    }
  }

}

void Frequency::mhUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  std::vector<double> newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1);
  double lnMetropRatio, lnU, tmpBinom, tmpBeta, tmpLik;

  for(int l = 0; l < vals.size(); l++){

    // propose  new values using normal RV centered at current val with s.d. equal to sqrt(tune).
    // Keep making proposals if the proposed new values is outside the range [0,1].
    while(newVals[l] < 0 || newVals[l] > 1){
      newVals[l] = r->normalRv(vals[l], tune);
    }

  }

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){
        tmpLik += gLiks[i*loci*ploidy + l*ploidy + a] * r->binomPdf(ploidy, a, newVals[l]);
      }
      newLogLiks[l] += log(tmpLik);
    }
  }

  for(int l = 0; l < loci; l++){

    lnMetropRatio = (newLogLiks[l] + (aa - 1)*log(newVals[l]) + (bb - 1)*log(1 - newVals[l]))
                    - (currLogLiks[l]  + (aa - 1)*log(vals[l]) + (bb - 1)*log(1 - vals[l]));
    lnU = log(r->uniformRv());

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

  double meanAcceptRatio = 0.0;

  for(int l = 0; l < acceptRatio.size(); l++){
    meanAcceptRatio += acceptRatio[l];
  }

  meanAcceptRatio /= (double) acceptRatio.size();

  std::cout << "Allele frequency acceptance ratio (mean): " << meanAcceptRatio << "\n";

}
