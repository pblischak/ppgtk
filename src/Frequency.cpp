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
  tune = 0.1;
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted.resize(loci, 0);
  nProposals.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran, alpha = 0.5, beta = 0.5;

  for(int l = 0; l < loci; l++){
    ran = r->uniformRv();
    vals.push_back(ran);
  }

}

void Frequency::getLogLiks(std::vector<double> &gLogLiks, int ind, int loci, int ploidy){

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){

        if(gLogLiks[i*loci*ploidy + l*ploidy + a] != -999){
          currLogLiks[l] += gLogLiks[i*loci*ploidy + l*ploidy + a] + r->lnBinomPdf(ploidy, a, vals[l]) + r->lnBetaPdf(alpha, beta, vals[l]);
        } else {
          continue;
        }
        
      }
    }
  }

}

void Frequency::mhUpdate(std::vector<double> &gLogLiks, int ind, int loci, int ploidy){

  std::vector<double> newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1);
  double lnMetropRatio, lnU;

  for(int l = 0; l < vals.size(); l++){

    // propose  new values using normal RV centered at current val with s.d. equal to sqrt(tune).
    // Keep making proposals if the proposed new values is outside the range [0,1].
    while(newVals[l] < 0 || newVals[l] > 1){
      newVals[l] = r->normalRv(vals[l], sqrt(tune));
    }

  }

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){
        newLogLiks[l] += gLogLiks[i*loci*ploidy + l*ploidy + a] + r->lnBinomPdf(ploidy, a, newVals[l]) + r->lnBetaPdf(alpha, beta, newVals[l]);
      }
    }
  }

  for(int l = 0; l < loci; l++){
    lnMetropRatio = newLogLiks[l] - currLogLiks[l];
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
