// System headers
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

// Boost headers


// Local headers
#include "Phi.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

Phi::Phi(int &loci){

  tune = 0.05;
  outFile = "phi.txt";
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted.resize(loci, 0);
  nProposals.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran;
  aa = 2, bb = 0.1;

  for(int l = 0; l < loci; l++){
    ran = r->gammaRv(aa, bb);
    vals.push_back(ran);
  }

}

void Phi::getLogLiks(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  int pos_lia;
  double indLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, vals[l]*freqs[l], (1-vals[l])*freqs[l]));

      }

      indLik = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      currLogLiks[l] += log(indLik);
    }
  }

}

std::vector<double> Phi::calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  int pos_lia;
  double indLik;
  std::vector<double> indLikVec(ploidy+1, 0.0), logLiks(loci, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, vals[l]*freqs[l], (1-vals[l])*freqs[l]));

      }

      indLik = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLiks[l] += log(indLik);

    }
  }

  return logLiks;

}

double Phi::calcLogLik(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loc, int ploidy){

  int pos_lia;
  double indLik, logLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){
    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, vals[loc]*freqs[loc], (1-vals[loc])*freqs[loc]));
    }

    indLik = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
    logLik += log(indLik);

  }

  return logLik;

}

void Phi::mhUpdate(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

}

void Phi::writePhi(const int &iter){

  std::ofstream outFileStream;
  outFileStream.open(outFile, std::ios::out | std::ios::app);

  if(outFileStream.is_open()){

    outFileStream << iter << "\t";

    for(int l = 0; l < vals.size(); l++){
      outFileStream << std::setw(8) << std::setprecision(8) << vals[l] << "\t";
    }

    outFileStream << "\n";

  } else {
    std::cout << "Failed to open file: " << outFile << "...\n";
    exit(1);
  }

}
