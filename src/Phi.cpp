// Standard headers
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

// Boost headers
// <none>

// OpenMP header
#ifdef _OPENMP
  #include <omp.h>
#endif

// Local headers
#include "Phi.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

/***************************************************************

  Defining Phi member functions in the Diseq namespace

****************************************************************/

Diseq::Phi::Phi(int &loci){

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

void Diseq::Phi::getLogLiks(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  int pos_lia;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, freqs[l]*vals[l], (1-freqs[l])*vals[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      currLogLiks[l] += log(indLikSum);

    }
  }

}

std::vector<double> Diseq::Phi::calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  int pos_lia;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0), logLikVec(loci, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, freqs[l]*vals[l], (1-freqs[l])*vals[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLikVec[l] += log(indLikSum);

    }
  }

  return logLikVec;

}

double Diseq::Phi::calcLogLik(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loc, int ploidy, double p){

  int pos_lia;
  double indLikSum, logLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){
    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, freqs[loc]*p, (1-freqs[loc])*p));

    }

    indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
    logLik += log(indLikSum);

  }

  return logLik;

}

void Diseq::Phi::mhUpdate(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  std::vector<double> indLikVec(ploidy+1, 0.0), newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1.0);
  double lnMetropRatio, lnU, indLikSum;

  for(int l = 0; l < vals.size(); l++){

    while(newVals[l] < 0){
      newVals[l] = r->normalRv(vals[l], tune);
    }

    if(l == 0){
      //std::cout << vals[l] << " , " << r->normalRv(vals[l], tune) << "\n";
    }

  }

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        indLikVec[a] = exp(gLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] + r->lnBetaBinomPdf(ploidy, a, freqs[l]*newVals[l], (1-freqs[l])*newVals[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      newLogLiks[l] += log(indLikSum);

    }
  }

  for(int l = 0; l < loci; l++){

    lnMetropRatio = (newLogLiks[l] + r->lnGammaPdf(aa, bb, newVals[l]))
                    - (currLogLiks[l] + r->lnGammaPdf(aa, bb, vals[l]));
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

void Diseq::Phi::mhUpdateParallel(std::vector<double> &gLiks, std::vector<double> &freqs, int ind, int loci, int ploidy){

  std::vector<double> newVals(loci, -1.0), newLogLiks(loci, 0.0);
  double lnMetropRatio, lnU;

  #pragma omp parallel private(lnMetropRatio, lnU)
  {

    lnMetropRatio = 0.0;
    lnU = 0.0;

    #pragma omp for
    for(int l = 0; l < loci; l++){

      #pragma omp critical
      {
        while(newVals[l] < 0){
          newVals[l] = r->normalRv(vals[l], tune);
        }
      }

      newLogLiks[l] = calcLogLik(gLiks, freqs, ind, l, ploidy, newVals[l]);

      lnMetropRatio = (newLogLiks[l] + r->lnGammaPdf(aa, bb, newVals[l]))
                      - (currLogLiks[l] + r->lnGammaPdf(aa, bb, vals[l]));
                      
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

}

void Diseq::Phi::writePhi(const int &iter){

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

void Diseq::Phi::printMeanAcceptRatio(){

  std::cout << "Phi acceptance ratio: ";
  for(int l = 0; l < acceptRatio.size(); l++){
     std::cout << acceptRatio[l] << "\t";
  }

  //meanAcceptRatio /= (double) acceptRatio.size();

  std::cout  << "\n";

}

/***************************************************************

  Defining Phi member functions in the BetaMix namespace.

****************************************************************/

BetaMix::Phi::Phi(const int &loci){

}

/***************************************************************

  Defining Phi member functions in the PopAdmix namespace.

****************************************************************/

PopAdmix::Phi::Phi(const int &loci){

}
