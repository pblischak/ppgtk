// Standard headers
#include <vector>
#include <numeric>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

// Boost headers
// <none>

// OpenMP header
#ifdef _OPENMP
  #include <omp.h>
#endif

// Local headers
#include "Frequency.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

/***************************************************************

  Member function definitions for Frequency class within the
  Freqs model namespace.

****************************************************************/

Freqs::Frequency::Frequency(int &loci){

  outFile = "frequencies.txt";
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted.resize(loci, 0);
  nProposals.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran = 0.0;
  aa = 0.5, bb = 0.5;

  for(int l = 0; l < loci; l++){
    ran = r->betaRv(aa, bb);
    vals.push_back(ran);
  }

}

void Freqs::Frequency::getLogLiks(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        indLikVec[a] = exp(gLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] + r->lnBinomPdf(ploidy, a, vals[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      currLogLiks[l] += log(indLikSum);

    }
  }

}

std::vector<double> Freqs::Frequency::calcLogLik(std::vector<double> &gLiks, std::vector<int> &tot, std::vector<int> &ref, int ind, int loci, int ploidy, double f){

  std::vector<double> logLiks(loci, 0), indLikVec(ploidy + 1, 0);
  double indLikSum;

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        /*tmpVal1 = log(gLiks[l*ind*ploidy + i*ploidy + a]);
        tmpVal2 = r->lnBinomPdf(ploidy, a, f);
        indLik += exp(tmpVal1 + tmpVal2);*/
        indLikVec[a] = exp(gLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] + r->lnBinomPdf(ploidy, a, f));
        if(l == 0){
          //std::cout << l*ind*ploidy + i*ploidy + a << " " << tot[l*ind + i] << "\t"
          //<< ref[l*ind + i] << "\t" << gLiks[l*ind*ploidy + i*ploidy + a] << "\n";
        }

      }

      if(l == 0){
        //std::cout << "----\n";
      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLiks[l] += log(indLikSum);

    }
  }

  return logLiks;

}

std::vector<double> Freqs::Frequency::calcLogLik(std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err, int ind, int loci, int ploidy, double f){

  double gEpsilon, indLikSum;
  std::vector<double> logLiks(loci, 0), indLikVec(ploidy + 1, 0);

    for(int l = 0; l < loci; l++){
      for(int i = 0; i < ind; i++){

        /*if(l == 0){
          std::cout << tot[l*ind + i] << "\t" << ref[l*ind + i] << "\n";
        }*/

        for(int a = 0; a <= ploidy; a++){

          gEpsilon = (a / (double) ploidy) * (1 - err[l]) + (1 - (a / (double) ploidy)) * err[l];

          indLikVec[a] = exp(r->lnBinomPdf(tot[l*ind + i], ref[l*ind + i], gEpsilon) +
                             r->lnBinomPdf(ploidy, a, f));

          //std::cout << r->lnBinomPdf(tot[l*ind + i], ref[l*ind + i], gEpsilon) << "\t";


          /*if(l == 0){
            std::cout << std::setw(10) << std::setprecision(10) << log(indLikVec[a]) << "\t" << r->lnBinomPdf(tot[l*ind + i], ref[l*ind + i], gEpsilon) << "\t" << r->lnBinomPdf(ploidy, a, f) << "\n";
          }*/


        }

        indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);

        /*if(l == 0){
          std::cout <<std::setw(10) << std::setprecision(10) << indLik << "\n";
        }*/

        logLiks[l] += log(indLikSum);

      }
    }

    return logLiks;

}

double Freqs::Frequency::calcLogLik(std::vector<double> &gLiks, int ind, int loc, int ploidy, double f){

  int pos_lia;
  double indLikSum, logLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){

    if(gLiks[loc*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
      continue;
    }

    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gLiks[pos_lia] + r->lnBinomPdf(ploidy, a, f));

    }

    indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
    logLik += log(indLikSum);

  }

  return logLik;

}

void Freqs::Frequency::mhUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  std::vector<double> indLikVec(ploidy+1, 0), newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1);
  double lnMetropRatio, lnU, indLikSum;

  for(int l = 0; l < vals.size(); l++){

    // propose  new values using normal RV centered at current val with s.d. equal to sqrt(tune).
    // Keep making proposals if the proposed new values is outside the range [0,1].
    while(newVals[l] < 0 || newVals[l] > 1){
      newVals[l] = r->normalRv(vals[l], tune);
    }
    //std::cout << vals[l] << "," << newVals[l] << "\n";

  }

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        indLikVec[a] = exp(gLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] + r->lnBinomPdf(ploidy, a, newVals[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      newLogLiks[l] += log(indLikSum);

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

void Freqs::Frequency::mhUpdateParallel(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  std::vector<double> newVals(loci, -1.0), newLogLiks(loci, 0.0);//, lnMetropRatio(loci, 0.0), lnU(loci, 0.0);
  double lnMetropRatio, lnU;

  #pragma omp parallel private(lnMetropRatio, lnU)
  {

    lnMetropRatio = 0.0;
    lnU = 0.0;

    #pragma omp for
    for(int l = 0; l < loci; l++){

      #pragma omp critical
      {
        while(newVals[l] < 0 || newVals[l] > 1){
          newVals[l] = r->normalRv(vals[l], tune);
        }
      }

      newLogLiks[l] = calcLogLik(gLiks, ind, l, ploidy, newVals[l]);

      //std::cout << "\n";

      lnMetropRatio = (newLogLiks[l] + (aa - 1)*log(newVals[l]) + (bb - 1)*log(1 - newVals[l]))
                      - (currLogLiks[l]  + (aa - 1)*log(vals[l]) + (bb - 1)*log(1 - vals[l]));
      //std::cout << lnMetropRatio << " , " << l << "\n";
      lnU = log(r->uniformRv());

      //#pragma omp critical
      //{
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
    //}

    }
  }

}

void Freqs::Frequency::writeFrequency(const int &iter){

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

void Freqs::Frequency::printMeanAcceptRatio(){

  //double meanAcceptRatio = 0.0;
  std::cout << "Allele frequency acceptance ratio: ";
  for(int l = 0; l < acceptRatio.size(); l++){
     std::cout << acceptRatio[l] << "\t";
  }

  //meanAcceptRatio /= (double) acceptRatio.size();

  std::cout  << "\n";

}

/***************************************************************

  Member function definitions for Frequency class within the
  Diseq model namespace.

****************************************************************/

Diseq::Frequency::Frequency(const int &loci){

  outFile = "frequencies.txt";
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted.resize(loci, 0);
  nProposals.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran = 0.0;
  aa = 0.5, bb = 0.5;

  for(int l = 0; l < loci; l++){
    ran = r->betaRv(aa, bb);
    vals.push_back(ran);
  }

}

void Diseq::Frequency::getLogLiks(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy){

  int pos_lia;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, vals[l]*phi[l], (1-vals[l])*phi[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      currLogLiks[l] += log(indLikSum);

    }
  }

}

std::vector<double> Diseq::Frequency::calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy){

  int pos_lia;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0), logLikVec(loci, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, vals[l]*phi[l], (1-vals[l])*phi[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLikVec[l] += log(indLikSum);

    }
  }

  return logLikVec;

}

double Diseq::Frequency::calcLogLik(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loc, int ploidy, double f){

  int pos_lia;
  double indLikSum, logLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){

    if(gLiks[loc*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
      continue;
    }

    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gLiks[pos_lia] + r->lnBetaBinomPdf(ploidy, a, f*phi[loc], (1-f)*phi[loc]));

    }

    indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
    logLik += log(indLikSum);

  }

  return logLik;

}

void Diseq::Frequency::mhUpdate(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy){

  std::vector<double> indLikVec(ploidy+1, 0.0), newLogLiks(currLogLiks.size());
  std::vector<double> newVals(vals.size(), -1.0);
  double lnMetropRatio, lnU, indLikSum;

  for(int l = 0; l < vals.size(); l++){

    while(newVals[l] < 0 || newVals[l] > 1){
      newVals[l] = r->normalRv(vals[l], tune);
    }

  }

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        indLikVec[a] = exp(gLiks[l*ind*(ploidy+1) + i*(ploidy+1) + a] + r->lnBetaBinomPdf(ploidy, a, newVals[l]*phi[l], (1-newVals[l])*phi[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      newLogLiks[l] += log(indLikSum);

    }
  }

  for(int l = 0; l < loci; l++){

    lnMetropRatio = (newLogLiks[l] + (aa - 1)*log(newVals[l]) + (bb - 1)*log(1 - newVals[l]))
                    - (currLogLiks[l] + (aa - 1)*log(vals[l]) + (bb - 1)*log(vals[l]));

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

void Diseq::Frequency::mhUpdateParallel(std::vector<double> &gLiks, std::vector<double> &phi, int ind, int loci, int ploidy){

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
        while(newVals[l] < 0 || newVals[l] > 1){
          newVals[l] = r->normalRv(vals[l], tune);
        }
      }

      newLogLiks[l] = calcLogLik(gLiks, phi, ind, l, ploidy, newVals[l]);

      lnMetropRatio = (newLogLiks[l] + (aa - 1)*log(newVals[l]) + (bb - 1)*log(1 - newVals[l]))
                      - (currLogLiks[l] + (aa - 1)*log(vals[l]) + (bb - 1)*log(vals[l]));

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

void Diseq::Frequency::writeFrequency(const int &iter){

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

void Diseq::Frequency::printMeanAcceptRatio(){

  //double meanAcceptRatio = 0.0;
  std::cout << "Allele frequency acceptance ratio: ";
  for(int l = 0; l < acceptRatio.size(); l++){
     std::cout << acceptRatio[l] << "\t";
  }

  //meanAcceptRatio /= (double) acceptRatio.size();

  std::cout  << "\n";

}

/****************************************************************

  Frequency class in the AlloSNP namespace.

****************************************************************/

AlloSNP::Frequency::Frequency(const int &loci, std::vector<double> &theta, double anc){

  outFile1 = "frequencies1.txt";
  outFile2 = "frequencies2.txt";
  nRow = loci;
  size = loci;
  currLogLiks.resize(loci, 0.0);
  nAccepted1.resize(loci, 0);
  nAccepted2.resize(loci, 0);
  nProposals1.resize(loci, 0);
  nProposals2.resize(loci, 0);
  acceptRatio.resize(loci, 0.0);
  double ran1 = 0.0, ran2 = 0.0;

  for(int l = 0; l < loci; l++){
    ran1 = r->betaRv(anc * theta[0], (1 - anc) * theta[0]);
    vals1.push_back(ran1);

    ran2 = r->betaRv(anc * theta[1], (1 - anc) * theta[1]);
    vals2.push_back(ran2);
  }


}

void AlloSNP::Frequency::getLogLiks(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2){

  int pos_lia, ploidy = ploidy1 + ploidy2;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnTwoCatPoissBinomPdf(ploidy1, ploidy2, a, vals1[l], vals2[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      currLogLiks[l] += log(indLikSum);

    }

    currLogLiks[l] += r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], vals1[l]) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], vals2[l]);

  }

}

std::vector<double> AlloSNP::Frequency::calcLogLikVec(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2){

  int pos_lia, ploidy = ploidy1 + ploidy2;
  double indLikSum;
  std::vector<double> indLikVec(ploidy+1, 0.0), logLikVec(loci, 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnTwoCatPoissBinomPdf(ploidy1, ploidy2, a, vals1[l], vals2[l]));
      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLikVec[l] += log(indLikSum);

    }

    logLikVec[l] += r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], vals1[l]) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], vals2[l]);

  }

  return logLikVec;

}

double AlloSNP::Frequency::calcLogLik(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loc, int ploidy1, int ploidy2, double f1, double f2){

  int pos_lia, ploidy = ploidy1 + ploidy2;
  double indLikSum, logLik;
  std::vector<double> indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){

    if(gLiks[loc*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
      continue;
    }

    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gLiks[pos_lia] + r->lnTwoCatPoissBinomPdf(ploidy1, ploidy2, a, f1, f2));

    }

    indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
    logLik += log(indLikSum);

  }

  logLik += r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], f1) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], f2);

  return logLik;

}

void AlloSNP::Frequency::mhUpdate(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2){

  int pos_lia, ploidy = ploidy1 + ploidy2;
  std::vector<double> indLikVec(ploidy+1, 0.0), newLogLiks(currLogLiks.size());
  std::vector<double> newVals1(vals1.size(), -1), newVals2(vals2.size(),-1);
  double lnMetropRatio, lnU, indLikSum;

  for(int l = 0; l < loci; l++){

    while(newVals1[l] < 0 || newVals1[l] > 1){
      newVals1[l] = r->normalRv(vals1[l], tune);
    }

    while(newVals2[l] < 0 || newVals2[l] > 1){
      newVals2[l] = r->normalRv(vals2[l], tune);
    }

  }

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){

      if(gLiks[l*ind*(ploidy+1) + i*(ploidy+1)] == BADLIK){
        continue;
      }

      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        indLikVec[a] = exp(gLiks[pos_lia] + r->lnTwoCatPoissBinomPdf(ploidy1, ploidy2, a, newVals1[l], newVals2[l]));

      }

      indLikSum = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      newLogLiks[l] += log(indLikSum);

    }

    //newLogLiks[l] += r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], vals1[l]) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], vals2[l]);

  }

  for(int l = 0; l < loci; l++){

    lnMetropRatio = (newLogLiks[l] + r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], newVals1[l]) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], newVals2[l]))
                    - (currLogLiks[l] + r->lnBetaPdf(anc * theta[0], (1 - anc) * theta[0], vals1[l]) + r->lnBetaPdf(anc * theta[1], (1 - anc) * theta[1], vals2[l]));

    lnU = log(r->uniformRv());

    if(lnU < lnMetropRatio){
      vals1[l] = newVals1[l];
      vals2[l] = newVals2[l];
      currLogLiks[l] = newLogLiks[l];
      nAccepted1[l]++;
      nAccepted2[l]++;
      nProposals1[l]++;
      nProposals2[l]++;
      acceptRatio[l] = nAccepted1[l] / (double) nProposals1[l];
    } else {
      nProposals1[l]++;
      nProposals2[l]++;
      acceptRatio[l] = nAccepted1[l] / (double) nProposals1[l];
    }

  }

}

void AlloSNP::Frequency::mhUpdateParallel(std::vector<double> &gLiks, std::vector<double> &theta, double anc, int ind, int loci, int ploidy1, int ploidy2){

}

void AlloSNP::Frequency::writeFrequency(const int &iter){

  std::ofstream outFileStream1, outFileStream2;
  outFileStream1.open(outFile1, std::ios::out | std::ios::app);
  outFileStream2.open(outFile2, std::ios::out | std::ios::app);

  if(outFileStream1.is_open() && outFileStream2.is_open()){

    outFileStream1 << iter << "\t";
    outFileStream2 << iter << "\t";

    for(int l = 0; l < vals1.size(); l++){
      outFileStream1 << std::setw(8) << std::setprecision(8) << vals1[l] << "\t";
      outFileStream2 << std::setw(8) << std::setprecision(8) << vals2[l] << "\t";
    }

    outFileStream1 << "\n";
    outFileStream2 << "\n";

  } else {
    std::cout << "Failed to open files: " << outFile1 << " and " << outFile2 << "...\n";
    exit(1);
  }

}

void AlloSNP::Frequency::printMeanAcceptRatio(){

}
