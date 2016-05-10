// System headers
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

// Boost headers

#ifdef _OPENMP
  #include<omp.h>
#endif

// Local headers
#include "Frequency.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

//


Frequency::Frequency(int &loci){

  outFile = "frequencies.txt";
  tune = 0.2;
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

  //#pragma omp parallel for collapse(3)
  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){

        if(a == 0) tmpLik = 0.0;

        if(gLiks[i*loci*ploidy + l*ploidy + a] != -9999.0){
          tmpVal1 = log(r->binomPdf(ploidy, a, vals[l]));
          tmpVal2 = log(gLiks[i*loci*ploidy + l*ploidy + a]);
          //std::cout << tmpVal1 << "," << tmpVal2 << "\t";
          tmpLik += exp(tmpVal1 + tmpVal2);
        } else {
          continue;
        }
        //std::cout << tmpLik << "\n";

        if(a == ploidy) currLogLiks[l] += log(tmpLik);

      }
      //std::cout << "----\n";
    }
  }

}

std::vector<double> Frequency::calcLogLik(std::vector<double> &gLiks, std::vector<int> &tot, std::vector<int> &ref, int ind, int loci, int ploidy, double f){

  std::vector<double> logLiks(loci, 0), indLikVec(ploidy + 1, 0);
  double indLik;

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
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

      indLik = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);
      logLiks[l] += log(indLik);

    }
  }

  return logLiks;

}

std::vector<double> Frequency::calcLogLik(std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err, int ind, int loci, int ploidy, double f){

  double gEpsilon, indLik;
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

        indLik = std::accumulate(indLikVec.begin(), indLikVec.end(), 0.0);

        /*if(l == 0){
          std::cout <<std::setw(10) << std::setprecision(10) << indLik << "\n";
        }*/

        logLiks[l] += log(indLik);

      }
    }

    return logLiks;
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
    //std::cout << vals[l] << "," << newVals[l] << "\n";

  }

  //#pragma omp parallel for collapse(3)
  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){

        if(a == 0) tmpLik = 0.0;

        tmpVal1 = log(r->binomPdf(ploidy, a, newVals[l]));
        tmpVal2 = log(gLiks[i*loci*ploidy + l*ploidy + a]);
        tmpLik +=  exp(tmpVal1 + tmpVal2);

        if(a == ploidy) newLogLiks[l] += log(tmpLik);

      }
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

std::vector<double> Frequency::emUpdate(std::vector<double> &gLiks, int ind, int loci, int ploidy){

  std::vector<double> newLogLiks(currLogLiks.size()), indSum(vals.size(), 0);
  std::vector<double> newVals(vals.size(), -1);
  double tmpVal1, tmpVal2, tmpLik;
  int ploidyCount = 0;

  for(int i = 0; i < ind; i++){
    ploidyCount += ploidy;
    for(int l = 0; l < loci; l++){
      for(int a = 0; a <= ploidy; a++){

        if(a == 0){
          tmpVal1 = 0.0;
          tmpVal2 = 0.0;
        }

        tmpVal1 += a * gLiks[i*loci*ploidy + l*ploidy + a] * r->binomPdf(ploidy, a, vals[l]);
        tmpVal2 += gLiks[i*loci*ploidy + l*ploidy + a] * r->binomPdf(ploidy, a, vals[l]);

      }

      indSum[l] += tmpVal1 / tmpVal2;
      newLogLiks[l] += log(tmpVal2);

    }
  }

  for(int l = 0; l < loci; l++){
    vals[l] = indSum[l] / (double) ploidyCount;
  }

  return(newLogLiks);

}

void Frequency::brentUpdate(){

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
