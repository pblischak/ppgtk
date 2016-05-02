// System headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

// Local headers
#include "DataParser.hpp"

void DataParser::getReadData(const std::string &totFilename, const std::string &refFilename, const std::string &errFilename, const int &ind, const int &loci){

  int val;
  double err_val;
  std::ifstream totStream(totFilename);
  std::ifstream refStream(refFilename);
  std::ifstream errStream(errFilename);

  if(totStream.is_open()){
    while(totStream >> val){
      totMat.push_back(val);
    }
  } else {
    std::cout << "Could not open total reads file: " << totFilename << "...\n";
    exit(1);
  }

  if(refStream.is_open()){
    while(refStream >> val){
      refMat.push_back(val);
    }
  } else {
    std::cout << "Could not open reference reads file: " << refFilename << "...\n";
    exit(1);
  }

  if(errStream.is_open()){
    while(errStream >> err_val){
      err.push_back(err_val);
    }
  } else {
    std::cout << "Could not open error rates file: " << errFilename << "...\n";
    exit(1);
  }

  checkReadMats(ind, loci);

}

void DataParser::printMat(const int &ind, const int &loci){

  int pos;

  for(int i = 0; i < ind; i++){
    for(int l = 0; l < loci; l++){
      pos = i * loci + l;
      std::cout << refMat[pos] << "," << totMat[pos] << "\t";
    }
    std::cout << "\n";
  }
}

void DataParser::checkReadMats(const int a, const int b){
 int prod = a * b;
  if(totMat.size() != prod || refMat.size() != prod){
    std::cout << "\nThe size of the read matrices is not equal to nInd * nLoci.\n\n";
    exit(EXIT_FAILURE);
  }
}
