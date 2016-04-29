// System headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

// Local headers
#include "DataParser.hpp"

void DataParser::getReadData(const std::string &totFilename, const std::string &refFilename, const int &ind, const int &loci){

  int val;
  std::ifstream totStream(totFilename);
  std::ifstream refStream(refFilename);

  while(totStream >> val){
    totMat.push_back(val);
  }

  while(refStream >> val){
    refMat.push_back(val);
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
