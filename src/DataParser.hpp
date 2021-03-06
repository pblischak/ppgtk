#ifndef DataParser_HPP
#define DataParser_HPP

class DataParser {
public:
  DataParser(){};
  ~DataParser(){};
  void getReadData(const std::string &totFilename, const std::string &refFilename, const std::string &errFilename, const int &ind, const int &loci); /*!< . */
  void checkReadMats(const int, const int);
  void printMat(const int &ind, const int &loci);
  std::vector<int> totMat;
  std::vector<int> tTotMat; // transposed totMat: loci x ind
  std::vector<int> refMat;
  std::vector<int> tRefMat; // transposed refMat: loci x ind
  std::vector<double> err;
private:
  void t(const int a, const int b);
};

#endif
