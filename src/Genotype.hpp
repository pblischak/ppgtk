#ifndef Genotype_HPP
#define Genotype_HPP

class Genotype {
public:
  Genotype(int ind, int loci, int ploidy, std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err);
  ~Genotype(){};
  std::vector<double> liks;
};

#endif
