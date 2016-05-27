#ifndef Genotype_HPP
#define Genotype_HPP

class Genotype {

public:
  Genotype(int ind, int loci, int ploidy, std::vector<int> &tot, std::vector<int> &ref, std::vector<double> &err);
  ~Genotype(){};
  std::vector<double> liks, tLiks; // genotype likelihoods: ind x locus x allele,
                                   // and transposed likelihoods: locus x ind x allele

private:
  void t(const int ii, const int ll, const int pp);

};

namespace Freqs {

}

namespace Diseq {

}

#endif
