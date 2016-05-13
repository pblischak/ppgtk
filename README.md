# PPGtk: the polyploid pop-gen toolkit

PPGtk is a C++ program for the analysis of population genomic data collected in polyploids (or mixtures of ploidies). It uses genotype likelihoods to integrate over genotype uncertainty. It relies on the Boost C++ libraries.

<a href="http://www.boost.org/" target="_blank">Boost</a> and <a href="http://openmp.org/wp/" target="_blank">OpenMP</a>.

## Installing Boost

The only Boost library that the PPGtk uses is the Program Options library for command line parsing. Boost is a very large set of libraries for C++ programming and can take a long time to install everything. Since we only need Program Options, use the following commands to build the libraries.

```
cd boost_1_60/
./bootstrap.sh --prefix=/usr/local --with-libraries=program_options
sudo ./b2 install
```

If you want to build all of them, just remove the `--with-libraries=program_options` flag.
