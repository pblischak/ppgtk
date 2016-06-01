# PPGtk: the polyploid pop-gen toolkit

PPGtk is a C++ program for the analysis of population genomic data collected in polyploids (or mixtures of ploidies). It uses genotype likelihoods to integrate over genotype uncertainty. It relies on the Boost C++ libraries.

[Boost](http://www.boost.org/) and [OpenMP](http://openmp.org/wp/).

## Installing Boost

The only Boost library that the PPGtk uses is the Program Options library for command line parsing. Boost is a very large set of libraries for C++ programming and can take a long time to install everything. Since we only need Program Options, use the following commands to build the libraries.

```
cd boost_1_60/
./bootstrap.sh --prefix=/usr/local --with-libraries=program_options
sudo ./b2 install
```

If you want to build all of them, just remove the `--with-libraries=program_options` flag.

## Setting the number of threads

You can set the number of threads that you want to use with the `OMP_NUM_THREADS` environmental variable in your Bash shell.

```
export OMP_NUM_THREADS=4

# now running with 4 threads
ppgtk --model freqs --config freqs.txt
```
