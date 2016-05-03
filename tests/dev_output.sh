#!/bin/bash

# Testing output of std::cout statements from initial development versions of PPGtk.
# Need to make sure that reading in data works, and parameter initialization
# in constructors.

ppgtk_mhFreqs -n 25 -l 2 -p 4 -t ../example/total.txt -r ../example/reference.txt -e ../example/error.txt \
  -m 50000 -b 5000 --thin 100
