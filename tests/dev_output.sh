#!/bin/bash

# Testing output of std::cout statements from initial development versions of PPGtk.
# Need to make sure that reading in data works, and parameter initialization
# in constructors.

ppgtk --model freqs -n 25 -l 2 -p 4 -t ../example/total.txt -r ../example/reference.txt -e ../example/error.txt \
  -m 5000 -b 500 --thin 10 -q
