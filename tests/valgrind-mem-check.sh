#!/bin/bash

# Using valgrind to check for memory leaks/errors

valgrind ppgtk --model freqs -n 25 -l 2 -p 6 -t ../example/total.txt -r ../example/reference.txt -e ../example/error.txt \
  -m 5000 -b 500 --thin 10 -q
