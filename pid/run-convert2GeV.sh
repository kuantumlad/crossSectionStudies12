#!/bin/bash

# works with root 6.12

g++ -o convert2GeV filter_clas12Main.cxx `root-config --cflags --glibs`


# $1 input file
# $2 output file
# $3 run number 

./convert2GeV $1 $2 $3
