#!/bin/bash

g++ invert.cpp -std=c++14 -O0 -I/home/german/chemfiles/build/install/include -L/home/german/chemfiles/build/install/lib/ -lchemfiles -lnetcdf -o einvert

exit 0
