#!/bin/bash

g++ $1 -std=c++14 -O0 -I/home/german/chemfiles/build/install/include -L/home/german/chemfiles/build/install/lib -lchemfiles -lnetcdf -o $2

exit 0
