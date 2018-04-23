#!/bin/bash

g++ invert.cpp -O0 -I/home/german/labo/18/chemfiles/build/install/include -L/home/german/labo/18/chemfiles/build/install/lib -lchemfiles -lnetcdf -o einvert

exit 0
