#!/bin/bash
#

compiler='gfortran'
rm -f *.o *.mod example *~
$compiler -O3 -c MAH.f90
$compiler -O3 -o example example.f90 MAH.o

