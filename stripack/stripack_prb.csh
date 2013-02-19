#!/bin/csh
#
F90 -c -g stripack_prb.f90 >& compiler.txt
if ( $status != 0 ) then
  echo "Errors compiling stripack_prb.f90"
  exit
endif
rm compiler.txt
#
F90 stripack_prb.o -L$HOME/lib/$ARCH -lstripack
if ( $status != 0 ) then
  echo "Errors linking and loading stripack_prb.o"
  exit
endif
rm stripack_prb.o
#
mv a.out stripack_prb
./stripack_prb > stripack_prb_output.txt
if ( $status != 0 ) then
  echo "Errors running stripack_prb"
  exit
endif
rm stripack_prb
#
if ( -e stripack_prb_del.eps ) then
  convert stripack_prb_del.eps stripack_prb_del.png
  rm stripack_prb_del.eps
endif
#
if ( -e stripack_prb_vor.eps ) then
  convert stripack_prb_vor.eps stripack_prb_vor.png
  rm stripack_prb_vor.eps
endif
#
echo "Program output written to stripack_prb_output.txt"
