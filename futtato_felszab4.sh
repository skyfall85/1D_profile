#!/bin/bash
for i in  100 1000 10000
do
 j=100
   mkdir 'H'$i'd0'$j
   cp parameter_maker.cpp 'H'$i'd0'$j
   cp motion8.cpp 'H'$i'd0'$j

   cd 'H'$i'd0'$j
   
   g++ parameter_maker.cpp -o tm
   time ./tm $i $j
   
   g++ motion8.cpp -o tmo
   time ./tmo
   
   rm motion8.cpp
   rm parameter_maker.cpp
   rm tm 
   rm tmo
   cd ..
done

 