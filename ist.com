#!/bin/csh 

set F = $1

set in = $F:t
set out = $F:t:r.phf

echo $in $out 

hmsIST -dim 2 -incr 1 -autoN 1 -user 1  \
    -itr 400 -verb 0 -ref 0 -vlist nuslist \
    < ./yzx/${in} >! ./yzx_ist/${out}

