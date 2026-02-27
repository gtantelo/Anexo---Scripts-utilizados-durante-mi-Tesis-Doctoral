#!/bin/csh -f

rm -rf xyz #clean up

xyz2pipe -in fid/test%03d.fid -x \
| nmrPipe -fn SOL \
#| nmrPipe  -fn EXT -x1 0 -xn 600 \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5  \
| nmrPipe  -fn ZF -auto                         \
| nmrPipe  -fn FT -verb                             \
| nmrPipe  -fn PS -p0 0 -p1 0.0 -di              \
| nmrPipe  -fn EXT -left -sw           \
#| nmrPipe  -fn EXT -x1 10.5ppm -xn 5.5ppm -sw \
| pipe2xyz -ov -out xyz/test%03d.dat -x

echo
echo Press any key to continue to check phasing with nmrDraw.
echo If you change the phase, rerun ft1.com to check again.
echo Otherwise, correct the phase in the reconstruct.com
echo also consider correcting the EXT extraction in the same way
echo and then run reconstruct.com.
setenv Anykey "$<"
nmrDraw xyz/test001.dat




