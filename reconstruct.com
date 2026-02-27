#!/bin/csh
rm -rf yzx # clean up
rm -rf yzx_ist # clean up
mkdir yzx
mkdir yzx_ist

xyz2pipe -in fid/test%03d.fid -x \
| nmrPipe -fn SOL \
#| nmrPipe  -fn EXT -x1 0 -xn 600 \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5   \
| nmrPipe  -fn ZF -auto                      \
| nmrPipe  -fn FT -verb                             \
| nmrPipe  -fn PS -p0 0 -p1 0.0 -di              \
| nmrPipe  -fn EXT -left -sw           \
#| nmrPipe  -fn EXT -x1 10.5ppm -xn 5.5ppm -sw \
| pipe2xyz -ov -out yzx/test%03d.dat -z

parallel -j 100% './ist.com {} > /dev/null ; echo {}' ::: yzx/test*.dat

#!/bin/csh

xyz2pipe -in yzx_ist/test%03d.phf | phf2pipe -user 1 -xproj xz.ft1 -yproj yz.ft1 \
| pipe2xyz -out rec/test%03d.ft1


