#!/bin/csh -f

#
# 3D States-Mode HN-Detected Processing.

xyz2pipe -in rec/test%03d.ft1 -x \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 \
| nmrPipe  -fn ZF -auto                       \
| nmrPipe  -fn FT -auto                          \
| nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di              \
| nmrPipe -fn REV -verb \
| nmrPipe  -fn TP \
| nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 \
| nmrPipe  -fn ZF -auto  \
| nmrPipe  -fn FT -alt -verb                       \
| nmrPipe  -fn PS -p0 0  -p1 0.0 -di              \
| nmrPipe  -fn POLY -auto -ord 1 \
| nmrPipe  -fn TP \
| nmrPipe  -fn POLY -auto -ord 1 \
| nmrPipe  -fn ZTP \
| nmrPipe  -fn POLY -auto -ord 1 \
> rec/test.pipe
proj3D.tcl -in rec/test.pipe
pipe2ucsf rec/test.pipe 3D.ucsf

