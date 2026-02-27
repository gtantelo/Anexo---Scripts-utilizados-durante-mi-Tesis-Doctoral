#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 2560 -dspfvs 21 -grpdly 76  \
  -xN              1536  -yN               256  \
  -xT              696  -yT               128  \
  -xMODE            DQD  -yMODE        Echo-AntiEcho  \
  -xSW         7812.500  -ySW         1701.838  \
  -xOBS         599.803  -yOBS          60.784  \
  -xCAR           4.628  -yCAR         116.931 \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D          Complex  \
  -out ./test.fid -verb -ov

sleep 5

nmrPipe -in test.fid                                  \
| nmrPipe -fn SOL      \
| nmrPipe -fn SP -off 0.45 -end 0.98 -c 1 -pow 2  \
| nmrPipe -fn ZF -auto      \
| nmrPipe -fn FT          \
| nmrPipe -fn PS -p0 -90  -p1 0 -di    \
#| nmrPipe -fn EXT -left -sw \
| nmrPipe -fn EXT -x1 6ppm -xn 11ppm -sw    \
| nmrPipe -fn TP     \
#| nmrPipe -fn LP -fb      \
| nmrPipe -fn SP -off 0.45 -end 0.98 -c 1 -pow 1   \
| nmrPipe -fn ZF -auto  \
| nmrPipe -fn FT         \
| nmrPipe -fn PS -p0 90 -p1 0  -di   \
#| nmrPipe -fn EXT -x1 102ppm -xn 133ppm -sw    \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn TP \
| nmrPipe -fn POLY -auto \
| nmrPipe -fn MULT -r 1000 -inv \
-out test.ft2 -ov -verb 
pipe2ucsf test.ft2  test.ucsf
nmrDraw
