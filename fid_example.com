#!/bin/csh

bruk2pipe -in ../ser \
  -bad 0.0 -aswap -DMX -decim 1792 -dspfvs 20 -grpdly 67.9841766357422  \
  -xN              2048  -yN                 4  -zN               818  \
  -xT              1024  -yT                 2  -zT               409  \
  -xMODE            DQD  -yMODE           Real  -zMODE           Real  \
  -xSW        11160.714  -ySW         2919.649  -zSW         4444.933  \
  -xOBS         800.284  -yOBS          81.101  -zOBS         201.266  \
  -xCAR           4.754  -yCAR         119.556  -zCAR         178.683  \
  -xLAB              HN  -yLAB             15N  -zLAB             13C  \
  -ndim               3  -aq2D          States                         \
| nmrPipe -fn MAC -macro $NMRTXT/ranceY.M -noRd -noWr   \
| pipe2xyz -x -out ./fid/data%03d.fid -verb -ov -to 0


