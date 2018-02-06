#!/usr/bin/env python

import math
import sys

def GetScaleFactorAndError(ee,see,emu,semu):
  ee/=1000.0
  see/=1000.0
  emu/=1000.0
  semu/=1000.0
  r = (ee/emu)
  err=r*math.sqrt(pow(see/ee,2)+pow(semu/emu,2))
  print 'R_ee,emu =',r,'+/-',err
  return r,err

params = [float(num) for num in sys.argv[1:] ]
if len(params) < 4:
    print 'ERROR: must specify eejj and emujj event yields from MC.'
    print 'Usage: calc_ReeEmu.py eejj err_eejj emujj err_emujj'
    exit(-1)
GetScaleFactorAndError(params[0],params[1],params[2],params[3])

