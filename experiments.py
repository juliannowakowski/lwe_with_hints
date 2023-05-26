from lwe_with_hints import *

import numpy as np
from random import randrange

import argparse
from multiprocessing import Pool
import time
import traceback


def parseArguments():
  parser = argparse.ArgumentParser()

  parser.add_argument("scheme", type=str)
  parser.add_argument("-hints", type=str, default="0:0:1", help="String of the form \"min:max:step\" specifying the number of hints to integrate.")
  parser.add_argument("-trials", type=int, default=1, help="Number of instances to generate.")
  parser.add_argument("-file", type=str, default="output_experiments.txt", help="Output is written into file with name -file.")
  parser.add_argument("-hints_centered", action="store_true", help="If set, then hints are drawn from {-(q-1)/2,...,(q-1)/2}^n. Otherwise from {0,...,q-1}^n.")
  parser.add_argument("-modular", action="store_true", help="If set, then generate mod-q hints. Otherwise perfect hints.")
  parser.add_argument("-verbose", action="store_true")

  args, unknown = parser.parse_known_args()

  if len(unknown) > 0:
    print("Unknown arguments " + str(unknown) + " will be ignored.")

  return vars(args)

def generateHints(s, q, k, centered):
  n = len(s)
  
  V = []
  L = []

  for i in range(k):
    if centered:
      v = np.array([ randrange( -int((q-1)/2), int((q+1)/2) ) for _ in range(n) ])
    else:
      v = np.array([ randrange(q) for _ in range(n) ])

    l = v.dot(s)
    
    V.append(v)
    L.append(l)
  
  return V,L

def experiment( A,b,q,hints,modular,fileName,verbose ):
  
  lattice = LWELattice(A,b,q,verbose=verbose)
  V, L = hints
  ctrHints = len(V)
  
  try:
    start = time.time()
    
    for i in range(ctrHints):
      if modular:
        lattice.integrateModularHint( V[i], L[i] % q, q )
      else:
        lattice.integratePerfectHint( V[i], L[i] )
    
    lattice.reduce(maxBlocksize=40)
    stop = time.time()
    
    output = "Finished experiment.\tHints: %d\tBlocksize: %d\tTime: %fs" % (ctrHints, lattice.successBlocksize, (stop-start))
    print( "\033[94m" + output + "\033[0m" )
    print(lattice.s)
    with open(fileName, "a+") as f:
      print(output, file=f)
      
        
  except Exception:
    #Print exceptions manualy, because Pool may hide them.
    output = traceback.format_exc()
    print(output)
    with open(fileName, "a+") as f:
      print(output, file=f)

args = parseArguments()

scheme = args["scheme"]

hints_min, hints_max, hints_step = args["hints"].split(":")
hints_min = int(hints_min)
hints_max  = int(hints_max)
hints_step = int(hints_step)
range_hints = range( hints_min, hints_max + 1, hints_step )

hints_centered = args["hints_centered"]
modular = args["modular"]

trials = args["trials"]
fileName = args["file"]
verbose = args["verbose"]


instances = []

for i in range(trials):
  if scheme=="test":
    n = hints_max*2
    A,b,q,s,e = generateToyInstance( n = n, m = n, q = 3329, eta = 3 )
  else:
    A,b,q,s,e = generateLWEInstance( scheme )
  
  print(s)
  
  hints = generateHints(s, q, hints_max + hints_step, hints_centered)
  
  for j in range_hints:
    hints = tuple( [  y for y in x[:-hints_step] ] for x in hints )
    instances.append( (A, b, q, hints, modular, fileName,verbose) )

pool = Pool()
pool.starmap( experiment, instances )