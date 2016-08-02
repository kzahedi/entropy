#!/usr/bin/env python


# /analysis erzeugen
# /plots    erzeugen
# boolean flags zum an-/abschalten von funktionen

import math
import os
import shutil
import argparse
from decimal import *
from convertfunctions import *

parser = argparse.ArgumentParser(description="Arguments:")
parser.add_argument("-d", type=str, default=None, help="Parent directory.")
args = parser.parse_args()

c_states = "control.states.csv"
s_states = "hand.sofastates.txt"
s_csv    = "hand.sofastates.csv"

# subdirectories = [x[0] for x in os.walk(args.d)]
subdirectories = os.listdir(args.d)
subdirectories = filter(lambda x: os.path.isdir(args.d + "/" + x) == True, subdirectories)

for s in subdirectories:
    try:
        print "creating " + args.d + s + "/analysis"
        os.mkdir(args.d + "/" + s + "/analysis")
    except:
        shutil.rmtree(args.d + "/" + s + "/analysis")
        os.mkdir(args.d + "/" + s + "/analysis")

subdirectories = [args.d + "/"  + v             for v in subdirectories]
controlstates  = [v + "/raw/" + c_states        for v in subdirectories]
sofastates     = [v + "/raw/" + s_states        for v in subdirectories]
sofacsvstates  = [v + "/analysis/" + s_csv      for v in subdirectories]
analysisdirs   = [v + "/analysis/"              for v in subdirectories]
wdomains       = [v + "/analysis/W.domains.csv" for v in subdirectories]
adomains       = [v + "/analysis/A.domains.csv" for v in subdirectories]

# convert sofa state files
for i in sofastates:
    if os.path.exists(i):
        print "working on " + i.split("/")[-3]
        fd = open(i, "r")
        lines = fd.readlines()
        fd.close()

        data = [get_position(line) for line in lines if "T=" not in line]

        o = i.replace(".txt",".csv")
        o = o.replace("raw","analysis")
        print "writing " + "/".join(o.split("/")[-3:])
        fd = open(o,"w")
        for d in data:
            fd.write(convert_to_csv(d) + "\n")
        fd.close()

# creating sofa domain file
sofa_min_values = None
sofa_max_values = None
for c in sofacsvstates:

    print "reading " + c 
    fd = open(c,"r")
    lines = fd.readlines()
    fd.close()

    for line in lines:
        values = [Decimal(v) for v in line.split(",")]
        if sofa_min_values is None:
            sofa_min_values = [v for v in values]
            sofa_max_values = [v for v in values]
        else:
            for i in range(0, len(values)):
                if values[i] > sofa_max_values[i]:
                    sofa_max_values[i] = values[i]
                if values[i] < sofa_min_values[i]:
                    sofa_min_values[i] = values[i]

print str(sofa_max_values[4])

for w in wdomains:
    fd = open(w,"w")
    fd.write(",".join([str(v) for v in sofa_min_values]) + "\n")
    fd.write(",".join([str(v) for v in sofa_max_values]) + "\n")
    fd.close()

# creating control domain file
control_min_values = None
control_max_values = None
for c in controlstates:
    fd = open(c,"r")
    lines = fd.readlines()[1:]
    fd.close()

    for line in lines:
        values = [Decimal(v) for v in line.split(",")]
        if control_min_values is None:
            control_min_values = [v for v in values]
            control_max_values = [v for v in values]
        else:
            for i in range(0, len(values)):
                if values[i] > control_max_values[i]:
                    control_max_values[i] = values[i]
                if values[i] < control_min_values[i]:
                    control_min_values[i] = values[i]


for a in adomains:
    fd = open(a,"w")
    fd.write(",".join([str(v) for v in control_min_values]) + "\n")
    fd.write(",".join([str(v) for v in control_max_values]) + "\n")
    fd.close()


# copy control states
for c in controlstates:
    shutil.copyfile(c, c.replace("raw","analysis"))


