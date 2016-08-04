#!/usr/bin/env python

# /analysis erzeugen
# /plots    erzeugen
# boolean flags zum an-/abschalten von funktionen

import math
import os
import re
import shutil
import argparse
from decimal import *
from convertfunctions import *

parser = argparse.ArgumentParser(description="Arguments:")
parser.add_argument("-d",  type=str, default=None, help="Parent directory.")
parser.add_argument("-wi", type=str, default=None, help="W indices.")
parser.add_argument("-ai", type=str, default=None, help="A indices.")
parser.add_argument("-wb", type=int, default=300, help="W bins.")
parser.add_argument("-ab", type=int, default=300, help="A bins.")
parser.add_argument("-csv", type=bool, default=True, help="Export data.")
parser.add_argument("-analyseonly", type=int, default=1, help="only run analyse scripts.")
parser.add_argument("-rbo", type=str,
        default="/Users/zahedi/projects/builds/entropy-build/bin/rbo_mc",
        help="RBO Binary.")
args = parser.parse_args()

option_string = "--wbins " + str(args.wb) + " --abins " + str(args.ab)
if args.wi is not None:
    option_string = option_string + " --wi " + args.wi
if args.ai is not None:
    option_string = option_string + " --ai " + args.ai
if args.csv is True:
    option_string = option_string + " -csv"

c_states = "control.states.csv"
s_states = "hand.sofastates.txt"
s_csv    = "hand.sofastates.csv"

# subdirectories = [x[0] for x in os.walk(args.d)]
subdirectories = os.listdir(args.d)
subdirectories = filter(lambda x: os.path.isdir(args.d + "/" + x) == True and os.path.isdir(args.d + "/" + x + "/raw") == True, subdirectories)

parents = [re.search('^(.+)-object', v).group(1) for v in subdirectories]
uniqueparents = list(set(parents))

if args.analyseonly is 0:
    try:
        os.mkdir(args.d + "/clustered/")
    except:
        shutil.rmtree(args.d + "/clustered/")
        os.mkdir(args.d + "/clustered/")

    for p in parents:
        try:
            os.mkdir(args.d + "/clustered/" + p)
        except:
            pass
        try:
            os.mkdir(args.d + "/clustered/" + p + "/analysis")
        except:
            pass

    for i in range(0, len(parents)):
        if os.path.exists(args.d + "/" + "/clustered/" + p + "/analysis/hand.sofastates.csv"):
            os.remove(args.d + "/" + "/clustered/" + p + "/analysis/hand.sofastates.csv")
        if os.path.exists(args.d + "/" + "/clustered/" + p + "/analysis/control.states.csv"):
            os.remove(args.d + "/" + "/clustered/" + p + "/analysis/control.states.csv")

    for i in range(0, len(parents)):
        d = subdirectories[i]
        p = parents[i]

        print "reading " + d

        fd = open(args.d + "/" + d + "/analysis/hand.sofastates.csv","r")
        lines = fd.read()
        fd.close()

        fd = open(args.d + "/" + "/clustered/" + p + "/analysis/" + "/hand.sofastates.csv","a")
        fd.write(lines)
        fd.close()

        fd = open(args.d + "/" + d + "/analysis/control.states.csv","r")
        lines = fd.read()
        fd.close()

        fd = open(args.d + "/" + "/clustered/" + p + "/analysis/" + "/control.states.csv","a")
        fd.write(lines)
        fd.close()

        
    sofa_min_values = None
    sofa_max_values = None
    for p in uniqueparents:
        c = args.d + "/clustered/" + p + "/analysis/" + "/hand.sofastates.csv"
        if os.path.exists(c):
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

    for p in uniqueparents:
        w = args.d + "/clustered/" + p + "/analysis/" + "/W.domains.csv"
        print "writing " + w
        fd = open(w,"w")
        fd.write(",".join([str(v) for v in sofa_min_values]) + "\n")
        fd.write(",".join([str(v) for v in sofa_max_values]) + "\n")
        fd.close()

    control_min_values = None
    control_max_values = None
    for p in uniqueparents:
        c = args.d + "/clustered/" + p + "/analysis/" + "/control.states.csv"
        if os.path.exists(c):
            fd = open(c,"r")
            lines = fd.readlines()
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

    for p in uniqueparents:
        a = args.d + "/clustered/" + p + "/analysis/" + "/A.domains.csv"
        fd = open(a,"w")
        fd.write(",".join([str(v) for v in control_min_values]) + "\n")
        fd.write(",".join([str(v) for v in control_max_values]) + "\n")
        fd.close()



for p in uniqueparents:
    print "./rbo_mc " + option_string + " -d " + args.d + "/clustered/" + p + "/analysis"
    os.system("./rbo_mc " + option_string + " -d " + args.d + "/clustered/" + p + "/analysis")



