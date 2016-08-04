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
parser.add_argument("-d",  type=str, default=None, help="Parent directory.")
parser.add_argument("-wi", type=str, default=None, help="W indices.")
parser.add_argument("-ai", type=str, default=None, help="A indices.")
parser.add_argument("-wb", type=int, default=300, help="W bins.")
parser.add_argument("-ab", type=int, default=300, help="A bins.")
parser.add_argument("-start", type=int, default=0, help="start time index")
parser.add_argument("-wf", type=bool, default=True, help="Convert to wrist frame.")
parser.add_argument("-csv", type=bool, default=True, help="Export data.")
parser.add_argument("-redoall", type=bool, default=False, help="Redo all steps.")
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

for s in subdirectories:
    try:
        print "creating " + args.d + s + "/analysis"
        os.mkdir(args.d + "/" + s + "/analysis")
    except:
        print args.d + s + "/analysis already existed"
        # shutil.rmtree(args.d + "/" + s + "/analysis")
        # os.mkdir(args.d + "/" + s + "/analysis")

subdirectories = [os.path.normpath(args.d + "/"  + v)             for v in subdirectories]
controlstates  = [os.path.normpath(v + "/raw/" + c_states)        for v in subdirectories]
sofastates     = [os.path.normpath(v + "/raw/" + s_states)        for v in subdirectories]
sofacsvstates  = [os.path.normpath(v + "/analysis/" + s_csv)      for v in subdirectories]
analysisdirs   = [os.path.normpath(v + "/analysis/")              for v in subdirectories]
wdomains       = [os.path.normpath(v + "/analysis/W.domains.csv") for v in subdirectories]
adomains       = [os.path.normpath(v + "/analysis/A.domains.csv") for v in subdirectories]

# convert sofa state files
converted_a_file = False
for i in sofastates:
    o = i.replace(".txt",".csv")
    o = o.replace("raw","analysis")
    if os.path.exists(i) == True and (os.path.exists(o) == False or args.redoall == True):
        converted_a_file = True
        print "working on " + i.split("/")[-3]
        fd = open(i, "r")
        lines = fd.readlines()
        fd.close()

        data = None

        if args.wf:
            data = [get_position_local_wrist(line) for line in lines if "T=" not in line]
        else:
            data = [get_position(line) for line in lines if "T=" not in line]

        print "writing " + "/".join(o.split("/")[-3:])
        fd = open(o,"w")
        for d in data[0+args.start:]:
            fd.write(convert_to_csv(d) + "\n")
        fd.close()

# creating sofa domain file
if converted_a_file is True:
    sofa_min_values = None
    sofa_max_values = None
    for c in sofacsvstates:

        print "reading " + "/".join(c.split("/")[-3:])

        if os.path.exists(c):
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

    for w in wdomains:
        fd = open(w,"w")
        fd.write(",".join([str(v) for v in sofa_min_values]) + "\n")
        fd.write(",".join([str(v) for v in sofa_max_values]) + "\n")
        fd.close()

    # creating control domain file
    control_min_values = None
    control_max_values = None
    for c in controlstates:

        if os.path.exists(c):
            fd = open(c,"r")
            lines = fd.readlines()[1+args.start:]
            fd.close()
            lines = lines[1+args.start:len(lines)-1]

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

    # copy control states - remove first line
    for c in controlstates:
        fd = open(c,"r")
        lines = fd.readlines()
        fd.close()
        fd = open(c.replace("raw","analysis"),"w")
        for line in lines[1+args.start:len(lines)-1]:
            fd.write(line)
        fd.close()
        # shutil.copyfile(c, c.replace("raw","analysis"))

for a in analysisdirs:
    print "./rbo_mc " + option_string + " -d " + "/".join(a.split("/")[-3:])
    os.system("./rbo_mc " + option_string + " -d " + a)
