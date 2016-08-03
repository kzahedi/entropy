#!/usr/bin/env python

import sys
sys.path.append("../../py/")

from entropy import container, mc_w 
from entropy.state import mc_w as mc_w_d

fd = open("data/dcmot.csv","r")
data = fd.readlines()
fd.close()
dcmot = []
for line in data:
    values = [float(v) for v in line.split(",")]
    dcmot.append([values[1], values[2], values[3], values[9]])

fd = open("data/muslin.csv","r")
data = fd.readlines()
fd.close()
muslin = []
for line in data:
    values = [float(v) for v in line.split(",")]
    muslin.append([values[1], values[2], values[3], values[9], values[4]])

fd = open("data/musfib.csv","r")
data = fd.readlines()
fd.close()
musfib = []
for line in data:
    values = [float(v) for v in line.split(",")]
    musfib.append([values[1], values[2], values[3], values[9], values[4]])

dcmotW = container.Container(len(dcmot), 3)
dcmotS = container.Container(len(dcmot), 2)
dcmotA = container.Container(len(dcmot), 1)

muslinW = container.Container(len(muslin), 3)
muslinS = container.Container(len(muslin), 1)
muslinA = container.Container(len(muslin), 1)

musfibW = container.Container(len(musfib), 3)
musfibS = container.Container(len(musfib), 1)
musfibA = container.Container(len(musfib), 1)

for i in range(0, len(dcmot)):
    dcmotW.add(dcmot[i][0])
    dcmotW.add(dcmot[i][1])
    dcmotW.add(dcmot[i][2])
    dcmotS.add(dcmot[i][0])
    dcmotS.add(dcmot[i][1])
    dcmotA.add(dcmot[i][3])

for i in range(0, len(muslin)):
    muslinW.add(muslin[i][0])
    muslinW.add(muslin[i][1])
    muslinW.add(muslin[i][2])
    muslinS.add(muslin[i][4])
    muslinA.add(muslin[i][3])

for i in range(0, len(musfib)):
    musfibW.add(musfib[i][0])
    musfibW.add(musfib[i][1])
    musfibW.add(musfib[i][2])
    musfibS.add(musfib[i][4])
    musfibA.add(musfib[i][3])

p_min = min(dcmotW.data[:,0].min(), muslinW.data[:,0].min(),
        musfibW.data[:,0].min())
p_max = max(dcmotW.data[:,0].max(), muslinW.data[:,0].max(),
        musfibW.data[:,0].max())

v_min = min(dcmotW.data[:,1].min(), muslinW.data[:,1].min(),
        musfibW.data[:,1].min())
v_max = max(dcmotW.data[:,1].max(), muslinW.data[:,1].max(),
        musfibW.data[:,1].max())

a_min = min(dcmotW.data[:,2].min(), muslinW.data[:,2].min(),
        musfibW.data[:,2].min())
a_max = max(dcmotW.data[:,2].max(), muslinW.data[:,2].max(),
        musfibW.data[:,2].max())

mi_min = min(muslinS.data[:,0].min(), musfibS.data[:,0].min())
mi_max = max(muslinS.data[:,0].max(), musfibS.data[:,0].max())

ac_min = dcmotA.data[:,0].min()
ac_max = dcmotA.data[:,0].max()

print "Domains:"
print "  Position:          " + str(p_min)  + " " + str(p_max)
print "  Velocity:          " + str(v_min)  + " " + str(v_max)
print "  Acceleration:      " + str(a_min)  + " " + str(a_max)
print "  Action (DC Motor): " + str(ac_min) + " " + str(ac_max)


dcmotW.normaliseColumn(0,  p_min,  p_max)
dcmotW.normaliseColumn(1,  v_min,  v_max)
dcmotW.normaliseColumn(2,  a_min,  a_max)
dcmotS.normaliseColumn(0,  p_min,  p_max)
dcmotS.normaliseColumn(1,  v_min,  v_max)
dcmotA.normaliseColumn(0,  ac_min, ac_max)

dcmotW.setBins([300, 300, 300])
dcmotW.setDomains([ (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)])
dcmotS.setBins([300, 300])
dcmotS.setDomains([ (0.0, 1.0), (0.0, 1.0)])
dcmotA.setBins([300])
dcmotA.setDomains([ (0.0, 1.0) ])


muslinW.normaliseColumn(0, p_min,  p_max)
muslinW.normaliseColumn(1, v_min,  v_max)
muslinW.normaliseColumn(2, a_min,  a_max)
muslinS.normaliseColumn(0, mi_min, mi_max)

muslinW.setBins([300, 300, 300])
muslinW.setDomains([ (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)])
muslinS.setBins([300, 300])
muslinS.setDomains([ (0.0, 1.0), (0.0, 1.0)])
muslinA.setBins([300])
muslinA.setDomains([ (0.0, 1.0) ])


musfibW.setBins([300, 300, 300])
musfibW.setDomains([ (0.0, 1.0), (0.0, 1.0), (0.0, 1.0)])
musfibS.setBins([300, 300])
musfibS.setDomains([ (0.0, 1.0), (0.0, 1.0)])
musfibA.setBins([300])
musfibA.setDomains([ (0.0, 1.0) ])
musfibW.normaliseColumn(0, p_min,  p_max)
musfibW.normaliseColumn(1, v_min,  v_max)
musfibW.normaliseColumn(2, a_min,  a_max)
musfibS.normaliseColumn(0, mi_min, mi_max)


dcmotW.discretise()
dcmotS.discretise()
dcmotA.discretise()

muslinW.discretise()
muslinS.discretise()
muslinA.discretise()

musfibW.discretise()
musfibS.discretise()
musfibA.discretise()

dcmotW2 = dcmotW.drop(1)
dcmotW1 = dcmotW.drop(-1)
dcmotA1 = dcmotA.drop(-1)

muslinW2 = muslinW.drop(1)
muslinW1 = muslinW.drop(-1)
muslinA1 = muslinA.drop(-1)

musfibW2 = musfibW.drop(1)
musfibW1 = musfibW.drop(-1)
musfibA1 = musfibA.drop(-1)

mc_w_dcmot  = mc_w.MC_W(dcmotW2,  dcmotW1,  dcmotA1)
mc_w_muslin = mc_w.MC_W(muslinW2, muslinW1, muslinA1)
mc_w_musfib = mc_w.MC_W(musfibW2, musfibW1, musfibA1)

print "MC_W: "
print " dcmot:  " + str(mc_w_dcmot)
print " muslin: " + str(mc_w_muslin)
print " musfib: " + str(mc_w_musfib)

mc_w_dcmot_state  = mc_w_d.MC_W(dcmotW2,  dcmotW1,  dcmotA1)
mc_w_muslin_state = mc_w_d.MC_W(muslinW2, muslinW1, muslinA1)
mc_w_musfib_state = mc_w_d.MC_W(musfibW2, musfibW1, musfibA1)

print "MC_W: "
print " dcmot:  " + str(mc_w_dcmot_state)
print " muslin: " + str(mc_w_muslin_state)
print " musfib: " + str(mc_w_musfib_state)

fd = open("dcmot.csv","w")
for i in range(0, mc_w_dcmot_state.shape[0]):
    fd.write(str(mc_w_dcmot_state[i]) + "\n")
fd.close()

fd = open("muslin.csv","w")
for i in range(0, mc_w_muslin_state.shape[0]):
    fd.write(str(mc_w_muslin_state[i]) + "\n")
fd.close()

fd = open("musfib.csv","w")
for i in range(0, mc_w_musfib_state.shape[0]):
    fd.write(str(mc_w_musfib_state[i]) + "\n")
fd.close()

