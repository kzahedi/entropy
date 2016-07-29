import numpy
import math
from entropy import container
from entropy import cmi


x = container.Container(1000, 1)
y = container.Container(1000, 1)
z = container.Container(1000, 1)
for i in range(0, 1000):
    x.add(math.cos(i/10.0))
    y.add(math.cos(i/5.0))
    z.add(math.cos(i/5.0) * math.sin(i/5.0))
x.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
y.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
z.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
x.setBinsSizes([ 100, 100, 100 ])
y.setBinsSizes([ 100, 100, 100 ])
z.setBinsSizes([ 100, 100, 100 ])
x.discretise()
y.discretise()
z.discretise()
print cmi.CMI(x,y,z)

