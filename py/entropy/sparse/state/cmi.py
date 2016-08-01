import numpy
import math

def CMI(x,y,z):

    pxyz = numpy.zeros((x.data.max()+1, y.data.max()+1, z.data.max()+1)) # p(x,y,z)
    for index in range(0, x.data.shape[0]):
        xv = x.data[index,0]
        yv = y.data[index,0]
        zv = z.data[index,0]
        pxyz[xv, yv, zv] = pxyz[xv, yv, zv] + 1.0

    pxyz = pxyz / float(x.data.shape[0])

    pz  = numpy.sum(pxyz, axis=(0,1))
    pxz = numpy.sum(pxyz, axis=1)
    pyz = numpy.sum(pxyz, axis=0)

    r = 0.0

    r = numpy.zeros(x.data.shape[0])
    for index in range(0, x.data.shape[0]):
        xv = x.data[index,0]
        yv = y.data[index,0]
        zv = z.data[index,0]
        if (pxyz[xv,yv,zv] > 0.0 and pz[zv] > 0.0 and pxz[xv,zv] > 0.0 and pyz[yv,zv]):
                pxy_c_z = pxyz[xv,yv,zv] / pz[zv]
                px_c_z = pxz[xv,zv] / pz[zv] 
                py_c_z = pyz[yv,zv] / pz[zv] 
                r[index] = math.log(pxy_c_z,2.0) - math.log(px_c_z * py_c_z, 2.0)

    return r
