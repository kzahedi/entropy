import scipy
import numpy

def Matrix3D(class):

    data = None
    xdim = -1
    ydim = -1
    zdim = -1

    def __init__(self, x, y, z):

        self.data = scipy.dok_matrix((x*y, z), dtype=numpy.float32)
        self.xdim = x
        self.ydim = y
        self.zdim = z

    def __getitem__(self, x, y, z):
        return self.data[x + y * self.xdim, z]

    def __setitem__(self, x, y, z, value):
        self.data[x + y * self.xdim, z] = value
