"""
    >>> import numpy
    >>> from entropy import container

Create a container

    >>> c = container.Container(3, 5)
    >>> input = c.shape()
    >>> output = (3,5)

Filling data by row

    >>> for i in range(1,16): c.add(i)
    >>> c.data
    array([[  1.,   2.,   3.,   4.,   5.],
           [  6.,   7.,   8.,   9.,  10.],
           [ 11.,  12.,  13.,  14.,  15.]])

Filling data by column

    >>> c.fillIndex = 0
    >>> c.setAddByColumn()
    >>> for i in range(1,16): c.add(i)
    >>> c.data
    array([[ 1.,  4.,  7.,  10.,  13.],
           [ 2.,  5.,  8.,  11.,  14.],
           [ 3.,  6.,  9.,  12.,  15.]])

Testing discretisation by column

    >>> c.setBinsSizes( [10, 20, 30, 40, 50] )
    >>> c.setDomains( [ (1,3), (4,6), (7,9), (10,12), (13,15)] )
    >>> c.uniformDiscretisationByColumn()
    >>> c.data
    array([[ 0.,  0.,  0., 0.,  0.],
           [ 5., 10., 15., 20., 25.],
           [ 9., 19., 29., 39., 49.]])

Relabel the data

    >>> c.relabel()
    >>> c.data
    array([[ 0., 0., 0., 0., 0.],
           [ 1., 1., 1., 1., 1.],
           [ 2., 2., 2., 2., 2.]])

Combine into a single column data matrix

    >>> c.combineDiscretisedColumns()
    >>> c.data
    array([[      0.],
           [ 246211.],
           [ 492422.]])

"""

__author__ = 'Keyan Ghazi-Zahedi, keyan.zahedi@gmail.com'

from testsuites import runModuleTestSuite

if __name__ == "__main__":
    runModuleTestSuite(__import__('__main__'))

