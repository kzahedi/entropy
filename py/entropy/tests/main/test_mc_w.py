"""
    >>> from entropy import container
    >>> from entropy import cmi

Create a container

    >>> x = container.Container(1000, 1)
    >>> y = container.Container(1000, 1)
    >>> z = container.Container(1000, 1)
    >>> for i in range(0, 1000):
    ...    x.add(math.cos(i/10.0))
    ...    y.add(math.cos(i/5.0))
    ...    z.add(math.cos(i/5.0) * math.sin(i/5.0))
    >>> x.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
    >>> y.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
    >>> z.setDomains([ (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0) ])
    >>> x.setBins([ 100, 100, 100 ])
    >>> y.setBins([ 100, 100, 100 ])
    >>> z.setBins([ 100, 100, 100 ])
    >>> x.discretise()
    >>> y.discretise()
    >>> z.discretise()
    >>> cmi.CMI(x,y,z)
    2.4802395350017234  


"""

__author__ = 'Keyan Ghazi-Zahedi, keyan.zahedi@gmail.com'

from testsuites import runModuleTestSuite

if __name__ == "__main__":
    runModuleTestSuite(__import__('__main__'))

