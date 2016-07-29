import numpy

class Container:
    data          = None
    fillIndex     = 0
    addByRow      = True
    rows          = -1
    cols          = -1
    isDiscretised = False

    def __init__(self, rows, cols):
        self.fillIndex     = 0
        self.data          = numpy.zeros((rows, cols))
        self.rows          = rows
        self.cols          = cols
        self.addByRow      = True
        self.isDiscretised = False

    def setAddByRow(self):
        self.addByRow = True

    def setAddByColumn(self):
        self.addByRow = False
        
    def add(self, valum):
        c = -1
        r = -1
        if self.addByRow is True:
            c = self.fillIndex % self.cols
            r = self.fillIndex / self.cols
        else:
            c = self.fillIndex / self.rows
            r = self.fillIndex % self.rows

        self.data[r,c] = valum
        self.fillIndex = self.fillIndex + 1

    def shape(self):
        return self.data.shape

    def setBinsSizes(self, lst):
        self.bins = lst

    def setDomains(self, lst_of_tuples):
        self.domains = lst_of_tuples

    def discretise(self):
        self.discretiseByColumn()
        self.relabel() # count from 0,1,..., i.e. rename the bins
        self.combineDiscretisedColumns()
        self.relabel() # count from 0,1,..., i.e. rename the bins

    def discretiseByColumn(self):
        # currently only uniform binning is supported
        self.uniformDiscretisationByColumn()

    def uniformDiscretisationByColumn(self):
        for c in range(0, self.cols):
            for r in range(0, self.rows):
                assert (self.domains[c][0] <= self.data[r,c] and self.data[r,c] <= self.domains[c][1])
                mapped  = int(((self.data[r,c] - self.domains[c][0]) / (self.domains[c][1] - self.domains[c][0])) * self.bins[c])
                cropped = int(min(self.bins[c]-1, mapped))
                self.data[r,c] = cropped
        self.isDiscretised = True

    def relabel(self):
        for c in range(0, self.cols):
            values = []
            for r in range(0, self.rows):
                try:
                    index = values.index(self.data[r,c])
                except ValueError:
                    index = len(values)
                    values.append(self.data[r,c])
                self.data[r,c] = index


    def combineDiscretisedColumns(self):
        new_data = numpy.zeros((self.rows,1))
        for r in range(0, self.rows):
            v = 0
            f = 1
            for c in range(0, self.cols):
                if c > 0:
                    f  = f * self.bins[c-1]
                v += f * self.data[r,c]
            new_data[r,0] = v
        self.data = new_data
        self.cols = 1
