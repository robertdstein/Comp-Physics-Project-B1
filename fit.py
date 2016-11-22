import numpy as np
import math
import matplotlib.pyplot as plt
import csv
import minimise as m
import time
import readin as r
import plotting as p

#~ p.plot1d()
#~ p.plot2d()

sv = [0.2, 0.5, 1.0]
tolerance = 10**-10

xsol, limits, meansigma = m.run(sv, tolerance)
print "Mean sigma", meansigma

#~ p.plotNsigma(sv, tolerance)


sv = np.matrix([[1.0], [0.5]])
lims = np.matrix([[0.0, 3.0], [0.0, 1.0]])
print sv, lims
m.quasinewton(sv, tolerance, lims)
