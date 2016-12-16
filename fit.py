import numpy as np
import math
import matplotlib.pyplot as plt
import csv
import minimise as m
import time
import readin as r
import plotting as p

#~ p.plot1d()

tolerance = 10**-4
sv1d = [np.matrix([0.2]), np.matrix([0.5]), np.matrix([1.0])]


xsol, limits, meansigma = m.run1d(sv1d, tolerance)
print "Mean sigma", meansigma

#~ p.plotNsigma(sv1d, tolerance)

sv2d = np.matrix([[0.5], [0.99]])
lims = np.matrix([[0.0, 3.0], [0.0, 1.0]])
print sv2d, lims
path, limits=m.run2d(sv2d, tolerance, lims)
#~ print limits

#~ print path

#~ p.plot2d(path)

#~ def f(v):
	#~ return math.exp((1.5*v.item(0)) + (2*v.item(1)))

#~ v = np.matrix([[0.],[0.]])
#~ m.getd2v(v, f)
