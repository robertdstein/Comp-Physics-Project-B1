import numpy as np
import math
import matplotlib.pyplot as plt
import csv
import minimise as m
import time
import readin as r
import plotting as p

#~ p.plot1d()

sv = [0.2, 0.5, 1.0]
tolerance = 10**-10

#~ xsol, limits = m.run(sv, tolerance)
#~ meansigma = 0.5*abs(limits[0]-limits[1])
#~ print "Mean sigma", meansigma

p.plotNsigma(sv, tolerance)
