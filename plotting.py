import numpy as np
import math
import matplotlib.pyplot as plt
import readin as r
import minimise as m
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

data = r.run()

def plot1d():
	times = data[:,0]
	
	plt.figure()
	plt.hist(times, bins=50)
	plt.xlabel("Reconstructed time t (ps)")
	plt.ylabel("Count")
	plt.savefig("graphs/decaytimes.pdf")
	plt.close()
	
	taurange = np.linspace(0.1, 3.0, 100)
	
	llvals =[]
	print taurange
	for tau in taurange:
		llvector = m.nll(np.matrix([tau]), data)
		llvals.append(llvector)
		
	plt.figure()
	plt.plot(taurange, llvals)
	plt.xlabel(r"$\tau (10^{-12} s)$")
	plt.ylabel("-Log Likelihood")
	plt.title(r"$\tau$ Optimisation")
	plt.savefig("graphs/nll.pdf")
	plt.close()
	
def plot2d(path=None):
	times = data[:,0]
	
	nbins=10
	fig = plt.figure()
		
	x1 = np.linspace(0.1, 3.0, nbins)
	x2 = np.linspace(0.40, 0.42, nbins)
	
	xs = [x1,x2]
	
	y1 = np.linspace(0.01, 1.0, nbins)
	y2 = np.linspace(0.97, 1.0, nbins)
	
	ys=[y1,y2]
	
	for i in range(2):
		
		plt.subplot(2,1,i+1)
		
		x=xs[i]
		y=ys[i]
		
		grid=[]
		for a in y:
			llvals =[]
			for tau in x:
				v = np.matrix([[tau],[a]])
				llvector = m.nll(v, array=data)
				llvals.append(llvector)
			grid.insert(0,llvals)
		
		cmap = cm.jet_r
		grid = np.array(grid)
		
		plt.imshow(grid, aspect="auto", extent=(x[0], x[-1], y[0], y[-1]), interpolation='bilinear', cmap=cmap)
		plt.colorbar()
		if type(path) != type(None):
			plt.plot(path[0], path[1], color="white")
			plt.scatter(path[0][0], path[1][0], color="white", marker="x")
			plt.scatter(path[0][-1], path[1][-1], color="white", marker="*")
			
		plt.ylim(y[0], y[-1])
		plt.xlim(x[0], x[-1])
		plt.xlabel(r"$\tau$")
		plt.ylabel("a")
	
	plt.suptitle(r"$\tau$ and a Optimisation")
	plt.savefig("graphs/nll2d.pdf")
	plt.close()

def plotNsigma(startingvals, tolerance):
	
	lognrange = np.linspace(1.0, math.log10(len(data)), 10)
	print lognrange
	
	logsigmas=[]
	yuperr = []
	ydownerr = [] 
	
	for logn in lognrange:
		n = int(10**logn)
		
		def f(x):
			return m.nll(x, array=data[:n])
		
		sol, lims, sigma = m.run1d(startingvals, tolerance, f=f)
		logsigmas.append(math.log10(sigma)-12)
		#~ yuperr.append(math.log(logn+math.sqrt(10**logn)))
		#~ ydownerr.append(math.log(logn-math.sqrt(10**logn)))
	
	y0 = logsigmas[0]
	y1 = logsigmas[-1]
	x0 = lognrange[0]
	x1 = lognrange[-1]
	
	guessgrad = (y1-y0)/(x1-x0)
	guessc = y0 - (guessgrad*x0)

		
	def residuals(v):
		res = 0
		
		def line(x):
			return v.item(0)*x + v.item(1)
		#~ print v
		
		for i in range(len(lognrange)):
			logn = lognrange[i]
			expected = line(logn)
			true = logsigmas[i]
			diff = expected-true
			res += (diff**2)/(10**(logn))**2
			#~ print logn, expected, true, diff, res
		
		return res
	
	sv= np.matrix([[guessgrad], [guessc]])
	sol, path = m.quasinewton(sv, tolerance=10**-6, lims=None, f=residuals)
	
	print sol
	
	grad = sol.item(0)
	c = sol.item(1)
	
	def h(x):
		return 10**(grad*x + c)
		
	def g(y):
		x= (y-c)/grad
		n = 10**x
		return x,n
		
	x,n = g(-15)
	
	yerr = []
	sigma = []
	
	for i in range(len(logsigmas)):
		sig = logsigmas[i]
		true = 10**sig
		logn = lognrange[i]
		err = true/math.sqrt(10**logn)
		yerr.append(err)
		sigma.append(true)
	
	print yerr
	print sigma
	
	line = r"y = " + str.format('{0:.1f}',grad) + "x + " + str.format('{0:.1f}', c) + " \n"
	line += "For $Log(\sigma)=$-15, we have $Log(N)=$" + str.format('{0:.1f}',x) + " \n"
	line += "Then $N=$"+str.format('{:.2G}',int(n))
	
	plt.figure()
	plt.yscale("log")
	plt.annotate(line, xy=(0.1, 0.1), xycoords="axes fraction")
	plt.scatter(lognrange, sigma)
	plt.errorbar(lognrange, sigma, yerr=yerr, xerr=0.0, fmt='o', color="red")
	plt.plot(lognrange, h(lognrange), label="Fit")
	plt.legend()
	plt.xlabel(r"$Log(N)$")
	plt.ylabel(r"$Log(\sigma)$" )
	plt.title(r"$Log(\sigma)$ against $Log(N)$")
	plt.savefig("graphs/nsigma.pdf")
	plt.close()
	
