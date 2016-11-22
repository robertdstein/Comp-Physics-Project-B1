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
	plt.savefig("graphs/decaytimes.pdf")
	plt.close()
	
	taurange = np.linspace(0.1, 3.0, 100)
	
	llvals =[]
	print taurange
	for tau in taurange:
		llvector = m.nll([tau], data)
		llvals.append(llvector)
		
	plt.figure()
	plt.plot(taurange, llvals)
	plt.xlabel(r"$\tau$")
	plt.ylabel("-Log Likelihood")
	plt.title(r"$\tau$ Optimisation")
	plt.savefig("graphs/nll.pdf")
	plt.close()
	
def plot2d():
	times = data[:,0]
	
	nbins=50
	fig = plt.figure()
		
	x1 = np.linspace(0.1, 3.0, nbins)
	x2 = np.linspace(0.2, 1.0, nbins)
	
	xs = [x1,x2]
	
	y1 = np.linspace(0.01, 1.0, nbins)
	y2 = np.linspace(0.5, 1.0, nbins)
	
	ys=[y1,y2]
	
	for i in range(2):
		
		plt.subplot(2,1,i+1)
		
		x=xs[i]
		y=ys[i]
		
		grid=[]
		for a in y:
			llvals =[]
			for tau in x:
				llvector = m.nll(tau, array=data, a=a)
				llvals.append(llvector)
			grid.insert(0,llvals)
		
		cmap = cm.jet_r
		grid = np.array(grid)
		
		print grid
		
		
		plt.imshow(grid, aspect="auto", extent=(x[0], x[-1], y[0], y[-1]), interpolation='bilinear', cmap=cmap)	
		plt.xlabel(r"$\tau$")
		plt.ylabel("a")
		plt.colorbar()
	
	plt.suptitle(r"$\tau$ and a Optimisation")
	plt.savefig("graphs/nll2d.pdf")
	plt.close()

def plotNsigma(startingvals, tolerance):
	
	lognrange = np.linspace(1.0, math.log10(len(data)), 10)
	print lognrange
	
	logsigmas=[]
	
	for logn in lognrange:
		n = int(10**logn)
		sol, lims, sigma = m.run(startingvals, tolerance, data[:n])
		logsigmas.append(math.log10(sigma))
	
	y0 = logsigmas[0]
	y1 = logsigmas[-1]
	x0 = lognrange[0]
	x1 = lognrange[-1]
	
	grad = (y1-y0)/(x1-x0)
	c = y0 - (grad*x0)
	
	def f(x):
		return grad*x + c
		
	def g(y):
		x= (y-c)/grad
		n = 10**x
		return x,n
		
	x,n = g(-15)
	
	line = r"y = " + str.format('{0:.2f}',grad) + "x + " + str.format('{0:.2f}', c) + " \n"
	line += "For $Log(\sigma)=$-15, we have $Log(N)=$" + str.format('{0:.2f}',x) + " \n"
	line += "Then $N=$"+str.format('{:.3G}',int(n))
	
	plt.figure()
	plt.annotate(line, xy=(0.1, 0.1), xycoords="axes fraction")
	plt.scatter(lognrange, logsigmas)
	plt.plot(lognrange, f(lognrange), label="Fit")
	plt.legend()
	plt.xlabel(r"$Log(N)$")
	plt.ylabel(r"$Log(\sigma)$" )
	plt.title(r"$Log(\sigma)$ against $Log(N)$")
	plt.savefig("graphs/nsigma.pdf")
	plt.close()
	
