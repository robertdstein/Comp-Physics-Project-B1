import math
import numpy as np
import time
import readin as r

data = r.run()

def parabolic(startingvals, tolerance, array=data):
	start = time.time()
	[x0,x1, x2] = startingvals

	y0 = nll([x0], array)
	y1 = nll([x1], array)
	y2 = nll([x2], array)
	
	val_list = [[x0, y0], [x1,y1], [x2, y2]]

	oldx3 = (x0-tolerance)
	x3 = x2 + tolerance

	while math.fabs(x3 - oldx3) > tolerance:
		[[x0, y0], [x1,y1], [x2, y2]] = val_list
		oldx3 = x3
		x3 = 0.5 * ((x2**2 - x1**2)*y0 + (x0**2 - x2**2)*y1 + (x1**2 - x0**2)*y2)/((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2)
		y3 = nll([x3], array)
		val_list.append([x3, y3])
		val_list.sort(key=lambda x: x[1])
		val_list.pop()
		val_list.sort(key=lambda x: x[0])
	
	end = time.time()
	print "Solution found! x =", x3, "leading to y =", y3, "Elapsed time", end-start
	return x3
	
def secant(f, a, b, tolerance):
	while abs(0.5*(b-a)) > tolerance:
		tmp = b - f(b)*(b-a)/(f(b)-f(a))
		a=b
		b=tmp
	return b
	
def bisect(f, a, b, tolerance):
	c=0.5*(a+b)
	i=0
	while abs(0.5*(b-a)) > tolerance:
		if f(c) == 0:
			return c
		elif f(a) * f(c) < 0:
			b=c
		else:
			a=c
		c = 0.5*(a+b)
		i+=1
	return c
	
def nll(v,array=data):
	if len(v)> 1:
		tau = v.item(0)
		a = v.item(1)
	else:
		tau = v[0]
		a = 1.0
		
	ll = 0
	for entry in array:
		t = entry[0]
		sigma = entry[1]
		erfcval = math.erfc(math.sqrt(0.5)*((sigma/tau) - (t/sigma)))
		l = a*math.exp(0.5*(sigma/tau)**2 - (t/tau))*erfcval/(2*tau)
		l += (1.0-a)*math.exp(-0.5*(t/sigma)**2)/(sigma*math.sqrt(2*math.pi))
		ll -= math.log(l)
	return ll

def limitsll(sol, tolerance, array=data):
	truell = nll([sol], array)
	intervals = [-1, 1]
	lims = []
	
	def f(x):
		return nll([x], array) - (truell + 0.5)	
	
	for val in intervals:
		scale = (1+(0.01*val))
		
		start = time.time()
		
		while f(sol*scale) < 0:
			scale = scale**2
			
		second = scale*sol
		points = [sol, second]
		points.sort()
		
		[a,b] = points
		x = secant(f, a, b, tolerance)
		
		if (x - sol)*val < 0:
			print "Limit", x, "found does not lie on correct side of solution", sol
			print "Repeating with simple bisection method"
			x = bisect(f, a, b, tolerance)
	
		lims.append(x)
		end = time.time()
		print "Limit", x, "Value", nll([x]), "Sigma =", math.fabs(x -sol), "Elapsed time", end-start
	return lims
	
def run(startingvals, tolerance, array=data):
	sol = parabolic(startingvals, tolerance, array)
	lims = limitsll(sol, tolerance, array)
	meansigma = 0.5*abs(lims[0]-lims[1])
	return sol, lims, meansigma
	
def quasinewton(v, tolerance, lims=None):
	nvariables = len(v)
	G0 =  np.identity(nvariables)
	
	start = time.time()
	
	alpha = 0.00001
	
	v0 = v
	dv0 = getdv(v)
	i = 0
	
	tr = 2 * tolerance
	
	while tr > tolerance:
		newv = v0 - (alpha * G0.dot(dv0))
		#~ print i, newv
		
		if type(lims) != type(None):
			for j in range(len(newv)):
				val = newv.item(j)
				llim = lims.item(2*j)
				ulim = lims.item((2*j)+1)
				#~ print val, llim, ulim
				if (val > ulim) or (val < llim):
					newv[j] = v[j]
					print "Minimnisation exceeded limits, as", val, "does not lie between", llim, ulim
					print "Replacing alpha =", alpha, "with"
					alpha *= 0.1
					print alpha, "and beginning again with", newv

		delta = newv - v0
		gamma = (getdv(newv) - dv0)
		
		v0 = newv
		dv0 = getdv(newv)
		op=np.outer(delta, delta)
		
		G = G0 + (op/np.inner(gamma, delta).trace()) - G0.dot(op).dot(G0)/np.inner(gamma,G0.dot(gamma)).trace()
		G0 = G
		
		
		tr = math.sqrt(np.inner(delta, delta).trace())
		
		i+=1

	end = time.time()
		
	print "Solution found! tau =", v0.item(0), "and a =", v0.item(1), "leading to y =", nll(v0), "Elapsed time", end-start
		
def getdv(v, f=nll):
	delta = 0.01
	derivs = []
	oldy = f(v)
	
	for val in v:
		val += delta
		new = f(v)
		val -= delta
		diff = (new-oldy)/delta
		derivs.append([diff])
		
	derivs = np.matrix(derivs)
	return derivs
	
def outerprod(a, b):
	if len(a) == len(b):
		ab = []
		for i in range(len(a)):
			row = []
			for j in range(len(b)):
				ai = a.item(i)
				bj = b.item(j)
				abij = ai*bj
				row.append(abij)
			ab.append(row)
		ab = np.matrix(ab)
		return ab	
	else:
		raise Exception("Vectors are not the same length")	

