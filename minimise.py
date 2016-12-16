import math
import numpy as np
import time
from numpy.linalg import inv
import readin as r

data = r.run()

def nll(v,array=data):
	if not isinstance(v, np.matrix):
		np.matrix(v)
	if len(v)> 1:
		tau = v.item(0)
		a = v.item(1)
	else:
		tau = v.item(0)
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

def parabolic(startingvals, tolerance, f=nll, i=0):
	start = time.time()
	[a,b, c] = startingvals
	
	x0 = a.copy()
	x1 = b.copy()
	x2 = c.copy()

	y0 = f(x0)
	y1 = f(x1)
	y2 = f(x2)
	
	val_list = [[x0, y0], [x1,y1], [x2, y2]]

	oldx3 = x0.copy()
	oldx3.itemset(i, x0.item(i)-tolerance)
	x3 = x2.copy()
	x3.itemset(i, x2.item(i) + tolerance)

	while math.fabs(f(x3) - f(oldx3)) > tolerance:
		[[x0, y0], [x1,y1], [x2, y2]] = val_list
		oldx3 = x3.copy()
		newx3val = (0.5 * ((x2.item(i)**2 - x1.item(i)**2)*y0 + (x0.item(i)**2 - x2.item(i)**2)*y1 + (x1.item(i)**2 - x0.item(i)**2)*y2)/
					((x2.item(i)-x1.item(i))*y0 + (x0.item(i)-x2.item(i))*y1 + (x1.item(i)-x0.item(i))*y2))
		x3.itemset(i, newx3val)
		y3 = f(x3)
		val_list.append([x3.copy(), y3])
		val_list.sort(key=lambda x: x[1])
		val_list.pop()
		val_list.sort(key=lambda x: x[0].item(i))
	
	end = time.time()
	print "Solution found! x =", x3, "leading to y =", y3, "Elapsed time", end-start
	return x3
	
def secant(f, a, b, tolerance, i=0):
	va = a.copy()
	vb = b.copy()
	while abs(0.5*(vb.item(i)-va.item(i))) > tolerance:
		tmp = vb.item(i) - f(vb)*(vb.item(i)-va.item(i))/(f(vb)-f(va))
		va.itemset(i, vb.item(i))
		vb.itemset(i, tmp)
	return vb
	
def bisect(f, a, b, tolerance, i=0):
	va = a.copy()
	vb = b.copy()
	vc = a.copy()
	vc.itemset(i, 0.5*(va.item(i)+vb.item(i)))
	while abs(0.5*(va.item(i)-vb.item(i))) > tolerance:
		if f(vc) == 0:
			return vc
		elif f(va) * f(vc) < 0:
			vb=vc
		else:
			va=vc
		vc.itemset(i, 0.5*(va.item(i)+vb.item(i)))
	return vc

def limitsll(sol, tolerance, f):
	truell = f(sol)
	intervals = [-1, 1]
	alllims=[]
	
	def g(x):
		return f(x) - (truell + 0.5)	
	
	for i in range(len(sol)):
		lims = []
		trueval = sol.item(i)
		for val in intervals:
			copy = sol.copy()
			
			scale = (1+(0.01*val))
			var = copy.item(i)*scale
			
			copy.itemset(i, var)
			
			start = time.time()
			
			while g(copy) < 0:
				scale = scale**2
				var *= scale
				copy.itemset(i, var)
				
			points = [sol, copy]
			points.sort(key=lambda x: x[i])
			[a,b] = points
			
			x = secant(g, a, b, tolerance, i=i)
								
			if (x.item(i) - trueval)*val < 0:
				print "Limit", x, "found does not lie on correct side of solution", trueval
				print "Repeating with simple bisection method"
				x = bisect(g, a, b, tolerance)
		
			lims.append(x.item(i))
			end = time.time()
			print "Limit", x.item(i), "Value", f(x), "Sigma =", math.fabs(x.item(i) -trueval), "Elapsed time", end-start
		alllims.append(lims)
	return alllims
	
def run1d(startingvals, tolerance, f=nll):
	sol = parabolic(startingvals, tolerance, f)
	lims = limitsll(sol, tolerance, f)
	meansigma = 0.5*abs(lims[0][0]-lims[0][1])
	return sol, lims, meansigma
	
def run2d(v, tolerance, lims=None, f=nll):
	sol, path = quasinewton(v, tolerance, lims, f)
	print path
	lims = limitsll(sol, tolerance, f)
	return lims, path
	
def quasinewton(v, tolerance, lims=None,f=nll):
	nvariables = len(v)
	G0 =  np.identity(nvariables)
	
	start = time.time()
	
	alpha = 0.00001
	
	v0 = v
	dv0 = getdv(v, f)
	i = 0
	
	tr = 2 * tolerance
	
	path=[[],[]]
	
	while not (tolerance > tr > 0):
		
		if tr < 0:
			G0 *= -1
		newv = v0 - (alpha * G0.dot(dv0))
		print i, newv, G0
		
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
		gamma = (getdv(newv, f) - dv0)
		
		tr = f(v0) - f(newv)
		
		print tr
		
		v0 = newv
		dv0 = getdv(newv, f)
		op=np.outer(delta, delta)
		
		path[0].append(v0.item(0))
		path[1].append(v0.item(1))
		
		G = G0 + (op/np.inner(gamma, delta).trace()) - G0.dot(op).dot(G0)/np.inner(gamma,G0.dot(gamma)).trace()
		G0 = G
		
		i+=1

	end = time.time()
		
	print "Solution found in", i, "steps! variable_1 =", v0.item(0), "and variable_2 =", v0.item(1) 
	print "leading to y =", f(v0), "with tolerance", tolerance, "Elapsed time", end-start
	print "Covariance Matrix is", 
	
	covar = inv(hessian2d(v0, f))
	print covar,
	xsig = math.sqrt(covar.item(0))
	ysig = math.sqrt(covar.item(3))
	print "leaving sigmas", xsig, ysig
	
	return np.matrix([[path[0][-1]], [path[1][-1]]]), path
	
def getdv(v, f):
	delta = 0.001
	derivs = []
	oldy = f(v)
	
	for val in v:
		val += delta
		up = f(v)
		val -= 2*delta
		low = f(v)
		val += delta
		diff = 0.5*(up-low)/delta
		derivs.append([diff])
		
	derivs = np.matrix(derivs)
	return derivs

def hessian2d(v, f):
	delta = 0.001
	derivs = []
	func = f(v)
	v.itemset(0, v.item(0)-delta)
	fxm = f(v)
	v.itemset(0, v.item(0) + 2*delta)
	fxp = f(v)
	v.itemset(1, v.item(1) + delta)
	fxpyp  = f(v)
	v.itemset(0, v.item(0) - delta)
	fyp = f(v)
	v.itemset(1, v.item(1) - 2* delta)
	fym = f(v)
	v.itemset(1, v.item(1) + delta)
	
	print "values", func, fxp, fyp, fxm, fym, fxpyp
	
	d2fdx2 = (fxp + fxm - 2*func)/(delta**2)
	d2fdy2 = (fyp + fym - 2*func)/(delta**2)
	d2fdydx = (func + fxpyp - fyp - fxp)/(delta**2)
	
	hessian = np.matrix([[d2fdx2, d2fdydx], [d2fdydx, d2fdy2]])
	
	return hessian
	
#~ def getd2v(v, f=nll):
	#~ delta = 0.001
	#~ derivs = []
	#~ oldy = f(v)
	
	#~ for val1 in v:
		#~ derivrow = []
		#~ val1 += delta
		#~ new = f(v)
		#~ d1 = (new-oldy)/delta
		#~ for val2 in v:
			#~ val2 += delta
			#~ new = f(v)
			#~ val2 -= delta
			#~ d2 = (new-oldy)/delta
			#~ diff = (d2 - d1)/delta
			#~ derivrow.append(diff)
		#~ val1 -= delta
		#~ derivs.append(derivrow)
		
	#~ derivs = np.matrix(derivs)
	#~ print derivs
	#~ return derivs
	
#~ def outerprod(a, b):
	#~ if len(a) == len(b):
		#~ ab = []
		#~ for i in range(len(a)):
			#~ row = []
			#~ for j in range(len(b)):
				#~ ai = a.item(i)
				#~ bj = b.item(j)
				#~ abij = ai*bj
				#~ row.append(abij)
			#~ ab.append(row)
		#~ ab = np.matrix(ab)
		#~ return ab	
	#~ else:
		#~ raise Exception("Vectors are not the same length")	

