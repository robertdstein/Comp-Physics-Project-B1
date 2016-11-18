import math
import time
import readin as r

data = r.run()

def parabolic(startingvals, tolerance, array=data):
	start = time.time()
	[x0,x1, x2] = startingvals

	y0 = nll(x0, array)
	y1 = nll(x1, array)
	y2 = nll(x2, array)
	
	val_list = [[x0, y0], [x1,y1], [x2, y2]]

	oldx3 = (x0-tolerance)
	x3 = x2 + tolerance

	while math.fabs(x3 - oldx3) > tolerance:
		[[x0, y0], [x1,y1], [x2, y2]] = val_list
		oldx3 = x3
		x3 = 0.5 * ((x2**2 - x1**2)*y0 + (x0**2 - x2**2)*y1 + (x1**2 - x0**2)*y2)/((x2-x1)*y0 + (x0-x2)*y1 + (x1-x0)*y2)
		y3 = nll(x3, array)
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
	
def nll(tau,array=data):
	ll = 0
	for entry in array:
		t = entry[0]
		sigma = entry[1]
		erfcval = math.erfc(math.sqrt(0.5)*((sigma/tau) - (t/sigma)))
		ll -= math.log(0.5/tau) + 0.5*(sigma/tau)**2 - (t/tau) + math.log(erfcval)
	return ll

def limitsll(sol, tolerance, array=data):
	truell = nll(sol, array)
	intervals = [-1, 1]
	lims = []
	
	def f(x):
		return nll(x, array) - (truell + 0.5)	
	
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
		print "Limit", x, "Value", nll(x), "Sigma =", math.fabs(x -sol), "Elapsed time", end-start
	return lims
	
def run(startingvals, tolerance, array=data):
	sol = parabolic(startingvals, tolerance, array)
	lims = limitsll(sol, tolerance, array)
	meansigma = 0.5*abs(lims[0]-lims[1])
	return sol, lims, meansigma
	

