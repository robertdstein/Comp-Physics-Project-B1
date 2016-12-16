import numpy as np

def DerivOrder2():
	vals = ["f", "df/dx", "d2f/dx2"]
	f = [1,0,0]
	fp = [1,1,0.5]
	fm = [1,-1,0.5]
	fpp = [1,2,2]
	
	steps = ["f(x)","f(x+h)", "f(x-h)"]
	
	m = np.matrix([f, fp, fm])
	print m
	i = m.getI()
	print i
	
	for j in range(len(vals)):
		print vals[j], "=",
		for k in range(len(steps)):
			print "(", i.item(j, k), steps[k], ") +",
		print "..."

#~ DerivOrder2()

def BivarDerivOrder2():
	vals = ["f", "df/dx", "df/dy", "d2f/dx2", "d2f/dxdy","d2f/dy2"]
	f = [1,0,0,0,0,0]
	fxp = [1,1,0,0.5,0,0]
	fxm = [1,-1,0,0.5,0,0]
	fyp = [1,0,1,0,0,0.5]
	fym = [1,0,-1,0,0,0.5]
	fxpyp = [1,1,1,0.5,1,0.5]
	
	steps = ["f(x, y)","f(x+h, y)", "f(x-h, y)", "f(x, y+g)", "f(x, y-g)", "f(x+h, y+g)"]
	
	m = np.matrix([f, fxp, fxm, fyp, fym, fxpyp])
	print m
	i = m.getI()
	print i
	
	for j in range(len(vals)):
		print vals[j], "=",
		for k in range(len(steps)):
			print "(", i.item(j, k), steps[k], ") +",
		print "..."
		
#~ BivarDerivOrder2()
	
