import MDAnalysis as mda
import numpy as np
import math
import pandas as pd


u = mda.Universe('test.gro','test.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')

def calc_director():
	N = len(a)
	sum = 0

	for j in range(N):
		mag = vector_mag(a[j].position,b[j].position)
		unit_v = (a[j].position-b[j].position)/mag
		
		if unit_v[2] < 0:
			c=unit_v*-1
			sum += c
		else:
			sum += unit_v
	sum /= vector_mag(sum,[0,0,0])
	return sum


### Vector magnitude

def vector_mag(a,b):
	v = a-b
	v = v*v
	sum = np.sum(v)
	return math.sqrt(sum)


### Solve for the order parameter that includes angle of unit vector and molecular axis

def calc_orderparam(director):
	N = len(a)
	sum = 0

	for j in range(N):
		c1 = (a[j].position-b[j].position)
		dot = np.dot(c1,director)
		cos_theta=dot/vector_mag(a[j].position,b[j].position)
		sum += (3*(cos_theta*cos_theta)-1)/2
	sum /=N
	return sum


orderparams = []
for ts in u.trajectory:
	director = calc_director()
	orderparam = calc_orderparam(director)
	orderparams.append(orderparam)
	print(orderparam)





