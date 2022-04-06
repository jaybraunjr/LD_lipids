import MDAnalysis as mda
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


u = mda.Universe('whole.gro','whole.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')
c = u.select_atoms('resname CHYO and name O31')



def calc_director():
	N = len(a)
	sum = [0,0]

	for j in range(N):
		mag = vector_mag(a[j].position,b[j].position)
		unit_v = (a[j].position-b[j].position)/mag
		z_coord = c[j].position[2]
		if unit_v[2] < 0:
			cc=unit_v*-1
			if z_coord < 40:
				sum[0] += cc
			else:
				sum[1] += cc
		else:
			if z_coord < 40:
				sum[0] += unit_v
			else:
				sum[1] += unit_v
	sum[0] /= vector_mag(sum[0],[0,0,0])
	sum[1] /= vector_mag(sum[1],[0,0,0])
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
	sum = [0,0]
	nmol = [0,0]
	for j in range(N):
		z_coord = c[j].position[2]
		c1 = (a[j].position-b[j].position)
		if z_coord < 40:
			dot = np.dot(c1,director[0])
			cos_theta=dot/vector_mag(a[j].position,b[j].position)
			sum[0] += (3*(cos_theta*cos_theta)-1)/2
			nmol[0] += 1
		else:
			dot = np.dot(c1,director[1])
			cos_theta=dot/vector_mag(a[j].position,b[j].position)
			sum[1] += (3*(cos_theta*cos_theta)-1)/2
			nmol[1] +=1
	sum[0] /= nmol[0]
	sum[1] /= nmol[1]
	return sum


orderparams = [[],[]]
for ts in u.trajectory:
	director = calc_director()
	orderparam = calc_orderparam(director)
	orderparams[0].append(orderparam[0])
	orderparams[1].append(orderparam[1])

plt.plot(orderparams[0],'.')
plt.plot(orderparams[1],'.')
plt.ylabel('order')

plt.xlabel('frame')
plt.legend(['0','1'])
plt.savefig('order.png')




