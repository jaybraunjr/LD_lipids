import MDAnalysis as mda
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


u = mda.Universe('whole.gro','whole.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')
c = u.select_atoms('resname CHYO and name O31')


def wrap_z(z,zbox):
	if z < 0:
		z += zbox
	elif z > zbox:
		z -= zbox
	return z

def calc_director():
	N = len(a)
	sum = [0,0]

	for j in range(N):
		mag = vector_mag(a[j].position,b[j].position)
		unit_v = (a[j].position-b[j].position)/mag
		zbox = ts._unitcell[2]
		z_coord = wrap_z(c[j].position[2],zbox)
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
		zbox = ts._unitcell[2]
		z_coord = wrap_z(c[j].position[2],zbox)
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

combined = (np.add(orderparams[0],orderparams[1]))/2



plt.plot(orderparams[0],'.', label='bottom layer')
plt.plot(orderparams[1],'.', label='top layer')
plt.plot(combined, label='average')
plt.title('Orientational Order Parameter (per layer)', y=0.9)
plt.ylabel('order(S)')
plt.ylim(0.728,0.778)
plt.xlabel('frame')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig('order.png')




