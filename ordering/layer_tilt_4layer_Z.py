import MDAnalysis as mda
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


u = mda.Universe('whole.gro','whole.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')
c = u.select_atoms('resname CHYO and name C10')

def wrap_z(z,zbox):
	if z < 0:
		z += zbox
	elif z > zbox:
		z -= zbox
	return z

def calc_director():
	N = len(a)
	sum = [0,0,0,0]

	for j in range(N):
		mag = vector_mag(a[j].position,b[j].position)
		unit_v = (a[j].position-b[j].position)/mag
		zbox = ts._unitcell[2]
		z_coord = wrap_z(c[j].position[2],zbox)
		if unit_v[2] < 0:
			cc=unit_v*-1
			if z_coord < 40:
				sum[0] += cc
			elif 40 < z_coord < 80:
				 sum[1] += cc
			elif 80 < z_coord < 120:
				 sum[2] += cc
			else:
				sum[3] += cc
		else:
			if z_coord < 40:
				sum[0] += unit_v
			elif 40 < z_coord < 80:
				 sum[1] += unit_v
			elif 80 < z_coord < 120:
				 sum[2] += unit_v
			else:
				sum[3] += unit_v
	sum[0] /= vector_mag(sum[0],[0,0,0])
	sum[1] /= vector_mag(sum[1],[0,0,0])
	sum[2] /= vector_mag(sum[2],[0,0,0])
	sum[3] /= vector_mag(sum[3],[0,0,0])
	return sum


### Vector magnitude

def vector_mag(a,b):
	v = a-b
	v = v*v
	sum = np.sum(v)
	return math.sqrt(sum)


### Solve for the order parameter that includes angle of unit vector and molecular axis

def calc_orderparam(a,b):
	N = len(a)
	sum = [0,0,0,0]
	nmol = [0,0,0,0]
	for j in range(N):
		zbox = ts._unitcell[2]
		z_coord = wrap_z(c[j].position[2],zbox)
		c1 = (a[j].position-b[j].position)
		if z_coord < 40:

			mag = vector_mag(a[j].position,b[j].position)
			cos_theta=(a[j].position[2]-b[j].position[2])/mag
			
			sum[0] += (3*(cos_theta*cos_theta)-1)/2
			nmol[0] += 1
		elif 40 < z_coord < 80:
			mag = vector_mag(a[j].position,b[j].position)
			cos_theta=(a[j].position[2]-b[j].position[2])/mag
			sum[1] += (3*(cos_theta*cos_theta)-1)/2
			nmol[1] += 1
		elif 80 < z_coord < 120:
			mag = vector_mag(a[j].position,b[j].position)
			cos_theta=(a[j].position[2]-b[j].position[2])/mag
			sum[2] += (3*(cos_theta*cos_theta)-1)/2
			nmol[2] += 1
		else:
			mag = vector_mag(a[j].position,b[j].position)
			cos_theta=(a[j].position[2]-b[j].position[2])/mag
			sum[3] += (3*(cos_theta*cos_theta)-1)/2
			nmol[3] +=1
	sum[0] /= nmol[0]
	sum[1] /= nmol[1]
	sum[2] /= nmol[2]
	sum[3] /= nmol[3]
	return sum


orderparams = [[],[],[],[]]
for ts in u.trajectory:
	#director = calc_director()
	orderparam = calc_orderparam(a,b)
	orderparams[0].append(orderparam[0])
	orderparams[1].append(orderparam[1])
	orderparams[2].append(orderparam[2])
	orderparams[3].append(orderparam[3])

plt.plot(orderparams[0],'.', label='1')
plt.plot(orderparams[1],'.', label='2')
plt.plot(orderparams[2],'.', label='3')
plt.plot(orderparams[3],'.', label='4')
plt.legend()
plt.ylabel('order')
plt.xlabel('frame')
plt.savefig('testing.png')




