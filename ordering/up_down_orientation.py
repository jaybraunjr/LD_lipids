import MDAnalysis as mda
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

u = mda.Universe('test.gro','test.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')
c = u.select_atoms('resname CHYO and name O31')

#account for pbc
def wrap_z(z,zbox):
    if z < 0:
        z += zbox
    elif z > zbox:
        z -= zbox
    return z

#calculate molecules aiming down
def calc_director_down():
    N = len(a)
    sum = [0,0]
    for j in range(N):
        mag = vector_mag(a[j].position,b[j].position)
        unit_v = (a[j].position-b[j].position)/mag
        zbox = ts._unitcell[2]
        z_coord = wrap_z(c[j].position[2],zbox)
        if z_coord < 102:

            if unit_v[2] < 0:
                sum[0] += 1
            else:
                sum[1] += 1        
        else:
            pass
    return sum

#calculate molecules aiming up
def calc_director_up():
    N = len(a)
    sum1 = [0,0]
    for j in range(N):
        mag = vector_mag(a[j].position,b[j].position)
        unit_v = (a[j].position-b[j].position)/mag
        zbox = ts._unitcell[2]
        z_coord = wrap_z(c[j].position[2],zbox)
        if z_coord > 102:

            if unit_v[2] > 0:
                sum1[0] += 1
            else:
                sum1[1] +=1
        else:
            pass
    
    return sum1

### Vector magnitude
def vector_mag(a,b):
    v = a-b
    v = v*v
    sum = np.sum(v)
    return math.sqrt(sum)


directions=[]
directions1=[]
for ts in u.trajectory:
    direction_down = calc_director_down()
    direction_up = calc_director_up()
    directions.append([direction_down])

    directions1.append([direction_up])
arr=np.array(directions)
arr_reshaped = arr.reshape(arr.shape[0], -1)
np.savetxt("array.txt",arr_reshaped)

arr1=np.array(directions1)
arr_reshaped1 = arr1.reshape(arr1.shape[0], -1)
np.savetxt("array1.txt",arr_reshaped1)

