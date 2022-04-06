import MDAnalysis as mda
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


u = mda.Universe('nojump.gro','nojump.xtc')
a = u.select_atoms('resname CHYO and name C3')
b = u.select_atoms('resname CHYO and name C13')


print('Frames in trajecory:',len(u.trajectory))

# Magnitude of vector
def magnitude(a,b):
	v = a-b
	v = v*v
	sum = np.sum(v)
	return math.sqrt(sum)

# Solves the angle of the vectors, along with nematic order parameter
def solve(a,b):

	N = len(a)
	#print(len(a))
	sum = 0

	for j in range(N):
		if (a[j].resid != b[j].resid):
			assert(0)
		mag = magnitude(a[j].position,b[j].position)
		c_theta = (a[j].position[2]-b[j].position[2])/mag
		#print(c_theta)
		sum += (3*(c_theta*c_theta)-1)/2

	sum /= N
	return sum
print('Order parameter value:', solve(a,b))	
	
step = u.trajectory.time
#print(step)

# Calls on solve func, and appends the values into arrays
def test(a,b):
	
	mylist =[]

	for ts in u.trajectory:
		# calling function to solve order param in each timestep
		y = solve(a,b)
		# getting timestep
		x = u.trajectory.time
		mylist.append([x, y])
		#print(x)
		#print(y)
	return np.array(mylist)



# Plotting
p = test(a,b)
p.shape
p_plot = pd.DataFrame(p, columns=['Frame',
                                  'Smectic phase (cyclic)'])
p_plot.head()

# If you want txt file:

#p_file = open('order.txt', 'w+')
#p_file.write(p_plot.to_string())
#p_file.close()

# Plot with matlibplot

p_plot.plot(x='Frame')
plt.title('Orientational Order Parameter (S)')
plt.ylabel('Order') and plt.xlabel('ps')
plt.savefig('1.21_cyclic.png')

plt.draw()

print('***Success***')