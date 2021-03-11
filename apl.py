import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = []
with open('output.csv') as csvfile:
	#print(csvfile)
	read = csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)
	for row in read:
		data.append(row)
	
		#print(frame)

x = np.array(data)
df = pd.DataFrame(x)

df[1] = (df[1]**2)/135
df[0]/1000
df.to_csv('90.10.csv', index=False)
print(df)

#x = df[0]/1000
x = df[0]
y = df[1]


plt.plot(x,y, markersize=1)
plt.ylim(0.55,0.75)
plt.title('area per lipid (90:10)')
plt.xlabel('time (ns)')
plt.ylabel('nm$^2$')
plt.grid(axis='y')
plt.savefig('apl.png')
print('Finished')

