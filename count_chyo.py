import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

u  = mda.Universe('50.50_3us.part0010.gro', 'short.xtc')
halfz = u.dimensions[2]/2

C2 = ' '.join(['C2%d' %i for i in range(2, 22)])
C3 = ' '.join(['C3%d' %i for i in range(2, 22)])

### CHANGE POPC DOPE SAPI OF COURSE IF YOU HAVE DIFFERENT PL COMPOSITION!!!
### Select the PL tail groups of the upper leaflet
N_PL = u.select_atoms('resname POPC DOPE SAPI').n_residues

us = '(same residue as resname POPC DOPE SAPI \
        and name P and prop z > %f) \
        and name ' %halfz + C2 + C3

### Select the PL tail groups of the lower leaflet
ls = '(same residue as resname POPC DOPE SAPI \
        and name P and prop z < %f) \
        and name ' %halfz + C2 + C3

### Atomic groups organized into the dictionary
ag = {'PL': u.select_atoms('resname POPC DOPE SAPI'),
      'UT': u.select_atoms(us),
      'LT': u.select_atoms(ls),
      'CHYO': u.select_atoms('resname CHYO')}

### ag['TRIO'].names             <-- array(['C118', 'H18A', 'H18B', ..., 'H18X', 'H18Y', 'H18Z'], dtype=object) 
### ag['TRIO'].names.astype(str) <-- array(['C118', 'H18A', 'H18B', ..., 'H18X', 'H18Y', 'H18Z'], dtype='<U4')
names = (ag['CHYO'].names).astype(str)

times = []; Nsurfs = []; Rsurfs = []; n = [];
for ts in u.trajectory:
    times.append(ts.time/1000) #ps to ns
    utz = ag['UT'].center_of_mass()[2] # the COM of the PL upper tail
    ltz = ag['LT'].center_of_mass()[2] # the COM of the PL lower tail
    
    ### First select all the triolein atoms that are above utz or below ltz
    ### Then select only oxygen atoms among them
    
    bA1 = ((ag['CHYO'].positions[:,2] > utz) | \
          (ag['CHYO'].positions[:,2] < ltz)) & \
          np.char.startswith(names,'O')
          #np.char.startswith(names,'C318' and 'C317' and 'C316' and 'C315' and 'C314' and 'C313' and 'C312' and 'C311' and 'C310' and 'C39' and 'C38' and 'C37' and 'C36')
    
    ### Get the residue ids
    
    res = ag['CHYO'].resids[bA1]

   
    ### There are 6 oxygen atoms per a TRIO molecule.
    ### r[0] <-- resids; r[1] <-- counts
    ### select a resid (r[0]) in that all of its oxygen atoms are above utz or below ltz (r[1] == 6)
    
    r = np.unique(res, return_counts = True)

    #surf_chyo_resids = r[0][r[1]==2]
    
    surf_chyo_resids = r[0] ## MILD CONDITION: if at least one oxygen atom is above utz or below ltz.

    #print out the surf_tg_resids
   
    if surf_chyo_resids.size > 3:
        n.append(surf_chyo_resids)
        #print(n)
        np.savetxt("surf_chyo_O.txt",n,fmt = "%s")
  

    
    ### The number of SURF-TG (including both upper and lower leaflets)
    Nsurf = len(surf_chyo_resids)
    Nsurfs.append(Nsurf)
    
    ### The ratio of SURF-TG respective to PL
    Rsurf = Nsurf / (Nsurf + N_PL) * 100
    Rsurfs.append(Rsurf)

save = np.transpose([times, Nsurfs])
np.savetxt('surch_chyo_O.dat', save, fmt='%.3f %i', header='time (ns)   SURF-CHYO count   SURF-CHYO ratio')

plt.plot(Nsurfs)


plt.title('SURF-CHYO (O)', y=0.9)
plt.ylabel('amount SURF-CHYO')
plt.ylim(0,10)

plt.xlabel('frame')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig('SURF_CHYO_O.png')
