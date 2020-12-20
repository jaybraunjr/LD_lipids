import MDAnalysis as mda
import numpy as np

u  = mda.Universe('90.10_1us.tpr', '90.10_1us.xtc')
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
      'TRIO': u.select_atoms('resname TRIO')}

### ag['TRIO'].names             <-- array(['C118', 'H18A', 'H18B', ..., 'H18X', 'H18Y', 'H18Z'], dtype=object) 
### ag['TRIO'].names.astype(str) <-- array(['C118', 'H18A', 'H18B', ..., 'H18X', 'H18Y', 'H18Z'], dtype='<U4')
names = (ag['TRIO'].names).astype(str)

times = []; Nsurfs = []; Rsurfs = []; n = [];
for ts in u.trajectory:
    times.append(ts.time/1000) #ps to ns
    utz = ag['UT'].center_of_mass()[2] # the COM of the PL upper tail
    ltz = ag['LT'].center_of_mass()[2] # the COM of the PL lower tail
    
    ### First select all the triolein atoms that are above utz or below ltz
    ### Then select only oxygen atoms among them
    
    bA1 = ((ag['TRIO'].positions[:,2] > utz) | \
          (ag['TRIO'].positions[:,2] < ltz)) & \
          np.char.startswith(names,'O')
    
    ### Get the residue ids
    
    res = ag['TRIO'].resids[bA1]

   
    ### There are 6 oxygen atoms per a TRIO molecule.
    ### r[0] <-- resids; r[1] <-- counts
    ### select a resid (r[0]) in that all of its oxygen atoms are above utz or below ltz (r[1] == 6)
    
    r = np.unique(res, return_counts = True)

    surf_tg_resids = r[0][r[1]==6]
    
    #surf_tg_resids = r[0] ## MILD CONDITION: if at least one oxygen atom is above utz or below ltz.

    #print out the surf_tg_resids
   
    if surf_tg_resids.size > 8:
        n.append(surf_tg_resids)
        #print(n)
        np.savetxt("surf_tg_resids_90.10.txt",n,fmt = "%s")
  

    
    ### The number of SURF-TG (including both upper and lower leaflets)
    Nsurf = len(surf_tg_resids)
    Nsurfs.append(Nsurf)
    
    ### The ratio of SURF-TG respective to PL
    Rsurf = Nsurf / (Nsurf + N_PL) * 100
    Rsurfs.append(Rsurf)

save = np.transpose([times, Nsurfs, Rsurfs])
np.savetxt('tg_90.10.dat', save, fmt='%.3f %i %.3f', header='time (ns)   SURF-TG count   SURF-TG ratio')