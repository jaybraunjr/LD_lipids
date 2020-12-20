import MDAnalysis as mda
import numpy as np

def density_frame(pos, mass, pbc, bins):
    dz = bins[1] - bins[0]
    h, _ = np.histogram(pos, weights=mass, bins=bins)
    h /= pbc[0] * pbc[1] * dz * 0.602214
    return h

def cal_overlap(d1, d2):
    #ov = 4 * d1 * d2 / (d1 + d2)**2
    #ov[np.isnan(ov)] = 0
    thr = 0.1
    d_sum = d1 + d2
    d_mul = d1 * d2
    d_sum[d1 + d2 < thr] = 1
    d_mul[d1 + d2 < thr] = 0
    
    ov = 4 * d_mul / d_sum**2
    return ov


def cal_inter(ov, dz):
    interdigitation = np.sum(ov) * dz
    return interdigitation/10 #to nm

def interdigit(u, nbins=50, nblocks=5, b=0, e=None):
    halfz = u.dimensions[2]/2
    numP  = u.select_atoms('name P').n_atoms
    
    C2 = ' '.join(['C2%d' %i for i in range(2, 22)])
    C3 = ' '.join(['C3%d' %i for i in range(2, 22)])
    us = '(same residue as resname POPC DOPE SAPI and name P and prop z>{}) and name '.format(halfz) + C2 + C3
    ls = '(same residue as resname POPC DOPE SAPI and name P and prop z<{}) and name '.format(halfz) + C2 + C3
    #print(us)
    #print(ls)

    groups = {'memb':  u.select_atoms('resname POPC DOPE SAPI'),
              'umemb': u.select_atoms(us),
              'lmemb': u.select_atoms(ls),
              'trio':  u.select_atoms('resname TRIO')}
    
    assert groups['umemb'].n_atoms == groups['lmemb'].n_atoms, "# of umemb atoms != # of lmemb atoms"

    names  = groups['trio'].names #object type
    names  = names.astype(str)    #str type
    resids = groups['trio'].resids
    
    print("Number of phospholipid atoms: ", groups['memb'].n_atoms)
    print("Number of upper tail carbon atoms: ", groups['umemb'].n_atoms)
    print("Number of lower tail carbon atoms: ", groups['lmemb'].n_atoms)
    print("Number of trio atoms: ", groups['trio'].n_atoms)
    
    
    times = []
    zs    = []
    total_inter  = []
    strong_inter = []
    weak_inter   = []
    
    d0_densities = [] #PL
    d1_densities = [] #TRIO
    d2_densities = [] #SURF-TRIO
    d3_densities = [] #CORE-TRIO

    total_ov = []
    strong_ov = []
    weak_ov = []
    
    strong_num   = []
    

    for ts in u.trajectory[b:e]:
        if int(ts.time/1000) % 1000 == 0:
            print("analyzing %d us.... " %(ts.time/1000/1000))
        
        pbc = u.dimensions
        bins = np.linspace(0, pbc[2], nbins+1)
        dz = bins[1] - bins[0]
        trio_pos = groups['trio'].positions
    
        utz = np.average(groups['umemb'].positions[:,2])
        ltz = np.average(groups['lmemb'].positions[:,2])
        # print("upper tail z (A): ", utz)
        # print("lower tail z (A): ", ltz)
        
        ### d0: density of all membrane atoms
        d0 = density_frame(groups['memb'].positions[:,2],
                           groups['memb'].masses,
                           pbc  = pbc,
                           bins = bins)
        d0_densities.append(d0)
        
        ### TOTAL INTERDIGITATION
        ### d1: density of all trio atoms
        d1 = density_frame(groups['trio'].positions[:,2],
                           groups['trio'].masses,
                           pbc  = pbc,
                           bins = bins)
        d1_densities.append(d1)
        
        tov = cal_overlap(d0, d1) #total overlap
        tin = cal_inter(tov, dz)  #total interdigitation
        total_ov.append(tov)
        total_inter.append(tin)
    

        ### STRONG INTERDIGITATION
        # Select corresponding oxygen atoms => get resid => select molecules corresponding resid
        # 'All' oxygen atoms of one TRIO molecule should be upper than utz or lower than ltz
        boolArray  = ((trio_pos[:,2] > utz) | (trio_pos[:,2] < ltz)) & np.char.startswith(names, 'O')
        strong_resids = resids[boolArray]
        r = np.unique(strong_resids, return_counts=True)
        boolArray = (r[1] == 6)
        strong_resids = r[0][boolArray]
        boolArray = np.isin(resids, strong_resids)
        pp = trio_pos[boolArray]
        mm = groups['trio'].masses[boolArray]
        
        ### d2: density of SURF-TRIO
        d2 = density_frame(pp[:,2], mm, pbc, bins)
        d2_densities.append(d2)

        sov = cal_overlap(d0, d2) #strong overlap
        sin = cal_inter(sov, dz)  #strong interdigitation
        strong_ov.append(sov)
        strong_inter.append(sin)
    
        ### WEAK INTERDIGITATION
        #boolArray  = ((trio_pos[:,2] <= utz) & (trio_pos[:,2] >= ltz)) & np.char.startswith(names, 'O')
        #weak_resids = resids[boolArray]
        #boolArray = np.isin(resids, weak_resids)

        boolArray = np.isin(resids, strong_resids, invert=True)
        pp = trio_pos[boolArray]
        mm = groups['trio'].masses[boolArray]
        
        ### d3: density of CORE-TRIO 
        d3 = density_frame(pp[:,2], mm, pbc, bins)
        d3_densities.append(d3)

        wov = cal_overlap(d0, d3)
        win = cal_inter(wov, dz)
        weak_ov.append(wov)
        weak_inter.append(win)
        
        # print("resids: ", len(set(resids)))
        # print("strong resids: ", len(set(strong_resids)))
        # print("weak resids: ", len(set(weak_resids)))
        strong_num.append(len(strong_resids))

        times.append(ts.time/1000)
        zs.append(pbc[2]/10)
    
    strong_num = np.array(strong_num)
    XX = np.linspace(0, np.average(zs), nbins)
    results = {}

    results['inter']  = {}
    results['inter']['total']  = np.transpose([times, total_inter])
    results['inter']['strong'] = np.transpose([times, strong_inter])
    results['inter']['weak']   = np.transpose([times, weak_inter])
    
    results['ov'] = {}
    results['ov']['total']     = np.transpose([XX, np.average(total_ov, axis=0)])
    results['ov']['strong']    = np.transpose([XX, np.average(strong_ov, axis=0)])
    results['ov']['weak']      = np.transpose([XX, np.average(weak_ov, axis=0)])
 
    results['ratio'] = {}
    results['ratio']['num']             = np.transpose([times, strong_num])
    results['ratio']['trio-to-pl']      = np.transpose([times, strong_num/numP])
    results['ratio']['trio-to-pl+trio'] = np.transpose([times, strong_num/(numP + strong_num)])
    
    results['density'] = {}
    results['density']['PL']        = np.transpose([XX, np.average(d0_densities, axis=0)])
    results['density']['TRIO']      = np.transpose([XX, np.average(d1_densities, axis=0)])
    results['density']['SURF-TRIO'] = np.transpose([XX, np.average(d2_densities, axis=0)])
    results['density']['CORE-TRIO'] = np.transpose([XX, np.average(d3_densities, axis=0)])

    print("units: Z (nm), interdigitation (nm), time (ns), density (g/m3)")
    return results


u = mda.Universe('50.50_1us.gro', 'last_ns.xtc')
results = interdigit(u, nbins=100, b=0, e=None)
for key1 in results.keys():
    for key2 in results[key1].keys():
        np.savetxt('interdigit_new3.' + key1 + '.' + key2 + '.dat', results[key1][key2])