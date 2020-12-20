import SMDAnalysis as smda
import density_help
import MDAnalysis as mda
import time
import numpy as np
import os
from   SMDAnalysis.common.frame import Frame

mass = {'H': 1.008, 'O': 15.999, 'C': 12.011, 'N': 14.007, 
        'P': 30.976, 'S':32.06, 'SOD':22.99, 'CLA':35.45}

def running_average(array, N=100):
    t = array[:,0]
    x = array[:,1]
    new_x = np.convolve(x, np.ones((N,))/N, mode='valid')
    new_t = t[:len(new_x)]
    new_array = np.transpose(np.array([new_t, new_x]))
    return new_array

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

def get_mass(group):
    m = []
    for atom in group:
        if atom.name == 'SOD' or atom.name == 'CLA':
            name = atom.name
        else:
            name = atom.name[0]
        m.append(mass[name])
    m = np.array(m)
    return m


def interdigit(u, nbins=100, nblocks=5, b=0, e=10000000):
    d = smda.Density()
    bframe, eframe = Frame().frame(u, b, e)
    print("bframe: %d, eframe: %d" %(bframe, eframe))

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
    
    masses = {'memb':  get_mass(groups['memb']), 
              'umemb': get_mass(groups['umemb']),
              'lmemb': get_mass(groups['lmemb']),
              'trio':  get_mass(groups['trio'])}
    
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
    d2_densities = [] #MEMB-TRIO
    d3_densities = [] #CORE-TRIO

    total_ov = []
    strong_ov = []
    weak_ov = []
    
    strong_num   = []
    

    for ts in u.trajectory[bframe:eframe]:
        if int(ts.time/1000) % 1000 == 0:
            print("analyzing %d us.... " %(ts.time/1000/1000))

        pbcx, pbcy, pbcz = u.dimensions[0:3]
        dz = pbcz/nbins
        trio_pos = groups['trio'].positions
    
        utz = np.average(groups['umemb'].positions[:,2])
        ltz = np.average(groups['lmemb'].positions[:,2])
        # print("upper tail z (A): ", utz)
        # print("lower tail z (A): ", ltz)
        
        ### d0: density of all membrane atoms
        d0 = d.density_frame(pbcx, pbcy, pbcz, nbins, groups['memb'].positions, masses['memb'])
        d0_densities.append(d0)
        
        ### TOTAL INTERDIGITATION
        ### d1: density of all trio atoms
        d1 = d.density_frame(pbcx, pbcy, pbcz, nbins, trio_pos, masses['trio'])
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
        mm = masses['trio'][boolArray]
        ### d2: density of MEMB-TRIO
        d2 = d.density_frame(pbcx, pbcy, pbcz, nbins, pp, mm)
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
        mm = masses['trio'][boolArray]
        ### d3: density of CORE-TRIO 
        d3 = d.density_frame(pbcx, pbcy, pbcz, nbins, pp, mm)
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
        zs.append(pbcz/10)
    
    strong_num = np.array(strong_num)
    XX = np.linspace(0, np.average(zs), nbins)
    results = {}

    results['inter']  = {}
    results['inter']['total']  = np.transpose([times, total_inter])
    results['inter']['strong'] = np.transpose([times, strong_inter])
    results['inter']['weak']   = np.transpose([times, weak_inter])
    
    results['ov'] = {}
    results['ov']['total']     = np.insert(smda.Block().block(total_ov, nblocks), 0, XX, axis=0).T
    results['ov']['strong']    = np.insert(smda.Block().block(strong_ov, nblocks), 0, XX, axis=0).T
    results['ov']['weak']      = np.insert(smda.Block().block(weak_ov, nblocks), 0, XX, axis=0).T
 
    results['ratio'] = {}
    results['ratio']['num']             = np.transpose([times, strong_num])
    results['ratio']['trio-to-pl']      = np.transpose([times, strong_num/numP])
    results['ratio']['trio-to-pl+trio'] = np.transpose([times, strong_num/(numP + strong_num)])
    
    results['density'] = {}
    results['density']['PL']        = np.insert(smda.Block().block(d0_densities, nblocks), 0, XX, axis=0).T
    results['density']['TRIO']      = np.insert(smda.Block().block(d1_densities, nblocks), 0, XX, axis=0).T
    results['density']['MEMB-TRIO'] = np.insert(smda.Block().block(d2_densities, nblocks), 0, XX, axis=0).T
    results['density']['CORE-TRIO'] = np.insert(smda.Block().block(d3_densities, nblocks), 0, XX, axis=0).T

    print("units: Z (nm), interdigitation (nm), time (ns), density (kg/m3)")
    return results


folders = ['../../TRIO_8nm/', '../../TRIO_16nm/']

for f in folders:
    basename = os.path.basename(f[:-1])
    u = mda.Universe(f + 'md.gro', f + 'md.1ns.xtc')
    print('analyzing %s' %basename)
    results = interdigit(u, nbins=100, nblocks=5, b=5000)

    for key1 in results.keys():
        for key2 in results[key1].keys():
            np.savetxt(basename + '.' + key1 + '.' + key2 + '.dat', results[key1][key2])




