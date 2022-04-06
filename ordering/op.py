import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class OrderParameters:
    def __init__(self):
        pass
    
    def compute_OP(self, u, atomlists = [], resname = False):
        assert resname, 'provide a resname'
        
        C_numbers = []
        Cs = []
        Hs = []
        repeat = []
        for atoms in atomlists:
            C_number = atoms[0][2:]
            C_numbers.append(int(C_number))
            
            Cs.append(atoms[0])
            Hs.append(atoms[1:])
            ## How many hydrogen atoms per center carbon atom
            repeat.append(len(atoms)-1) 


        Hs_f = [item for sublist in Hs for item in sublist]
        assert int(np.sum(repeat)) == len(Hs_f), "wrong in repeats"
        ## Cs =   [C32, C33, ..., C318]
        ## Hs =   [[H2X, H2Y], [H3X, H3Y], ... [H18X, H18Y, H18Z]]
        ## Hs_f = [H2X, H2Y, H3X, ... , H18Z]
        ## repeat = [2, 2, 2, ..., 3]
        
        ## select carbon atoms (g1) and hydrogen atoms (g2)
        g1 = "resname %s" %resname + " and name " + " ".join(Cs)
        g2 = "resname %s" %resname + " and name " + " ".join(Hs_f)
        group1 = u.select_atoms(g1)
        group2 = u.select_atoms(g2)
        
        natoms       = len(Cs) ## How many center atoms
        nmols        = int(len(group1.positions)/natoms) ## How many molecules
        repeats      = repeat * nmols ## [2, 2, 2, ..., 3, | 2, 2, 2, ..., 3, | ...]
        splits       = np.cumsum(repeats)

        print('# of mols: %d' %nmols)
        print('# of Carbons per molecule: %d' %natoms)

        output = []
        for ts in u.trajectory:
            p1 = np.repeat(group1.positions, repeats, axis=0)
            p2 = group2.positions
            dp = p2 - p1
            norm = np.sqrt(np.sum(np.power(dp, 2), axis=-1))
            cos_theta = dp[...,2]/norm
            S = -0.5 * (3 * np.square(cos_theta) - 1)

            new_S = self._repeat(S, repeats)
            new_S.shape = (nmols, natoms)
            results = np.average(new_S, axis=0)
            output.append(results)
        
        print(np.array(output))
        avg = np.average(output, axis=0)
        np.transpose([C_numbers, avg])
        a = np.transpose([C_numbers, avg])
        return a
        print(a)



    def _repeat(self, x, reps):
        ''' 
        x = [1, 3, 5, 7, 9]
        reps = [2, 2, 1]
        out = [2, 6, 9]
        '''
        assert len(x) == int(np.sum(reps)), 'repeats wrong'

        i = 0
        out = []
        for rep in reps:
            tmp = []
            for r in range(rep):
                tmp.append(x[i])
                i += 1
            out.append(np.average(tmp))
    
        b = np.array(out)
        return b
        print(b)




import tail_data
import MDAnalysis as mda
u = mda.Universe('1.21_100ns.gro', '1.21_100ns.xtc')
#POPC1 = OrderParameters().compute_OP(u, opc.POPC1, resname = 'POPC')
#POPC2 = OrderParameters().compute_OP(u, opc.POPC2, resname = 'POPC')
#DOPE1 = OrderParameters().compute_OP(u, opc.DOPE1, resname = 'DOPE')
#DOPE2 = OrderParameters().compute_OP(u, opc.DOPE2, resname = 'DOPE')
CHYO = OrderParameters().compute_OP(u, tail_data.CHYO, resname = 'CHYO')
#np.savetxt('POPC1.dat', POPC1)
#np.savetxt('POPC2.dat', POPC2)
#np.savetxt('DOPE1.dat', DOPE1)
#np.savetxt('DOPE2.dat', DOPE2)
np.savetxt('1.21_100ns.dat',CHYO)

p=np.array(CHYO)
p.shape
p_plot = pd.DataFrame(p, columns=['Scd',
                                  'Smectic'])
p_plot.head()


p_plot.plot(x='Scd')
plt.title('Tail Order Parameters (Scd)')
plt.ylabel('Scd') and plt.xlabel('Carbon')
plt.savefig('order_smec.png')

plt.draw()

print('***Success***')





