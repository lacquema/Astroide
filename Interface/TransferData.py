#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import numpy as np


### --- Plot Window Generating --- ###

class TransferDataClass():

    def __init__(self):
        super().__init__()

    # Open followbodies.dat data
    def OpenFollowbodies(FileFollowbodies):

        # Header
        HeaderFollowbodies = np.loadtxt(FileFollowbodies, max_rows = 1, dtype = int)
        NbSteps, NbParams, NbBodies_f = HeaderFollowbodies

        # Data
        DataFollowbodies = np.loadtxt(FileFollowbodies, skiprows = 1, max_rows = 2).reshape(NbBodies_f, NbParams, NbSteps)

        t_f, a_f, e_f, i, W, w, M = [np.zeros((NbBodies_f, NbSteps)) for k in range(7)]

        for j in range(NbBodies_f):
            t_f[j], a_f[j], e_f[j], i[j], W[j], w[j], M[j] = [DataFollowbodies[j][k][:] for k in range(7)]


        return NbSteps, NbBodies_f, t_f, a_f, e_f, i, W, w, M
    

                    
    # Open mextract.dat data
    def OpenMextract(FileMextract):

        # Header
        HeaderMextract = np.loadtxt(FileMextract, max_rows = 1, dtype = int)
        NbLines, NbColumns = HeaderMextract[:2].astype(int)

        # Data
        DataMextract = np.loadtxt(FileMextract, skiprows = 1, max_rows = 2).reshape(NbColumns, NbLines)

        NbSnapshots, NbBodies0, NbParticles0 = DataMextract[6][:3].astype(int)
        
        t_m, NbBodies_m, NbParticles = [np.zeros((NbSnapshots)).astype(int) for k in range(3)]
        
        a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z = [[] for k in range(11)]

        R = []

        for j in range(NbSnapshots):
            if j == 0:
                indexLine = 2*NbBodies0
            else:
                indexLine += 1 + NbBodies_m[j-1] + NbParticles[j-1]

            t_m[j] = DataMextract[0][indexLine]
            NbBodies_m[j] = DataMextract[1][indexLine]
            NbParticles[j] = DataMextract[2][indexLine]

            a_m.append(DataMextract[0][indexLine+1:indexLine+1+NbBodies_m[j]+NbParticles[j]])
            e_m.append(DataMextract[1][indexLine+1:indexLine+1+NbBodies_m[j]+NbParticles[j]])
            Ex.append(DataMextract[2][indexLine+1:indexLine+1+NbBodies_m[j]])
            Ey.append(DataMextract[3][indexLine+1:indexLine+1+NbBodies_m[j]])
            Ez.append(DataMextract[4][indexLine+1:indexLine+1+NbBodies_m[j]])
            Epx.append(DataMextract[5][indexLine+1:indexLine+1+NbBodies_m[j]])
            Epy.append(DataMextract[6][indexLine+1:indexLine+1+NbBodies_m[j]])
            Epz.append(DataMextract[7][indexLine+1:indexLine+1+NbBodies_m[j]])
            X.append(DataMextract[8][indexLine+1:indexLine+1+NbBodies_m[j]+NbParticles[j]])
            Y.append(DataMextract[9][indexLine+1:indexLine+1+NbBodies_m[j]+NbParticles[j]])
            Z.append(DataMextract[10][indexLine+1:indexLine+1+NbBodies_m[j]+NbParticles[j]])

            R.append(DataMextract[6][indexLine+1+NbBodies_m[j]+1:indexLine+1+NbBodies_m[j]+1+NbParticles[j]])
            
        
        return NbSnapshots, t_m, NbBodies_m, NbParticles, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R

if __name__=="__main__":
    TransferDataClass.OpenFollowbodies('/Users/lacquema/ByeGildas/Data/followbodies.dat')
    TransferDataClass.OpenMextract('/Users/lacquema/ByeGildas/Data/mextract.dat')