# This program is used to produce a 3-dimension jet event which contains two showers,
# and then use the clustering algorithm to reconstruct these final state particles into jets

################################
## Jet Production Starts Here ##
################################

import random
import math
import matplotlib.pyplot as plt
import numpy as np

Delta_t = 0.5   # Here we assume that the time between a splitting and another is a constant
Ecrit = 10      # Particles with energy below Ecrit will not split
alpha = 0.01
# The actual energy of one of the emitted partons is random and distributed according to an 
# exponentially falling distribution Exp(−α*Eparton) where α is a physically meaningful parameter 
# that characterizes collisions at the LHC and Eparton is the randomly chosen energy of the parton
            

# Set loop times
JetNumberList = []
JetMassList = []
Ntime = 500

# Set the values of n and R
n = -1.0
R = 0.5

class Particle(object):
    
    def __init__(self, Px, Py, Pz, Xi, Yi, Zi, theta, phi):
        self._Px = Px
        self._Py = Py
        self._Pz = Pz
        self._Pt = math.sqrt(Py*Py + Pz*Pz)
        self._E = math.sqrt(Px*Px + Py*Py + Pz*Pz)
        self._theta = theta
        self._phi = phi
        self._M = 1.0
        self._Xi = Xi
        self._Yi = Yi
        self._Zi = Zi
        self._s = True   # s is a parameter which indicate that the particle have already splitted or not
    
    def Split(self):
        # The final position of a particle is:
        self._Xf = self._Xi + self._Px/self._M * Delta_t
        self._Yf = self._Yi + self._Py/self._M * Delta_t
        self._Zf = self._Zi + self._Pz/self._M * Delta_t
        # Here we set PDF(z) = 1/(1+z), 0 <= z <= 1, and theta has PDF(theta) = 1 / theta, where theta is between [0.001, Pi]   
        iran = 0
        while iran < 1:
            z = random.random()
            u = random.random()
            if u < 1/(1+z):
                self._z = z
                iran = iran+1
        iran = 0
        while iran < 1:
            newtheta = random.random()*math.pi
            if newtheta < 0.001:
                continue
            v = random.random()*1000
            if v < 1/newtheta:
                self._newtheta = newtheta
                iran = iran+1
        # Phi has a uniform distribution
        newphi = random.random()
        newphi = newphi * 2 * math.pi
        self._newphi = newphi
        self._s = False   # This means that this particle have already splitted into two particles

def ParAdd(par1, par2):
    newpx = par1._Px+par2._Px
    newpy = par1._Py+par2._Py
    newpz = par1._Pz+par2._Pz
    newE = math.sqrt(newpx*newpx + newpy*newpy + newpz*newpz)
    newtheta = np.arccos(newpx/newE)
    newphi = np.arctan(newpy/newpz)
    par12 = Particle(newpx, newpy, newpz, 0.0, 0.0, 0.0, newtheta, newphi)
    return par12

def PartonGenerate():
    # Firstly let's get a Eparton
    iE0 = 0
    while iE0 < 1:
        u = random.random()
        E0 = -math.log(math.e, (1.0-u))/alpha
        if E0 >= 300:
            continue
        iE0 = iE0+1
    theta0 = random.random()*math.pi
    phi0 = random.random()*math.pi*2.0    # According to the assumption, theta and phi are both distributed evenly
    
    Pz = E0*math.cos(theta0)
    Px = E0*math.sin(theta0)*math.cos(phi0)
    Py = E0*math.sin(theta0)*math.sin(phi0)
    X0 = 0; Y0 = 0; Z0 = 0;
    Par1 = Particle(Px, Py, Pz, X0, Y0, Z0, theta0, phi0)
    Par2 = Particle(-Px, -Py, -Pz, X0, Y0, Z0, -theta0, 2*math.pi-phi0)
    return Par1, Par2

def DeltaCalculate(par1, par2):    # Calculate Delta_ij
    
    Delta_theta = abs(par1._theta-par2._theta)
    if Delta_theta > math.pi:
        Delta_theta = Delta_theta-math.pi
    Delta_phi = abs(par1._phi-par2._phi)
    if Delta_phi > 2*math.pi:
        Delta_phi = Delta_phi-2*math.pi
               
    Delta_ij = math.sqrt(pow(Delta_theta, 2) + pow(Delta_phi, 2))  
    return Delta_ij

def draw_hist(myList, Title, Xlabel, Ylabel, Nbins, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, bins=Nbins)
    plt.xlabel(Xlabel)
    plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    #plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()
    
def draw_hist_test(myList, Title, Xlabel, Ylabel, Nbins, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, bins=Nbins)
    plt.xlabel(Xlabel)
    #plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    #plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()

#####################################
## Loop Shower and Reconstruction ###
#####################################

for iloop in range(Ntime):
    if iloop % 50 == 0:
        print (iloop, " finished")
    # Now we can generate two initial particles
    ParticleList = []
    Par01, Par02 = PartonGenerate()
    # print(iloop, Par01._E)
    ParticleList.append(Par01)
    ParticleList.append(Par02)

    # The splitting starts here
    while True:
        for pars in ParticleList:
            if pars._E >= Ecrit and pars._s == True:
                pars.Split()
            
                # The first particle's information
                Px = pars._E*pars._z*math.sin(pars._theta+pars._newtheta)*math.cos(pars._phi+pars._newphi)
                Py = pars._E*pars._z*math.sin(pars._theta+pars._newtheta)*math.sin(pars._phi+pars._newphi)
                Pz = pars._E*pars._z*math.cos(pars._theta)
                newpar1 = Particle(Px, Py, Pz, pars._Xf, pars._Yf, pars._Zf, pars._theta+pars._newtheta, pars._phi+pars._newphi)
                ParticleList.append(newpar1)
            
                # The second particle's information
                Px = pars._E*(1-pars._z)*math.sin(pars._theta-pars._newtheta)*math.cos(pars._phi-pars._newphi)
                Py = pars._E*(1-pars._z)*math.sin(pars._theta-pars._newtheta)*math.sin(pars._phi-pars._newphi)
                Pz = pars._E*(1-pars._z)*math.cos(pars._theta)
                newpar2 = Particle(Px, Py, Pz, pars._Xf, pars._Yf, pars._Zf, pars._theta-pars._newtheta, pars._phi-pars._newphi)
                ParticleList.append(newpar2)
            
        # If all the particles' energy are below Ecrit, the splitting will stop
        ipars = 0
        for pars in ParticleList:
            if pars._E < Ecrit or pars._s == False:
                ipars = ipars+1
        if ipars >= len(ParticleList):
            break

###############################
### Jet Production Ends Here ##
###############################

###############################
## Jet Reconstruction Starts ##
###############################

    # Find all the final state particles
    FinalParList = []
    # Prepare for the constituents of the reconstructed jets
    JetConsList = []

    JetList = []
    Njet = 0
    JetSubList = []
    for pars in ParticleList:
        if pars._s == True:
            FinalParList.append(pars)
            JetConsList.append([pars])
        
    N0 = len(FinalParList)
    N = N0

    # Now let's reconstruct the particles and count the number of jets
    while N > 0:
        Dij = np.zeros((N, N))
        DiB = np.zeros(N)

        i = 0
        for pars1 in FinalParList:
            DiB[i] = pow(pars1._Pt, 2*n)             # Calculate D_iB
            j = 0
            for pars2 in FinalParList:                
                Delta_ij = DeltaCalculate(pars1, pars2)
                Dij[i,j] = min(pow(pars1._Pt, 2*n), pow(pars2._Pt, 2*n)) * Delta_ij / R     # Calculate D_ij
                j = j+1
            
            i = i+1

        min_dij = Dij.max()
        min_dib = DiB.max()
        min_dij_i = 0; min_dij_j = 0
        min_dib_i = 0

        # Find the minimum of D_ij and D_iB
        for irow in range(N):
            for icol in range(N):
                if irow!=icol and Dij[irow, icol] < min_dij:
                    min_dij = Dij[irow, icol]
                    min_dij_i = irow
                    min_dij_j = icol
        for icol in range(N):
            if DiB[icol] < min_dib:
                min_dib = DiB[icol]
                min_dib_i = icol

        if N != 2:                  # If there are only two particles left, the situation will be different
            if min_dij < min_dib:     
                NewParticle = ParAdd(FinalParList[min_dij_i], FinalParList[min_dij_j])        # The two vectors are added here
                FinalParList[min_dij_i] = NewParticle
                JetConsList[min_dij_i].extend(JetConsList[min_dij_j])       # The two constituents of a jet are collected here
                del JetConsList[min_dij_j]
                del FinalParList[min_dij_j]
            else:
                # print (min_dij, min_dib)
                JetList.append(FinalParList[min_dib_i])                     # A jet is counted here
                JetSubList.append(JetConsList[min_dib_i])
                del JetConsList[min_dib_i]
                del FinalParList[min_dib_i]
                Njet = Njet+1
        elif N == 2:
            if min_dij < min_dib:
                JetList.append(ParAdd(FinalParList[min_dij_i], FinalParList[min_dij_j]))
                JetConsList[min_dij_i].extend(JetConsList[min_dij_j])
                JetSubList.append(JetConsList[min_dij_i])
                del JetConsList[min_dij_i]
                del JetConsList[min_dij_j]
                del FinalParList[min_dij_i]
                del FinalParList[min_dij_j]
                Njet = Njet+1
            else:
                min_dij_i = 1
                min_dij_j = 0
                JetList.append(FinalParList[min_dij_i])
                JetList.append(FinalParList[min_dij_j])
                JetSubList.append(JetConsList[min_dij_i])
                JetSubList.append(JetConsList[min_dij_j])
                del JetConsList[min_dij_i]
                del JetConsList[min_dij_j]
                del FinalParList[min_dij_i]
                del FinalParList[min_dij_j]
                Njet = Njet+2
        N = len(FinalParList)

    
    ## Calculate Simple Mass ##

    max_E = 0
    max_E_i = 0
    ijet = 0
    for jets in JetList:             # The jet with the highest energy is found here
        if max_E < jets._E:
            max_E = jets._E
            max_E_i = ijet
        ijet = ijet+1                

    Jetcons1 = JetSubList[max_E_i][0]     # Then we can find the two constituents with highest Pt
    ConsNumber1 = 0
    icons = 0
    for pars in JetSubList[max_E_i]:
        if Jetcons1._Pt < pars._Pt:
            Jetcons1 = pars
            ConsNumber1 = icons
        icons = icons+1
    
    Jetcons2 = JetSubList[max_E_i][0]
    ConsNumber2 = 0
    icons = 0
    for pars in JetSubList[max_E_i]:
        if Jetcons2._Pt < pars._Pt and icons != ConsNumber1:
            Jetcons2 = pars
            ConsNumber2 = icons
        icons = icons+1

    Mass = Jetcons1._E * Jetcons2._E * DeltaCalculate(Jetcons1, Jetcons2)
  
    JetNumberList.append(Njet)
    JetMassList.append(Mass)

JetMassListNew = []
for i in range(Ntime):
    if JetMassList[i] < 700.0:
        JetMassListNew.append(JetMassList[i])

draw_hist(JetNumberList, 'Number of Jet Distribution', 'Number of Jet', 'Times', 40, 0, 80, 0.0, 50.0)
draw_hist(JetMassListNew, 'Jet "Pseudo Mass" Distribution', 'Pseudo Mass', 'Times', 35, 0, 700, 0.0, 300.0)
#draw_hist(JetMassList, 'Jet "Pseudo Mass" Distribution', 'Pseudo Mass', 'Times', 35, 0, 700, 0.0, 300.0)
