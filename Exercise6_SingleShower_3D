# This program is used to produce a 3-dimension parton shower with an input energy of the initial state particle
import random
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Delta_t = 0.5   # Here we assume that the time between a splitting and another is a constant
Ecrit = 20      # Particles with energy below Ecrit will not split

class Particle(object):
    
    def __init__(self, Px, Py, Pz, Xi, Yi, Zi, theta, phi):
        self._Px = Px
        self._Py = Py
        self._Pz = Pz
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
        # Here we set PDF(z) = 1/（1+z）, 0 <= z <= 1, and theta has the PDF = 1/theta where theta is between[0.001, Pi/2]   
        iran = 0
        while iran < 1:
            z = random.random()
            u = random.random()
            if u < 1/(1+z):
                self._z = z
                iran = iran+1
        iran = 0
        while iran < 1:
            newtheta = random.random()*math.pi/2.0
            v = random.random()*1000
            if newtheta < 0.001:
                continue
            if v < 1/newtheta:
                self._newtheta = newtheta
                iran = iran+1
        # Phi has a uniform distribution
        newphi = random.random()
        newphi = newphi * 2 * math.pi
        self._newphi = newphi
        self._s = False   # This means that this particle have already splitted into two particles
        
            
user_input = input("Enter the initial particle with a Energy(GeV)\nPlease Input here: ")
inputline = user_input

ParticleList = []
E = float(inputline)
Px = 0.0; Py = 0.0; Pz = E
x = 0.0; y = 0.0; z = 0.0
Par0 = Particle(Px, Py, Pz, x, y, z, 0.0, 0.0)
ParticleList.append(Par0)

# The splitting starts here
while True:
    for pars in ParticleList:
        if pars._E > Ecrit and pars._s == True:
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
        
# Drawing the 3-D plot of the shower
fig = plt.figure()
ax1 = plt.axes(projection='3d')

for pars in ParticleList:
    if pars._E > Ecrit:
        ax1.plot3D([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], [pars._Zi, pars._Zf], color='r')
        ax1.scatter3D([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], [pars._Zi, pars._Zf], color='b')
    if pars._E < Ecrit:
        pars._Xf = pars._Xi + pars._Px*8*Delta_t      # This step is just for showing the final state particles' derection
        pars._Yf = pars._Yi + pars._Py*8*Delta_t
        pars._Zf = pars._Zi + pars._Pz*8*Delta_t
        ax1.plot3D([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], [pars._Zi, pars._Zf], color='r')
        
ax1.scatter3D(0.0, 0.0, 0.0, color = 'y')
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
