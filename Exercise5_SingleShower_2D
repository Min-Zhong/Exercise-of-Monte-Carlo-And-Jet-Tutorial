# This program is used to produce a 2-dimension parton shower with an input energy of initial state particle
import random
import math
import matplotlib.pyplot as plt
import numpy as np

Delta_t = 0.5   # Here we assume that the time between a splitting and another is a constant
Ecrit = 20      # Particles with energy below Ecrit will not split

class Particle(object):
    
    def __init__(self, Px, Py, Xi, Yi, theta):
        self._Px = Px
        self._Py = Py
        self._E = math.sqrt(Px*Px + Py*Py)
        self._M = 1.0
        self._Xi = Xi
        self._Yi = Yi
        self._theta = theta
        self._s = True   # s is a parameter which indicate that the particle have already splitted or not
    
    def Split(self):
        # The final position of a particle is:
        self._Xf = self._Xi + self._Px/self._M * Delta_t
        self._Yf = self._Yi + self._Py/self._M * Delta_t
        # Here we set PDF(z) = 1/(1+z), 0 <= z <= 1, and theta has the PDF = 1/theta where theta is between [0.001, Pi/2]    
        u = random.random()*np.log(2)
        z = np.exp(u)-1
        pars._z = z
        iran = 0
        while iran < 1:
            newtheta = random.random()*math.pi/2.0
            v = random.random()*1000
            if newtheta < 0.001:
                continue
            if v < 1/newtheta:
                self._newtheta = newtheta
                iran = iran+1
        self._s = False   # This means that this particle have already splitted into two particles
        
            
user_input = input("Enter the initial particle with a Energy(GeV)\nPlease Input here: ")
inputline = user_input

ParticleList = []
E = float(inputline)
Px = E; Py = 0.0
x = 0.0; y = 0.0
theta0 = 0.0
Par0 = Particle(Px, Py, x, y, theta0)
ParticleList.append(Par0)

# The splitting starts here
while True:
    for pars in ParticleList:
        if pars._E > Ecrit and pars._s == True:
            pars.Split()
            # The first particle's information
            Px = pars._E*pars._z*math.cos(pars._theta)
            Py = pars._E*pars._z*math.sin(pars._theta+pars._newtheta)
            newpar1 = Particle(Px, Py, pars._Xf, pars._Yf, pars._theta+pars._newtheta)
            ParticleList.append(newpar1)
            
            # The second particle's information
            Px = pars._E*(1-pars._z)*math.cos(pars._theta)
            Py = pars._E*(1-pars._z)*math.sin(pars._theta-pars._newtheta)
            newpar2 = Particle(Px, Py, pars._Xf, pars._Yf, pars._theta-pars._newtheta)
            ParticleList.append(newpar2)
                    
    # If all the particles' energy are below Ecrit, the splitting will stop
    ipars = 0
    for pars in ParticleList:
        if pars._E < Ecrit or pars._s == False:
            ipars = ipars+1
    if ipars >= len(ParticleList):
        break
        
# Drawing the plot of the shower
for pars in ParticleList:
    if pars._E > Ecrit:
        plt.plot([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], color='r')
        plt.scatter([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], color='b')
    if pars._E < Ecrit:
        pars._Xf = pars._Xi + pars._Px*4*Delta_t      # In fact, here is n    ot Splitting, it's just for calculating Xf and Yf
        pars._Yf = pars._Yi + pars._Py*4*Delta_t
        plt.plot([pars._Xi, pars._Xf], [pars._Yi, pars._Yf], color='r')
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
