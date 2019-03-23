import random
import math
import numpy as np 
import matplotlib.pyplot as plt

def draw_hist(myList, Title, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, 30)
    plt.xlabel(Xlabel)
    #plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    #plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()

## First: Calculate the value of Pi for 1000 times
PiValueList = []
for i in range(1000):
    Naccept = 0
    Ntotal = 0
    
    for iran in range(1000):
        x = random.random()
        y = random.random()
        if x*x + y*y <= 1.0:
            Naccept = Naccept+1
        Ntotal = Ntotal+1
        
    PiValueList.append(4*Naccept/Ntotal)

PiMean = np.mean(PiValueList)
PiSigma = np.std(PiValueList, ddof=1)

PiMean = (PiMean//0.001)/1000
PiSigma = (PiSigma//0.0001)/10000

draw_hist(PiValueList, u"µ" + " = " + str(PiMean)+", "+ u"σ" + " = " + str(PiSigma), 'Value', 'Counts', 0.0, 6, 0.0, 35.0)

## Second: Increase the times of calculation to see how 𝜎 changes (every time with 1000 points)
SigmaList1 = []
NtimeList = [10, 20, 40, 80, 160, 240, 320, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500]
for i in NtimeList:
    Ntimes = i
    PiValueListNew = []

    for j in range(Ntimes):
        Naccept = 0
        Ntotal = 0
        for k in range(1000):
            x = random.random()
            y = random.random()
            if x*x + y*y <= 1.0:
                Naccept = Naccept+1
            Ntotal = Ntotal+1
        
        PiValueListNew.append(4*Naccept/Ntotal)
        
    PiSigma = np.std(PiValueListNew, ddof=1)
    print(Ntimes,PiSigma, len(PiValueListNew))
    SigmaList1.append(PiSigma)

plt.plot(NtimeList, SigmaList1)
plt.xlabel("Times of Calculation")
plt.ylabel(u"σ")
plt.title("")
plt.show()

## Third: Increase the number of points to see how 𝜎 changes (Calculation times remains 1000)
SigmaList2 = []
NumofPointList = []
for i in range(10):
    NumofPoint = (i+1)*300
    PiValueListNew = []
    NumofPointList.append(NumofPoint)

    for j in range(1000):
        Naccept = 0
        Ntotal = 0
        for k in range(NumofPoint):
            x = random.random()
            y = random.random()
            if x*x + y*y <= 1.0:
                Naccept = Naccept+1
            Ntotal = Ntotal+1
        
        PiValueListNew.append(4*Naccept/Ntotal)
        
    PiSigma = np.std(PiValueListNew, ddof=1)
    print(NumofPoint,PiSigma, len(PiValueListNew))
    SigmaList2.append(PiSigma)

plt.plot(NumofPointList, SigmaList2)
plt.xlabel("Number of Random Points")
plt.ylabel(u"σ")
plt.title("")
plt.show()
