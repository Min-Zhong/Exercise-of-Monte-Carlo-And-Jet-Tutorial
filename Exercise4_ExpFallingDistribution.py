# This program is used to produce random numbers distributed according to a falling distribution 
# of the form exp(−x) in the range[0,5] using Inverse Transform method
import random
import math
import numpy as np

Naccept1 = 0
Ntotal1 = 0
ExpNumberList1 = []
lamda = 1.0

for iran in range(10000):
    u = random.random()
    x = -np.log(1.0-u)/lamda
    if x > 5:
        Ntotal1 = Ntotal1+1
    else:
        ExpNumberList1.append(x)
        Naccept1 = Naccept1+1
        Ntotal1 = Ntotal1+1

import matplotlib.pyplot as plt

def draw_hist(myList, Title, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, 50)
    plt.xlabel(Xlabel)
    plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()

draw_hist(ExpNumberList1, 'Inverse Transform Method', 'RandomNumber', 'Counts', 0.0, 5.0, 0.0, 1100.0)

# This program is used to produce random numbers distributed according to a falling distribution 
# of the form exp(−x) in the range[0,5] using Accept-Reject method

lamda = 1.0
Naccept2 = 0
Ntotal2 = 0
ExpNumberList2 = []

for iran in range(10000):
    x = random.random()
    x = 5*x
    y = random.random()
    if y > lamda* math.exp(-lamda*x):
        Ntotal2 = Ntotal2+1
    elif y < lamda* math.exp(-lamda*x):
        ExpNumberList2.append(x)
        Naccept2 = Naccept2+1
        Ntotal2 = Ntotal2+1
    
draw_hist(ExpNumberList2, 'Accept-Reject Method', 'RandomNumber', 'Counts', 0.0, 5.0, 0.0, 200.0)

# So we can compare the efficiency of the two different method now
print("The efficiency of inverse-transform method: ", Naccept1/Ntotal1)
print("The efficiency of accept-reject method: ", Naccept2/Ntotal2)

# We can see that the efficiency of inverse-transform method is much larger than that of accept-reject method. Note that 
# the efficiency of inverse-transform method is supposed to be 100%, but since we choose the region[0,5], some numbers are 
# abandoned, so the efficiency is not 100%.
