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
