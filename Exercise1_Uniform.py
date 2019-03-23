import random

NumberList1 = []
for iran in range(1000):
    ran1 = random.random()
    NumberList1.append(ran1)

import matplotlib.pyplot as plt

def draw_hist(myList, Title, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, 50)
    plt.xlabel(Xlabel)
    plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()

draw_hist(NumberList1, 'Uniform Distribution between [0,1]', 'RandomNumber', 'Counts', 0.0, 1.0, 0.0, 35.0)

NumberList2 = []
for iran in range(1000):
    ran2 = random.random()
    ran3 = ran2*10 + 5.0
    NumberList2.append(ran3)

draw_hist(NumberList2, 'Uniform Distribution between [5,15]', 'RandomNumber', 'Counts', 5.0, 15.0, 0.0, 35.0)
