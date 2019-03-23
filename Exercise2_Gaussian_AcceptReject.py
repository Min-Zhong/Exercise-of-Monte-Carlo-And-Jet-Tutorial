import random
import math

GaussianList = []
mu = 5
sigma = 2
iran = 0
while iran < 1000:
    x = random.random()
    x = 10*x
    y = random.random()
    fx = math.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))/(math.sqrt(2*math.pi)*sigma)
    if y > fx :
        continue
    elif x < 0: 
        continue
    elif x > 10:
        continue
    GaussianList.append(x)
    iran = iran+1

import matplotlib.pyplot as plt

def draw_hist(myList, Title, Xlabel, Ylabel, Xmin, Xmax, Ymin, Ymax):
    plt.hist(myList, 50)
    plt.xlabel(Xlabel)
    plt.xlim(Xmin, Xmax)
    plt.ylabel(Ylabel)
    plt.ylim(Ymin, Ymax)
    plt.title(Title)
    plt.show()

draw_hist(GaussianList, 'Gaussian Distribution', 'RandomNumber', 'Counts', 0.0, 10.0, 0.0, 50.0)
