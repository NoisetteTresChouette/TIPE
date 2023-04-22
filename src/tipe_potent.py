import numpy as np
from math import log
import matplotlib.pyplot as plt

def distance(A,B):
    return ((A[0]-B[0])**2 + (A[1] - B[1])**2)**(1/2)

def potent_inv(pt,repulsifs,murs):
    if pt in murs:
        return np.inf
    V = 0
    for r in repulsifs:
        V += 1/(1+distance(pt,r))
    return V

def potent_dist(pt,repulsifs,murs):
    if pt in murs:
        return np.inf
    V = 0
    for r in repulsifs:
        V -= distance(pt,r)
    return V

def potent_plateau(pt,repulsifs,murs):
    if pt in murs:
        return np.inf
    else:
        return 1

def ex_potent_dist():
    x = np.linspace(0,1,100)
    def dist(a,b):
        return abs(a-b)
    y = dist(x,0.5)

    plt.plot(0.5,0,'go')
    plt.text(0.56,0,'arrivée',color = 'g', size='x-large')

    plt.plot(x,y)
    plt.show()

def ex_potent_plateau():
    plt.plot([0,0.5],[0,0],'b')
    plt.plot([0.5,1.5],[1,1],'r')
    plt.plot([1.5,2],[0,0],'b')

    plt.show()

def ex_potent_inv():
    x = np.linspace(0,1,100)
    def dist(a,b):
        return abs(a-b)
    y = 1/(1+dist(x,0.5))

    plt.plot(0.5,1,'ro')

    plt.plot(x,y)
    plt.show()