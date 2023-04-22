import numpy as np
import matplotlib.pyplot as plt
from random import randint

from tipe_potent import potent_inv as potent
from tipe_poids import poids_simple as poids_pot

importance_deplacement = 0#coût d'un déplacement dans dijkstra
importance_potent = 1

from time import time

def distance(A,B):
    return ((A[0]-B[0])**2 + (A[1] - B[1])**2)**(1/2)

def calcul_poids(dim,repulsifs,murs):# ne pas oublier de mettre les contours de l'expérience dans "murs" au préalable (y compris dans le abscisses et ordonnées -1)
#murs et repulsifs sont des ensembles
    poids = [[[[importance_potent*poids_pot((i,j),(i+l-1,j+c-1),repulsifs,murs,potent) + importance_deplacement*distance((i,j),(i+l-1,j+c-1)) for c in range(3)] for l in range(3)] for j in range(dim)] for i in range(dim)]
    return poids

def dijkstra(dim,depart,arrivees,poids):#"arrivees" est un ENSEMBLE de points
    poids_tot = [[np.inf for j in range(dim+1)] for i in range(dim+1)]# dim+1 pour inclure les contours de l'expérience et ne pas avoir un index out of range (l.32)
    chemins = [[None for j in range(dim)] for i in range(dim)]

    poids_tot[depart[0]][depart[1]] = 0
    pt = depart
    poids_atteint = 0
    chemin = [pt]
    pts_passes = {pt}

    while pt not in arrivees:
        x,y = pt[0]-1, pt[1]-1
        for l in range(3):
            for c in range(3):
                if poids_tot[x+l][y+c] >  poids_atteint + poids[x+1][y+1][l][c]:
                    poids_tot[x+l][y+c] = poids_atteint + poids[x+1][y+1][l][c]
                    chemins[x+l][y+c] = chemin + [(x+l,y+c)]
        poids_min = np.inf
        pt_min = None
        for i in range(dim):
            for j in range(dim):
                if (i,j) not in pts_passes:
                    if poids_tot[i][j] < poids_min:
                        poids_min = poids_tot[i][j]
                        pt_min = (i,j)
        pt = pt_min
        pts_passes.add(pt)
        poids_atteint = poids_min
        chemin = chemins[pt[0]][pt[1]]

    return chemin

def test_show_dijkstra(dim,depart,arrivees,repulsifs,listes_murs):
    """
    dim = 3
    repulsifs = {(1,1)}
    listes_murs = [[(i,-1) for i in range(-1,4)] + [(3,j) for j in range(4)] + [(i,3) for i in range(2,-2,-1)] + [(-1,j) for j in range(2,-2,-1)]]
    depart = (0,0)
    arrivees = {(2,2)}

    dim = 100
    repulsifs = set([(i,50) for i in range(40)] + [(i,50) for i in range(60,100)])
    listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]]
    depart = (0,0)
    arrivees = {(99,99)}
    """

    start = time()

    murs = set()
    abscisses_traces_murs = []
    ordonnees_traces_murs = []
    i = 0
    for liste_murs in listes_murs:
        if len(liste_murs) > 1:
            abscisses_traces_murs.append([])
            ordonnees_traces_murs.append([])
            for m in liste_murs:
                murs.add(m)
                abscisses_traces_murs[i].append(m[0])
                ordonnees_traces_murs[i].append(m[1])
                plt.plot(abscisses_traces_murs[i],ordonnees_traces_murs[i],'b')
            i += 1
        else:
            plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')

    map_potent = np.array([[potent((i,j),repulsifs,murs) for i in range(dim)] for j in range(dim)])
    plt.imshow(map_potent,cmap='hot')

    plt.plot(depart[0],depart[1],'go')

    for a in arrivees:
        plt.plot(a[0],a[1],'ro')

    poids = calcul_poids(dim,repulsifs,murs)
    chemin = dijkstra(dim,depart,arrivees,poids)

    abscisses_chemin = []
    ordonnees_chemin = []
    for pt in chemin:
        abscisses_chemin.append(pt[0])
        ordonnees_chemin.append(pt[1])
    plt.plot(abscisses_chemin,ordonnees_chemin,'g')

    print(time()-start)

    plt.show()