import matplotlib.pyplot as plt
import numpy as np
from tipe_directions_souhaitees import *

r = 1
h = 1

def vitesse_souhaitee(dim,position,directions):
    if position == None:
        return (0,0)
    coeff = int(dim/len(directions))
    x,y = position
    i = int(x//coeff)
    j = int(y//coeff)
    dx = directions[i][j][0] - i
    dy = directions[i][j][1] - j
    d = ((dx**2)+(dy**2))**(1/2)
    if d == 0:
        return (0,0)
    return (dx/d,dy/d)

def vitesses_souhaitees(dim,positions,directions):
    vitesses = []
    for i in range(len(positions)):
        vitesses.append([vitesse_souhaitee(dim,positions[i],directions)[0],vitesse_souhaitee(dim,positions[i],directions)[1]])
    return vitesses

def uzawa_initialisation(positions,listes_murs):
    global r

    n = len(positions)
    n_obst = len(listes_murs)

    D = []
    E = []
    for q1 in range(n):
        D.append([])
        E.append([])

        for q2 in range(n):
            d = distance(positions[q1],positions[q2])
            D[q1].append(d - 2*r)
            if d != 0:
                E[q1].append(((positions[q2][0] - positions[q1][0])/d , (positions[q2][1] - positions[q1][1])/d ))
            else:
                E[q1].append((0,0))

        for obst in range(n_obst):
            pt_proche = (np.inf,np.inf)
            d = distance(positions[q1],pt_proche)
            for mur in listes_murs[obst]:
                if distance(positions[q1],mur) < d:
                    pt_proche = mur
                    d = distance(positions[q1],pt_proche)
            D[q1].append(d-r)
            if d !=0:
                E[q1].append( ((pt_proche[0]-positions[q1][0])/d, (pt_proche[1]-positions[q1][1])/d ) )
            else:
                E[q1].append( (0, 0) )

    epsilon = 0.1*r
    iter_max = 5000
    rho = 1# constante sélectionnée après l'essai de plusieurs valeurs

    return D,E,n,n_obst,epsilon,iter_max,rho

def phi(v,E,n,n_obst):
    global h

    retour = []
    for i in range(n):
        retour.append([])
        for nul in range(i+1):
            retour[i].append(0)
        for j in range(i+1,n):
            G = [[0,0] for _ in range(i)] + [[-1*E[i][j][0], -1*E[i][j][1]]] + [[0,0] for _ in range(i+1,j)] + [[E[i][j][0],E[i][j][1]]] + [[0,0] for _ in range(j+1,n)]
            retour[i].append(dot(np.array(G),v))# dot correspond au produit scalaire canonique
        for l in range(n_obst):
            G = [[0,0] for _ in range(i)] + [[-1*E[i][n+l][0], -1*E[i][n+l][1]]] + [[0,0] for _ in range(i+1,n)]
            retour[i].append(dot(np.array(G),v))
    return -1*h*np.array(retour)

def dot(G,v):
    return sum([v[i,j]*G[i,j] for i in range(len(v)) for j in range(len(v[0]))])

def phi_star(v,E,n,n_obst):
    global h
    retour = np.array([[0.,0.] for i in range(n)])
    for i in range(n):
        for j in range(i+1,n):
            G = [[0,0] for _ in range(i)] + [[-1*E[i][j][0], -1*E[i][j][1]]] + [[0,0] for _ in range(i+1,j)] + [[E[i][j][0],E[i][j][1]]] + [[0,0] for _ in range(j+1,n)]
            retour -=  v[i,j] * np.array(G)
        for l in range(n_obst):
            G = [[0,0] for _ in range(i)] + [[-1*E[i][n+l][0], -1*E[i][n+l][1]]] + [[0,0] for _ in range(i+1,n)]
            retour -= v[i,n+l] * np.array(G)
    return h*retour

def projete(mu):
    mu_retour = mu.copy()
    for i in range(len(mu_retour)):
        for j in range(len(mu_retour[0])):
            if mu_retour[i,j] < 0:
                mu_retour[i,j] = 0
    return mu_retour

def uzawa(positions,u,D,E,n,n_obst,epsilon,iter_max,rho,listes_murs):
    global r
    global h

    positions_candidat = positions + h*u

    v = u.copy()
    k = 0

    liste_distances = [distance(tuple(positions_candidat[i]),tuple(positions_candidat[j])) -2*r for i in range(n) for j in range(i+1,n)]
    for i in range(n):
        for l in range(n_obst):
            liste_distances.append(min([distance(tuple(positions_candidat[i]), mur)-r for mur in listes_murs[l]]))
    Dmin = min(liste_distances)
    mu = np.array([[0 for j in range(n+n_obst)] for i in range(n)])

    while (k<iter_max) and (Dmin<-1*epsilon):
        v = u - phi_star(mu,E,n,n_obst)
        mu = projete(mu+rho*(phi(v,E,n,n_obst) - np.array(D)))

        positions_candidat = positions + h*v

        liste_distances = [distance(tuple(positions_candidat[i]),tuple(positions_candidat[j])) -2*r for i in range(n) for j in range(i+1,n)]
        for i in range(n):
            for l in range(n_obst):
                liste_distances.append(min([distance(tuple(positions_candidat[i]), mur)-r for mur in listes_murs[l]]))
        Dmin = min(liste_distances)
        k += 1

    print(Dmin)

    return positions_candidat

def mouvement_foule(dim,departs,listes_murs,directions):
    etapes_max = 1000
    nbr_individus = len(departs)
    positions = []
    etapes = []
    mouvement_fini = []
    for depart in departs:
        positions.append(depart)
    etapes.append(positions)
    vitesses_souhaitees_brutes = vitesses_souhaitees(dim,positions,directions)
    for i in range(len(positions)):
        if vitesses_souhaitees_brutes[i] == [0,0]:
            mouvement_fini.append(True)
        else:
            mouvement_fini.append(False)

    while False in mouvement_fini and len(etapes)< etapes_max:
        positions_nettoyees = []
        vitesses_souhaitees_nettes = []
        vitesses_souhaitees_brutes = vitesses_souhaitees(dim,positions,directions)
        for i in range(len(vitesses_souhaitees_brutes)):
            if vitesses_souhaitees_brutes[i] == [0,0]:
                mouvement_fini[i] = True
            else:
                positions_nettoyees.append(list(positions[i]))
                vitesses_souhaitees_nettes.append(vitesses_souhaitees_brutes[i])
        if len(positions_nettoyees)==0:
            break

        D,E,n,n_obst,epsilon,iter_max,rho = uzawa_initialisation(positions_nettoyees,listes_murs)
        vitesses_souhaitees_nettes = np.array(vitesses_souhaitees_nettes)
        positions_nettoyees = np.array(positions_nettoyees)

        nouvelles_positions = uzawa(positions_nettoyees,vitesses_souhaitees_nettes,D,E,n,n_obst,epsilon,iter_max,rho,listes_murs)

        positions = []
        j = 0
        for i in range(nbr_individus):
            if mouvement_fini[i] == True:
                positions.append(None)
            else:
                positions.append(tuple(nouvelles_positions[j]))
                j += 1

        etapes.append(positions)

    return etapes

def test_show_mouvement_foule_instants(dim,departs,listes_murs,directions,arrivees,instants):
    """
    dim = 100
    departs = {(20,20)}
    listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]] + [[(i,40) for i in range(10,29)] + [(29,j) for j in range(40,60)] + [(i,59) for i in range(29,9,-1)] + [(10,j) for j in range(59,39,-1)]] + [[(i,40) for i in range(40,65)] + [(64,j) for j in range(44,60)] + [(i,59) for i in range(64,39,-1)] + [(40,j) for j in range(59,39,-1)]] + [[(i,5) for i in range(50,65)] + [(64,j) for j in range(5,30)] + [(i,29) for i in range(64,49,-1)] + [(50,j) for j in range(29,4,-1)]] + [[(i,80) for i in range(65,80)] + [(79,j) for j in range(80,95)] + [(i,94) for i in range(79,64,-1)] + [(65,j) for j in range(94,79,-1)]]
    directions = [[(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 4), (0, 7), (0, 8), (0, 9), (0, 10), (1, 11), (1, 12), (1, 13), (1, 14), (1, 15), (1, 16), (1, 17), (1, 18), (1, 19), (1, 19)], [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 5), (1, 8), (1, 9), (1, 10), (1, 11), (2, 12), (2, 13), (2, 14), (2, 15), (2, 16), (2, 17), (2, 18), (2, 19), (2, 19)], [(3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 6), (2, 8), (2, 9), (2, 10), (2, 11), (3, 13), (3, 14), (3, 15), (3, 16), (3, 17), (3, 18), (3, 19), (3, 19)], [(4, 0), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 6), (3, 8), (3, 9), (3, 10), (3, 11), (4, 13), (4, 14), (4,15), (4, 16), (4, 17), (4, 18), (4, 19), (4, 19)], [(5, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 6), (4, 8), (4, 9), (4, 10), (4, 11), (5, 13), (5, 14), (5, 15), (5, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(6, 1), (6, 2), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 6), (5, 8), (5, 9), (5, 10), (5, 11), (6, 13), (6, 14), (6, 15), (6, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(7, 1), (7, 2), (7, 3), (7, 3), (7, 4), (7, 5), (7, 6), (7, 6), (7, 7), (7, 10), (7,11), (5, 12), (7, 13), (7, 14), (7, 15), (6, 16), (6, 17), (5, 18), (5, 19), (5, 19)], [(8, 1), (8, 2), (8, 3), (8, 4), (8, 4), (8, 5), (8, 6), (8, 7), (8, 7),(7, 10), (7, 11), (8, 12), (8, 13), (7, 14), (7, 15), (7, 16), (6, 17), (6, 18), (6, 19), (6, 19)], [(9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 5), (9, 6), (9, 7), (8, 8), (8, 9), (8, 10), (8, 11), (9, 12), (8, 14), (8, 15), (7, 16), (7, 17), (7, 18), (7, 19), (7, 19)], [(10, 0), (9, 2), (9, 3), (9, 4), (9, 5), (10, 6), (10, 6), (10, 7), (9, 8), (9, 9), (9, 10), (9, 11), (10, 12), (10, 13), (8, 15), (8, 16), (8, 17), (8, 18), (8, 19), (8, 19)], [(11, 0), (10, 1), (10, 2), (10, 3), (10, 4), (10, 5), (11, 6), (11, 7), (10, 8), (10, 9), (10, 10), (10, 11), (11, 12), (11, 13), (11, 14), (11, 15), (11, 17), (11, 18), (11, 19), (11, 19)], [(12, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5), (12, 7), (12, 7), (11,8), (11, 9), (11, 10), (11, 11), (12, 12), (12, 13), (12, 14), (12, 15), (12, 15), (12, 18), (12, 19), (12, 19)], [(13, 1), (12, 1), (12, 2), (12, 3), (12, 4), (12, 5), (13, 7), (13, 8), (12, 8), (12, 9), (12, 10), (12, 11), (13, 13), (13, 13), (13, 14), (13, 15), (13, 15), (12, 18), (13, 19), (13, 19)], [(14, 1), (14, 2), (14, 3), (14, 4), (14, 5), (14, 6), (14, 7), (14, 8), (14, 9), (14, 10), (14, 11), (14, 12), (14, 13), (14, 14), (14, 14), (14, 15), (13, 16), (13, 17), (13, 18), (14, 19)], [(15, 1), (15, 2), (15, 3), (15, 4), (15, 5), (15, 6), (15, 7), (15, 8), (15, 9), (15, 10), (15, 11), (15, 12), (15, 13), (15, 14), (15, 15), (15, 15), (14, 16), (14, 17), (14, 18), (15, 19)], [(16, 1), (16, 2), (16, 3), (16, 4), (16, 5), (16, 6), (16, 7), (16, 8), (16, 9), (16, 10), (16, 11), (16, 12), (16, 13), (16, 14), (16, 15), (16, 16), (15, 16), (15, 17), (15, 18), (16, 19)], [(17, 1), (17, 2), (17, 3), (17, 4), (17, 5), (17, 6), (17, 7), (17, 8), (17, 9), (17, 10), (17, 11), (17, 12), (17, 13), (17, 14), (17, 15), (17, 16), (17, 17), (17, 18), (17, 19), (17, 19)], [(18, 1),(18, 2), (18, 3), (18, 4), (18, 5), (18, 6), (18, 7), (18, 8), (18, 9), (18, 10), (18, 11), (18, 12), (18, 13), (18, 14), (18, 15), (18, 16), (18, 17), (18, 18), (18, 19), (18, 19)], [(19, 1), (19, 2), (19, 3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19,15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)], [(19, 1), (19, 2), (19,3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19, 15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)]]
    arrivees = {(20,99),(99,99)}
    """
    global r

    start = time()
    etapes = mouvement_foule(dim,departs,listes_murs,directions)
    print(time()-start)

    print(etapes)

    for t in instants:
        if t<len(etapes):
            print(etapes[t])

            figure,axes = plt.subplots()
            axes.set_aspect(1)
            for pt in (etapes[t]):
                if pt != None:
                    circle = plt.Circle(pt,r,facecolor = 'orange',edgecolor = 'k')
                    axes.add_artist(circle)

            abscisses_traces_murs = []
            ordonnees_traces_murs = []
            i =0
            for liste_murs in listes_murs:
                if len(liste_murs) > 1:
                    abscisses_traces_murs.append([])
                    ordonnees_traces_murs.append([])
                    for m in liste_murs:
                        abscisses_traces_murs[i].append(m[0])
                        ordonnees_traces_murs[i].append(m[1])
                        plt.plot(abscisses_traces_murs[i],ordonnees_traces_murs[i],'b')
                    i += 1
                else:
                    plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')

            for a in arrivees:
                plt.plot(a[0],a[1],'gv')

            plt.show()

def test_show_mouvement_foule_traj(dim,departs,listes_murs,directions,arrivees):
    """
    dim = 100
    departs = {(20,20),(2,80)}
    listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]] + [[(i,40) for i in range(10,29)] + [(29,j) for j in range(40,60)] + [(i,59) for i in range(29,9,-1)] + [(10,j) for j in range(59,39,-1)]] + [[(i,40) for i in range(40,65)] + [(64,j) for j in range(44,60)] + [(i,59) for i in range(64,39,-1)] + [(40,j) for j in range(59,39,-1)]] + [[(i,5) for i in range(50,65)] + [(64,j) for j in range(5,30)] + [(i,29) for i in range(64,49,-1)] + [(50,j) for j in range(29,4,-1)]] + [[(i,80) for i in range(65,80)] + [(79,j) for j in range(80,95)] + [(i,94) for i in range(79,64,-1)] + [(65,j) for j in range(94,79,-1)]]
    directions = [[(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 4), (0, 7), (0, 8), (0, 9), (0, 10), (1, 11), (1, 12), (1, 13), (1, 14), (1, 15), (1, 16), (1, 17), (1, 18), (1, 19), (1, 19)], [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 5), (1, 8), (1, 9), (1, 10), (1, 11), (2, 12), (2, 13), (2, 14), (2, 15), (2, 16), (2, 17), (2, 18), (2, 19), (2, 19)], [(3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 6), (2, 8), (2, 9), (2, 10), (2, 11), (3, 13), (3, 14), (3, 15), (3, 16), (3, 17), (3, 18), (3, 19), (3, 19)], [(4, 0), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 6), (3, 8), (3, 9), (3, 10), (3, 11), (4, 13), (4, 14), (4,15), (4, 16), (4, 17), (4, 18), (4, 19), (4, 19)], [(5, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 6), (4, 8), (4, 9), (4, 10), (4, 11), (5, 13), (5, 14), (5, 15), (5, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(6, 1), (6, 2), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 6), (5, 8), (5, 9), (5, 10), (5, 11), (6, 13), (6, 14), (6, 15), (6, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(7, 1), (7, 2), (7, 3), (7, 3), (7, 4), (7, 5), (7, 6), (7, 6), (7, 7), (7, 10), (7,11), (5, 12), (7, 13), (7, 14), (7, 15), (6, 16), (6, 17), (5, 18), (5, 19), (5, 19)], [(8, 1), (8, 2), (8, 3), (8, 4), (8, 4), (8, 5), (8, 6), (8, 7), (8, 7),(7, 10), (7, 11), (8, 12), (8, 13), (7, 14), (7, 15), (7, 16), (6, 17), (6, 18), (6, 19), (6, 19)], [(9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 5), (9, 6), (9, 7), (8, 8), (8, 9), (8, 10), (8, 11), (9, 12), (8, 14), (8, 15), (7, 16), (7, 17), (7, 18), (7, 19), (7, 19)], [(10, 0), (9, 2), (9, 3), (9, 4), (9, 5), (10, 6), (10, 6), (10, 7), (9, 8), (9, 9), (9, 10), (9, 11), (10, 12), (10, 13), (8, 15), (8, 16), (8, 17), (8, 18), (8, 19), (8, 19)], [(11, 0), (10, 1), (10, 2), (10, 3), (10, 4), (10, 5), (11, 6), (11, 7), (10, 8), (10, 9), (10, 10), (10, 11), (11, 12), (11, 13), (11, 14), (11, 15), (11, 17), (11, 18), (11, 19), (11, 19)], [(12, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5), (12, 7), (12, 7), (11,8), (11, 9), (11, 10), (11, 11), (12, 12), (12, 13), (12, 14), (12, 15), (12, 15), (12, 18), (12, 19), (12, 19)], [(13, 1), (12, 1), (12, 2), (12, 3), (12, 4), (12, 5), (13, 7), (13, 8), (12, 8), (12, 9), (12, 10), (12, 11), (13, 13), (13, 13), (13, 14), (13, 15), (13, 15), (12, 18), (13, 19), (13, 19)], [(14, 1), (14, 2), (14, 3), (14, 4), (14, 5), (14, 6), (14, 7), (14, 8), (14, 9), (14, 10), (14, 11), (14, 12), (14, 13), (14, 14), (14, 14), (14, 15), (13, 16), (13, 17), (13, 18), (14, 19)], [(15, 1), (15, 2), (15, 3), (15, 4), (15, 5), (15, 6), (15, 7), (15, 8), (15, 9), (15, 10), (15, 11), (15, 12), (15, 13), (15, 14), (15, 15), (15, 15), (14, 16), (14, 17), (14, 18), (15, 19)], [(16, 1), (16, 2), (16, 3), (16, 4), (16, 5), (16, 6), (16, 7), (16, 8), (16, 9), (16, 10), (16, 11), (16, 12), (16, 13), (16, 14), (16, 15), (16, 16), (15, 16), (15, 17), (15, 18), (16, 19)], [(17, 1), (17, 2), (17, 3), (17, 4), (17, 5), (17, 6), (17, 7), (17, 8), (17, 9), (17, 10), (17, 11), (17, 12), (17, 13), (17, 14), (17, 15), (17, 16), (17, 17), (17, 18), (17, 19), (17, 19)], [(18, 1),(18, 2), (18, 3), (18, 4), (18, 5), (18, 6), (18, 7), (18, 8), (18, 9), (18, 10), (18, 11), (18, 12), (18, 13), (18, 14), (18, 15), (18, 16), (18, 17), (18, 18), (18, 19), (18, 19)], [(19, 1), (19, 2), (19, 3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19,15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)], [(19, 1), (19, 2), (19,3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19, 15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)]]
    arrivees = {(20,99),(99,99)}
    """

    global r

    start = time()
    etapes = mouvement_foule(dim,departs,listes_murs,directions)
    print(time()-start)

    print(etapes)

    abscisses_traces_murs = []
    ordonnees_traces_murs = []
    i =0
    for liste_murs in listes_murs:
        if len(liste_murs) > 1:
            abscisses_traces_murs.append([])
            ordonnees_traces_murs.append([])
            for m in liste_murs:
                abscisses_traces_murs[i].append(m[0])
                ordonnees_traces_murs[i].append(m[1])
                plt.plot(abscisses_traces_murs[i],ordonnees_traces_murs[i],'b')
            i += 1
        else:
            plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')

    for a in arrivees:
        plt.plot(a[0],a[1],'gv')

    for i in range(len(departs)):
        abscisses_traj = []
        ordonnees_traj = []
        for t in range(len(etapes)):
            if etapes[t][i] == None:
                break
            else:
                abscisses_traj.append(etapes[t][i][0])
                ordonnees_traj.append(etapes[t][i][1])
        plt.plot(abscisses_traj,ordonnees_traj)

    plt.show()

def test_show_mouvement_foule(dim,departs,listes_murs,directions,arrivees,instants):
    """
    dim = 100
    departs = {(20,20),(2,80)}
    listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]] + [[(i,40) for i in range(10,29)] + [(29,j) for j in range(40,60)] + [(i,59) for i in range(29,9,-1)] + [(10,j) for j in range(59,39,-1)]] + [[(i,40) for i in range(40,65)] + [(64,j) for j in range(44,60)] + [(i,59) for i in range(64,39,-1)] + [(40,j) for j in range(59,39,-1)]] + [[(i,5) for i in range(50,65)] + [(64,j) for j in range(5,30)] + [(i,29) for i in range(64,49,-1)] + [(50,j) for j in range(29,4,-1)]] + [[(i,80) for i in range(65,80)] + [(79,j) for j in range(80,95)] + [(i,94) for i in range(79,64,-1)] + [(65,j) for j in range(94,79,-1)]]
    directions = [[(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 4), (0, 7), (0, 8), (0, 9), (0, 10), (1, 11), (1, 12), (1, 13), (1, 14), (1, 15), (1, 16), (1, 17), (1, 18), (1, 19), (1, 19)], [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 5), (1, 8), (1, 9), (1, 10), (1, 11), (2, 12), (2, 13), (2, 14), (2, 15), (2, 16), (2, 17), (2, 18), (2, 19), (2, 19)], [(3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 6), (2, 8), (2, 9), (2, 10), (2, 11), (3, 13), (3, 14), (3, 15), (3, 16), (3, 17), (3, 18), (3, 19), (3, 19)], [(4, 0), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 6), (3, 8), (3, 9), (3, 10), (3, 11), (4, 13), (4, 14), (4,15), (4, 16), (4, 17), (4, 18), (4, 19), (4, 19)], [(5, 1), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 6), (4, 8), (4, 9), (4, 10), (4, 11), (5, 13), (5, 14), (5, 15), (5, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(6, 1), (6, 2), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 6), (5, 8), (5, 9), (5, 10), (5, 11), (6, 13), (6, 14), (6, 15), (6, 16), (5, 17), (5, 18), (4, 19), (4, 19)], [(7, 1), (7, 2), (7, 3), (7, 3), (7, 4), (7, 5), (7, 6), (7, 6), (7, 7), (7, 10), (7,11), (5, 12), (7, 13), (7, 14), (7, 15), (6, 16), (6, 17), (5, 18), (5, 19), (5, 19)], [(8, 1), (8, 2), (8, 3), (8, 4), (8, 4), (8, 5), (8, 6), (8, 7), (8, 7),(7, 10), (7, 11), (8, 12), (8, 13), (7, 14), (7, 15), (7, 16), (6, 17), (6, 18), (6, 19), (6, 19)], [(9, 1), (9, 2), (9, 3), (9, 4), (9, 5), (9, 5), (9, 6), (9, 7), (8, 8), (8, 9), (8, 10), (8, 11), (9, 12), (8, 14), (8, 15), (7, 16), (7, 17), (7, 18), (7, 19), (7, 19)], [(10, 0), (9, 2), (9, 3), (9, 4), (9, 5), (10, 6), (10, 6), (10, 7), (9, 8), (9, 9), (9, 10), (9, 11), (10, 12), (10, 13), (8, 15), (8, 16), (8, 17), (8, 18), (8, 19), (8, 19)], [(11, 0), (10, 1), (10, 2), (10, 3), (10, 4), (10, 5), (11, 6), (11, 7), (10, 8), (10, 9), (10, 10), (10, 11), (11, 12), (11, 13), (11, 14), (11, 15), (11, 17), (11, 18), (11, 19), (11, 19)], [(12, 0), (11, 1), (11, 2), (11, 3), (11, 4), (11, 5), (12, 7), (12, 7), (11,8), (11, 9), (11, 10), (11, 11), (12, 12), (12, 13), (12, 14), (12, 15), (12, 15), (12, 18), (12, 19), (12, 19)], [(13, 1), (12, 1), (12, 2), (12, 3), (12, 4), (12, 5), (13, 7), (13, 8), (12, 8), (12, 9), (12, 10), (12, 11), (13, 13), (13, 13), (13, 14), (13, 15), (13, 15), (12, 18), (13, 19), (13, 19)], [(14, 1), (14, 2), (14, 3), (14, 4), (14, 5), (14, 6), (14, 7), (14, 8), (14, 9), (14, 10), (14, 11), (14, 12), (14, 13), (14, 14), (14, 14), (14, 15), (13, 16), (13, 17), (13, 18), (14, 19)], [(15, 1), (15, 2), (15, 3), (15, 4), (15, 5), (15, 6), (15, 7), (15, 8), (15, 9), (15, 10), (15, 11), (15, 12), (15, 13), (15, 14), (15, 15), (15, 15), (14, 16), (14, 17), (14, 18), (15, 19)], [(16, 1), (16, 2), (16, 3), (16, 4), (16, 5), (16, 6), (16, 7), (16, 8), (16, 9), (16, 10), (16, 11), (16, 12), (16, 13), (16, 14), (16, 15), (16, 16), (15, 16), (15, 17), (15, 18), (16, 19)], [(17, 1), (17, 2), (17, 3), (17, 4), (17, 5), (17, 6), (17, 7), (17, 8), (17, 9), (17, 10), (17, 11), (17, 12), (17, 13), (17, 14), (17, 15), (17, 16), (17, 17), (17, 18), (17, 19), (17, 19)], [(18, 1),(18, 2), (18, 3), (18, 4), (18, 5), (18, 6), (18, 7), (18, 8), (18, 9), (18, 10), (18, 11), (18, 12), (18, 13), (18, 14), (18, 15), (18, 16), (18, 17), (18, 18), (18, 19), (18, 19)], [(19, 1), (19, 2), (19, 3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19,15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)], [(19, 1), (19, 2), (19,3), (19, 4), (19, 5), (19, 6), (19, 7), (19, 8), (19, 9), (19, 10), (19, 11), (19, 12), (19, 13), (19, 14), (19, 15), (19, 16), (19, 17), (19, 18), (19, 19), (19, 19)]]
    arrivees = {(20,99),(99,99)}
    """

    global r

    start = time()
    etapes = mouvement_foule(dim,departs,listes_murs,directions)
    print(time()-start)

    print(etapes)

    for t in instants:
        if t < len(etapes):
            print(t)

            figure,axes = plt.subplots()
            axes.set_aspect(1)
            for pt in (etapes[t]):
                if pt != None:
                    circle = plt.Circle(pt,r,facecolor = 'orange',edgecolor = 'k')
                    axes.add_artist(circle)

            for i in range(len(departs)):
                abscisses_traj = []
                ordonnees_traj = []
                for t in range(len(etapes)):
                    if etapes[t][i] == None:
                        break
                    else:
                        abscisses_traj.append(etapes[t][i][0])
                        ordonnees_traj.append(etapes[t][i][1])
                plt.plot(abscisses_traj,ordonnees_traj)

            abscisses_traces_murs = []
            ordonnees_traces_murs = []
            i =0
            for liste_murs in listes_murs:
                if len(liste_murs) > 1:
                    abscisses_traces_murs.append([])
                    ordonnees_traces_murs.append([])
                    for m in liste_murs:
                        abscisses_traces_murs[i].append(m[0])
                        ordonnees_traces_murs[i].append(m[1])
                        plt.plot(abscisses_traces_murs[i],ordonnees_traces_murs[i],'b')
                    i += 1
                else:
                    plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')

            for a in arrivees:
                plt.plot(a[0],a[1],'gv')

            plt.show()

def show_etapes(dim,listes_murs,arrivees,etapes,instants):
    for t in instants:
        if t < len(etapes):
            print(t)

            figure,axes = plt.subplots()
            axes.set_aspect(1)
            for pt in (etapes[t]):
                if pt != None:
                    circle = plt.Circle(pt,r,facecolor = 'orange',edgecolor = 'k')
                    axes.add_artist(circle)

            for i in range(len(etapes[0])):
                abscisses_traj = []
                ordonnees_traj = []
                for t in range(len(etapes)):
                    if etapes[t][i] == None:
                        break
                    else:
                        abscisses_traj.append(etapes[t][i][0])
                        ordonnees_traj.append(etapes[t][i][1])
                plt.plot(abscisses_traj,ordonnees_traj)

            abscisses_traces_murs = []
            ordonnees_traces_murs = []
            i =0
            for liste_murs in listes_murs:
                if len(liste_murs) > 1:
                    abscisses_traces_murs.append([])
                    ordonnees_traces_murs.append([])
                    for m in liste_murs:
                        abscisses_traces_murs[i].append(m[0])
                        ordonnees_traces_murs[i].append(m[1])
                        plt.plot(abscisses_traces_murs[i],ordonnees_traces_murs[i],'b')
                    i += 1
                else:
                    plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')

            for a in arrivees:
                plt.plot(a[0],a[1],'gv')

            plt.show()
