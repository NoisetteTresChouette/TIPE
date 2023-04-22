import matplotlib.pyplot as plt
import numpy as np
from time import time

"""
dim = 20
repulsifs = {'inv': {(1,j) for j in range(15,20)}, 'plateau':{(i,j) for i in range(6,8) for j in range(7,13)} }
listes_murs = [[(i,-1) for i in range(-1,21)] + [(20,j) for j in range(21)] + [(i,20) for i in range(19,-2,-1)] + [(-1,j) for j in range(19,-2,-1)]] + [[(i,8) for i in range(2,6)] + [(5,j) for j in range(8,12)] + [(i,11) for i in range(5,1,-1)] + [(2,j) for j in range(11,7,-1)]] +[[(i,j)] for i in range(3,5) for j in range(9,11)]+ [[(i,8) for i in range(8,13)] + [(12,j) for j in range(8,12)] + [(i,11) for i in range(12,7,-1)] + [(8,j) for j in range(11,7,-1)]] + [[(i,j)] for i in range(9,12) for j in range(9,11)]+[[(i,1) for i in range(10,13)] + [(12,j) for j in range(1,6)] + [(i,5) for i in range(12,9,-1)] + [(10,j) for j in range(5,0,-1)]]+ [[(11,j)] for j in range(2,5)]+ [[(i,16) for i in range(13,16)] + [(15,j) for j in range(16,19)] + [(i,18) for i in range(15,12,-1)] + [(13,j) for j in range(18,15,-1)]]+[[(14,17)]]
arrivees = {(4,19),(19,19)}

dim = 100
repulsifs = {'inv': {(i,j) for j in range(75,100) for i in range(5,10)}, 'plateau':{(i,j) for i in range(30,40) for j in range(35,65)} }
listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]] + [[(i,40) for i in range(10,29)] + [(29,j) for j in range(40,60)] + [(i,59) for i in range(29,9,-1)] + [(10,j) for j in range(59,39,-1)]] +[[(i,j)] for i in range(11,29) for j in range(41,59)]+ [[(i,40) for i in range(40,65)] + [(64,j) for j in range(44,60)] + [(i,59) for i in range(64,39,-1)] + [(40,j) for j in range(59,39,-1)]] + [[(i,j)] for i in range(41,64) for j in range(41,59)]+[[(i,5) for i in range(50,65)] + [(64,j) for j in range(5,30)] + [(i,29) for i in range(64,49,-1)] + [(50,j) for j in range(29,4,-1)]]+ [[(i,j)] for j in range(6,29) for i in range(51,64)]+ [[(i,80) for i in range(65,80)] + [(79,j) for j in range(80,95)] + [(i,94) for i in range(79,64,-1)] + [(65,j) for j in range(94,79,-1)]]+[[(i,j)] for i in range(66,79) for j in range(81,94)]
arrivees = {(20,99),(99,99)}
"""

def distance(A,B):
    return ((A[0]-B[0])**2 + (A[1] - B[1])**2)**(1/2)

importance_potent = 0.75
importance_deplacement = 0.25/(2**(1/2))

def potent(pt,repulsifs,murs):
    if pt in murs:
        return np.inf
    V = 0
    for r in repulsifs['inv']:
        V += 1/(1+distance(r,pt))
    if pt in repulsifs['plateau']:
        V += 1
    return V

def calcul_poids(dim,repulsifs,murs):
    poids = [[[[importance_potent*potent((i+l-1,j+c-1),repulsifs,murs) + importance_deplacement*distance((i,j),(i+l-1,j+c-1)) for c in range(3)] for l in range(3)] for j in range(dim)] for i in range(dim)]
    return poids

def dijkstra(dim,depart,arrivees,poids):#"arrivees" est un ENSEMBLE de points
    poids_tot = [[np.inf for j in range(dim+1)] for i in range(dim+1)]# dim+1 pour inclure les contours de l'expérience et ne pas avoir un index out of range dans le while
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
            """plt.plot(liste_murs[0][0],liste_murs[0][1],'b+')"""
            murs.add(liste_murs[0])

    for r in repulsifs['inv']:
        plt.plot(r[0],r[1],'rv')

    for r in repulsifs['plateau']:
        plt.plot(r[0],r[1],'ro')

    map_potent = np.array([[potent((i,j),repulsifs,murs) for i in range(dim)] for j in range(dim)])
    plt.imshow(map_potent,cmap='hot')

    plt.plot(depart[0],depart[1],'go')

    for a in arrivees:
        plt.plot(a[0],a[1],'gv')

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

def directions_souhaitees(dim,arrivees,repulsifs,murs):
    poids = calcul_poids(dim,repulsifs,murs)
    directions = [[(i,j) for j in range(dim)] for i in range(dim)]
    points_fleches = arrivees.copy()
    classement = classement_distance(dim,arrivees,murs)
    for depart in classement:
        if depart not in points_fleches:
            chemin = dijkstra(dim,depart,arrivees,poids)
            for i in range(len(chemin)-1):
                if chemin[i] in points_fleches:
                    break
                else:
                    points_fleches.add(chemin[i])

                    print(len(points_fleches))

                    directions[chemin[i][0]][chemin[i][1]] = (chemin[i+1][0],chemin[i+1][1])
    return directions

def classement_distance(dim,arrivees,murs):
    map = []
    for i in range(dim):
        for j in range(dim):
            if (i,j) not in murs:
                map.append((i,j))
    map = tri_fusion(map,arrivees)
    return map

def tri_fusion(t,arrivees):
    if len(t) > 1:
        t1 = tri_fusion(t[:len(t)//2],arrivees)
        t2 = tri_fusion(t[len(t)//2:],arrivees)
        t = fusion(t1,t2,arrivees)
        return t
    else:
        return t

def fusion(t1,t2,arrivees):
    t_retour = []
    i=0
    j=0
    l1 = len(t1)
    l2 = len(t2)
    while i<l1 or j<l2:
        if j == l2 or (i<l1 and min([distance(t1[i],arrivee) for arrivee in arrivees]) > min([distance(t2[j],arrivee) for arrivee in arrivees])):
            t_retour.append(t1[i])
            i+=1
        else:
            t_retour.append(t2[j])
            j+=1
    return t_retour

def test_show_directions_souhaitees(dim,arrivees,repulsifs,listes_murs):
    """
    dim = 3
    repulsifs = {(1,1)}
    listes_murs = [[(i,-1) for i in range(-1,4)] + [(3,j) for j in range(4)] + [(i,3) for i in range(2,-2,-1)] + [(-1,j) for j in range(2,-2,-1)]]
    arrivees = {(2,2)}

    dim = 100
    repulsifs = set([(i,50) for i in range(40)] + [(i,50) for i in range(60,100)])
    listes_murs = [[(i,-1) for i in range(-1,101)] + [(100,j) for j in range(101)] + [(i,100) for i in range(99,-2,-1)] + [(-1,j) for j in range(99,-2,-1)]]
    arrivees = {(99,99)}
    """

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
            murs.add(liste_murs[0])

    for r in repulsifs['inv']:
        plt.plot(r[0],r[1],'rv')
    for r in repulsifs['plateau']:
        plt.plot(r[0],r[1],'ro')

    for a in arrivees:
        plt.plot(a[0],a[1],'gv')

    start = time()
    directions = directions_souhaitees(dim,arrivees,repulsifs,murs)
    print(time()-start)
    print(directions)

    abscisses_pt = []
    ordonnees_pt = []
    abscisses_vect = []
    ordonnees_vect = []
    for i in range(dim):
        for j in range(dim):
            abscisses_pt.append(i)
            ordonnees_pt.append(j)
            abscisses_vect.append(directions[i][j][0] - i)
            ordonnees_vect.append(directions[i][j][1] - j)

    plt.quiver(abscisses_pt,ordonnees_pt,abscisses_vect,ordonnees_vect,units='xy',scale=1.5)

    plt.show()

def set_murs(liste_murs):
    murs = set()
    for mur in liste_murs:
        for m in mur:
            murs.add(m)
    return murs

def show_flechage_inutil():
    murs = set()
    abscisses_traces_murs = []
    ordonnees_traces_murs = []
    listes_murs = [[(i,-1) for i in range(-1,11)] + [(10,j) for j in range(11)] + [(i,10) for i in range(9,-2,-1)] + [(-1,j) for j in range(9,-2,-1)]]

    i =0
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
            murs.add(liste_murs[0])
    arrivees = {(9,5)}
    repulsifs = {'inv':set(),'plateau':{}}
    dim = 10
    poids = calcul_poids(dim,repulsifs,murs)

    plt.plot(9,5,'gv')
    plt.text(8.5,4.3,'arrivée',color='g')

    chemin1 = dijkstra(dim,(5,9),arrivees,poids)
    chemin2 = dijkstra(dim,(0,9),arrivees,poids)

    plt.plot(5,9,'ko')
    plt.text(4.5,9.4,'départ 1',color = 'k')

    plt.plot(0,9,'bo')
    plt.text(-0.4,8.4,'départ 2',color = 'b')

    abscisses_pt1 = []
    ordonnees_pt1 = []
    abscisses_vect1 = []
    ordonnees_vect1 = []
    for i in range(len(chemin1)-1):
        abscisses_pt1.append(chemin1[i][0])
        ordonnees_pt1.append(chemin1[i][1])
        abscisses_vect1.append(chemin1[i+1][0] - chemin1[i][0])
        ordonnees_vect1.append(chemin1[i+1][1] - chemin1[i][1])

    abscisses_pt2 = []
    ordonnees_pt2 = []
    abscisses_vect2 = []
    ordonnees_vect2 = []
    for i in range(len(chemin2)-1):
        if chemin2[i] in chemin1:
            break
        else:
            abscisses_pt2.append(chemin2[i][0])
            ordonnees_pt2.append(chemin2[i][1])
            abscisses_vect2.append(chemin2[i+1][0] - chemin2[i][0])
            ordonnees_vect2.append(chemin2[i+1][1] - chemin2[i][1])

    plt.quiver(abscisses_pt1,ordonnees_pt1,abscisses_vect1,ordonnees_vect1,units='xy',scale=1.5)
    plt.quiver(abscisses_pt2,ordonnees_pt2,abscisses_vect2,ordonnees_vect2,units='xy',scale=1.5,color = 'b')

    plt.show()

def show_flechage_distance():
    murs = set()
    abscisses_traces_murs = []
    ordonnees_traces_murs = []
    listes_murs = [[(i,-1) for i in range(-1,11)] + [(10,j) for j in range(11)] + [(i,10) for i in range(9,-2,-1)] + [(-1,j) for j in range(9,-2,-1)]]

    i =0
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
            murs.add(liste_murs[0])
    arrivees = {(9,5)}
    repulsifs = {'inv':set(),'plateau':{}}
    dim = 10
    poids = calcul_poids(dim,repulsifs,murs)

    plt.plot(9,5,'gv')

    chemins1 = [dijkstra(dim,(5,9),arrivees,poids),dijkstra(dim,(5,5),arrivees,poids),dijkstra(dim,(5,1),arrivees,poids)]
    chemins2 = [dijkstra(dim,(0,9),arrivees,poids), dijkstra(dim,(0,5),arrivees,poids),dijkstra(dim,(0,1),arrivees,poids)]

    plt.plot(5,9,'ko')
    plt.plot(5,5,'ko')
    plt.plot(5,1,'ko')

    plt.plot(0,9,'bo')
    plt.plot(0,5,'bo')
    plt.plot(0,1,'bo')

    abscisses_pt1 = []
    ordonnees_pt1 = []
    abscisses_vect1 = []
    ordonnees_vect1 = []
    for chemin1 in chemins1:
        for i in range(len(chemin1)-1):
            abscisses_pt1.append(chemin1[i][0])
            ordonnees_pt1.append(chemin1[i][1])
            abscisses_vect1.append(chemin1[i+1][0] - chemin1[i][0])
            ordonnees_vect1.append(chemin1[i+1][1] - chemin1[i][1])

    abscisses_pt2 = []
    ordonnees_pt2 = []
    abscisses_vect2 = []
    ordonnees_vect2 = []
    for chemin2 in chemins2:
        for i in range(len(chemin2)-1):
            abscisses_pt2.append(chemin2[i][0])
            ordonnees_pt2.append(chemin2[i][1])
            abscisses_vect2.append(chemin2[i+1][0] - chemin2[i][0])
            ordonnees_vect2.append(chemin2[i+1][1] - chemin2[i][1])

    plt.quiver(abscisses_pt2,ordonnees_pt2,abscisses_vect2,ordonnees_vect2,units='xy',scale=1.5,color = 'b')
    plt.quiver(abscisses_pt1,ordonnees_pt1,abscisses_vect1,ordonnees_vect1,units='xy',scale=1.5)

    plt.show()

def show_map(dim,arrivees,repulsifs,murs):
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
    for r in repulsifs['inv']:
        plt.plot(r[0],r[1],'rv')
    for r in repulsifs['plateau']:
        plt.plot(r[0],r[1],'ro')

    for a in arrivees:
        plt.plot(a[0],a[1],'gv')
    map_potent = np.array([[potent((i,j),repulsifs,murs) for i in range(dim)] for j in range(dim)])
    plt.imshow(map_potent,cmap='hot')

    plt.show()

def show_directions(directions,repulsifs,arrivees,listes_murs):

    for r in repulsifs['inv']:
        plt.plot(r[0],r[1],'rv')
    for r in repulsifs['plateau']:
        plt.plot(r[0],r[1],'ro')

    for a in arrivees:
        plt.plot(a[0],a[1],'gv')

    abscisses_pt = []
    ordonnees_pt = []
    abscisses_vect = []
    ordonnees_vect = []
    for i in range(dim):
        for j in range(dim):
            abscisses_pt.append(i)
            ordonnees_pt.append(j)
            abscisses_vect.append(directions[i][j][0] - i)
            ordonnees_vect.append(directions[i][j][1] - j)

    plt.quiver(abscisses_pt,ordonnees_pt,abscisses_vect,ordonnees_vect,units='xy',scale=1.5)

    plt.show()