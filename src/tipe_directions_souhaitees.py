from tipe_dijkstra import *

def directions_souhaitees(dim,arrivees,repulsifs,murs):

    start = time()

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

                    directions[chemin[i][0]][chemin[i][1]] = (chemin[i+1][0]/distance(chemin[i],chemin[i+1]),chemin[i+1][1]/distance(chemin[i],chemin[i+1]))

    print(time()-start)

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

def test_show_classement(dim,arrivees,listes_murs):

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

    classement = classement_distance(dim,arrivees,murs)
    for i in range(len(classement)):
        plt.text(classement[i][0],classement[i][1],str(i))
    plt.show()

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

    for r in repulsifs:
        plt.plot(r[0],r[1],'rx')

    directions = directions_souhaitees(dim,arrivees,repulsifs,murs)

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

    print(abscisses_pt)
    print(ordonnees_pt)
    print(abscisses_vect)
    print(ordonnees_vect)

    plt.quiver(abscisses_pt,ordonnees_pt,abscisses_vect,ordonnees_vect,units='xy',scale=1.5)

    plt.grid()
    plt.show()