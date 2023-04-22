def poids_diff(pt1,pt2,repulsifs,murs,potent):
    return potent(pt2,repulsifs,murs)-potent(pt1,repulsifs,murs)

def poids_simple(pt1,pt2,repulsifs,murs,potent):
    return potent(pt2,repulsifs,murs)

