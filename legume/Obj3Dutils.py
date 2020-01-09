### set de fonctions utiles pour manipuler les geometry plantGL

from openalea.plantgl.all import *
from V3Dutils import *
from scipy import pi, array, sin , cos

def mesh(geometry):
    """ renvoie le mesh d'une geometry"""
    #d = Discretizer()
    #geometry.apply(d)
    #return d.result
    tessel = Tesselator()
    geometry.apply(tessel)
    mesh_ = tessel.triangulation
    return mesh_


def tri(p1, p2, p3):
    """ renvoie le TriangleSet d'un triangle a partir de ses 3 points """
    points= Point3Array([ Vector3(p1[0],p1[1],p1[2]),Vector3(p2[0],p2[1],p2[2]),Vector3(p3[0],p3[1],p3[2])]) 
    indices= Index3Array([ Index3(0,1,2)])
    return TriangleSet(points, indices)

def quadform(p1, p2, p3, p4, opt=None):
    """ renvoie le TriangleSet a 4 points ; si opt != None inverse les points pour representation des triangles"""
    points= Point3Array([ Vector3(p1[0],p1[1],p1[2]),Vector3(p2[0],p2[1],p2[2]),Vector3(p3[0],p3[1],p3[2]), Vector3(p4[0],p4[1],p4[2])]) 
    if opt == None:
        indices= Index3Array([ Index3(0,1,3), Index3(1,2,3)])
    else:
        indices= Index3Array([ Index3(0,1,2), Index3(0,2,3)])

    return TriangleSet(points, indices)


def turtle36():
    ## demi triangles hauts
    ## surface du turtle36 de 1 de rayon 2.9577
    #liste points
    ls_pt = []
    dazi = pi/6.
    for elv in [0, pi/6., pi/3.]:
        r = cos(elv)
        z = sin(elv)
        for i in range(12):
            teta = i*dazi
            x = r*cos(teta)
            y = r*sin(teta)
            ls_pt.append(Vector3(x,y,z))

    ls_pt.append(Vector3(0,0,1.))

    #liste d'index
    ls_id = []
    #1er ring
    for i in range(11):
        ls_id.append(Index3(i, i+12, i+13))

    ls_id.append(Index3(11, 11+12, 12))
    #2e ring
    for i in range(11):
        ls_id.append(Index3(i+12, i+24, i+25))

    ls_id.append(Index3(23, 23+12, 24))
    #3e ring
    for i in range(11):
        ls_id.append(Index3(i+24, i+25, 36))

    ls_id.append(Index3(35, 24, 36))
    
    # triangleset
    points= Point3Array(ls_pt) 
    indices= Index3Array(ls_id)
    return  TriangleSet(points, indices) #36 tiangles en 3 rings


def transformation(obj, sx, sy, sz, rx, ry, rz, tx, ty, tz ): 
    """ Return a scaled, rotated and translated 3D object - Similar to 'transformation' in PovRay """ 
    s_obj = Scaled (Vector3(sx,sy,sz), obj)
    r_obj = EulerRotated (rx, ry, rz, s_obj)
    t_obj = Translated (Vector3(tx,ty,tz), r_obj)
    return t_obj


def conv_cyl(p1, p2, r):
    """ calcule des parametres pour le positionnement d'un cylindre a partir des coordonnees 
    des deux points extremes et du rayon """
    vec = XyzToPol (p2-p1)
    return p1, vec[0], r, vec[1], vec[2] #p1, longueur, rayon, azi, incli

def euler_normal(AA, BB, CC):
    """ compute the normal at the plane defined by euler angles AA, BB, CC """
    tri0 = tri(scipy.array([0.,0.,0.]),scipy.array([1.,0.,0.]),scipy.array([0.,1.,0.]))
    v = mesh(transformation(tri0, 1,1,1,AA+pi/2.,CC,-BB,0,0,0))
    v.computeNormalList()
    return v.normalAt(0)

def tri_ortho(p1,p2,p3):
    """ calcul  orthocentre d'un triangle - meme resulat que .faceCenter(id) sur un triangle set"""
    #p4 milieu  p2-p3
    p4 = array(p2)+0.5*(array(p3)-array(p2))
    return array(p1)+2.*(p4-array(p1))/3.

def compute_ortho_list(ind_ls, pt_ls, epsilon = 0.001):
    """ calcule liste des orthocentre d'un triangle set """
    surf = compute_surface_list(ind_ls, pt_ls)
    ortho_ls = []
    for i in range(len(ind_ls)):
        if surf[i]>epsilon:
            p1, p2, p3 = pt_ls[ind_ls[i][0]], pt_ls[ind_ls[i][1]], pt_ls[ind_ls[i][2]] 
            ortho = tri_ortho(p1,p2,p3)
            ortho_ls.append(Vector3(ortho[0], ortho[1], ortho[2]))

    return ortho_ls

def compute_normal_list(ind_ls, pt_ls, epsilon = 0.001):
    """ calcule liste des normales d'un triangle set - contourne pb des nb normale different des nb faces avec .computeNormalList()"""
    surf = compute_surface_list(ind_ls, pt_ls)
    n_ls = []
    for i in range(len(ind_ls)):
        if surf[i]>epsilon:
            p1, p2, p3 = array(pt_ls[ind_ls[i][0]]), array(pt_ls[ind_ls[i][1]]), array(pt_ls[ind_ls[i][2]]) 
            pv = produit_vectoriel ((p2-p1), (p3-p1))
            n_ls.append(pv/norme_v(pv))

    return n_ls

def mesh_points(geometry):
    """ get the mesh points of a geometry """
    g = mesh(geometry)
    return map(array, g.pointList)

def triangle_area(p1, p2, p3):
    """ compute surface area of a triangle """
    u = p2-p1
    v = p3-p1
    return 0.5*norme_v(produit_vectoriel (u, v))

def compute_surface_list(ind_ls, pt_ls):
    """ calcule liste des surfaces d'un triangle set """
    s_ls = []
    for i in range(len(ind_ls)):
        p1, p2, p3 = array(pt_ls[ind_ls[i][0]]), array(pt_ls[ind_ls[i][1]]), array(pt_ls[ind_ls[i][2]]) 
        s_ls.append(triangle_area(p1, p2, p3))

    return s_ls


def leg_leaf(Lmax, largmax, alpha=0., gamma=0., unifol=0):
    gamma = gamma * 3.14 / 180  # en radians
    lf, la, pe, br, crois = 21. / 21., 6.5 / 21., 10. / 21., 3.6 / 21., 0.1 / 21.  # leaf Trudeau modifie
    leaf = quadform(array([0., 0., 0.]), array([0.5, 0.5, 0.]), array([0., 1., 0.]), array([-0.5, 0.5, 0.]),
                    opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    up = transformation(leaf, 1, 1, 1, 0, 0, gamma, 0, pe / lf * Lmax, 0)
    right = transformation(leaf, 1, 1, 1, -3.14 / 180 * 80, 0, gamma, br / 2. * Lmax, crois * Lmax, 0)
    left = transformation(leaf, 1, 1, 1, 3.14 / 180 * 80, 0, gamma, -br / 2. * Lmax, crois * Lmax, 0)
    if unifol == 0:
        return Group([up, right, left])  # groupe les differents geom
    else:
        return up  # 1 foliole pour la premiere feuille


def leg_leaf_lucas(Lmax, largmax, alpha=0., gamma=0., nfol=3, angfol=10., ecfol=6., anginit=45.,
                   geom=True):  # angfol : pi/nombre de rangs n?ssaires pour boucler le demi-cercle // ecfol : longueur de rachis entre chaque paire de folioles (mm).
    nr = (nfol - 1) / 2 if nfol % 2 == 1 else nfol / 2  # nombre de paires de folioles lateraux
    angfol = 10  # 180/nr #demi-cercle complet obtenu sur le nombre de rangs de paires de folioles. Attention, ne marche que si l'ecartement ecfol est constant!
    gamma = gamma * 3.14 / 180  # en radians
    anginit = anginit * 3.14 / 180
    angfol = angfol * 3.14 / 180
    ecfol = ((
                         1.88 * nfol ** -0.54) * Lmax) * 10  # relation empirique entre nombre de folioles et ecartement entre deux rang? pour le sainfoin.

    ls_pts = []

    lf, la, pe, br, crois = 21. / 21., 6.5 / 21., 10. / 21., 3.6 / 21., 0.1 / 21.  # leaf Trudeau modifie
    leaf = quadform(array([0., 0., 0.]), array([0.5, 0.5, 0.]), array([0., 1., 0.]), array([-0.5, 0.5, 0.]),
                    opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    angup = (angfol * -(nr - 1)) - 3.14 - anginit  # angle de placement du foliole central, au bout de la chaine
    up = transformation(leaf, 1, 1, 1, 0, 0, gamma, 0, (pe / lf * Lmax) + (ecfol * (cos(angup) + cos(anginit))),
                        ecfol * (sin(angup) - sin(anginit)))
    ls_pts.append(
        array([0, (pe / lf * Lmax) + (ecfol * (cos(angup) + cos(anginit))), ecfol * (sin(angup) - sin(anginit))]))
    # right = transformation(leaf, 1,1,1,-3.14/180*80,0,gamma, br/2.*Lmax, crois*Lmax,0)
    # left = transformation(leaf, 1,1,1,3.14/180*80,0,gamma, -br/2.*Lmax, crois*Lmax,0)
    listfol = [up] if nfol % 2 == 1 else []
    for i in range(nr):  # nombre de paires de folioles lateraux
        ang = (angfol * -i) - 3.14 - anginit
        ecfolopp = (sin(ang) - sin(anginit)) * ecfol
        ecfoladj = (cos(ang) + cos(anginit)) * ecfol
        listfol.append(transformation(leaf, 1, 1, 1, -3.14 / 180 * 80, 0, gamma, br / 2. * Lmax, ecfoladj, ecfolopp))
        listfol.append(transformation(leaf, 1, 1, 1, 3.14 / 180 * 80, 0, gamma, -br / 2. * Lmax, ecfoladj, ecfolopp))
        ls_pts.append(array([br / 2. * Lmax, ecfoladj, ecfolopp]))
        ls_pts.append(array([-br / 2. * Lmax, ecfoladj, ecfolopp]))

    if geom == True:
        return Group(listfol)  # groupe les differents geom
    else:
        return ls_pts


def geomstip(Lmax, largmax, alpha=0., gamma=0.):
    gamma = gamma * 3.14 / 180  # en radians
    stip = quadform(array([0., 0., 0.]), array([0.5, 0.5, 0.]), array([0., 1., 0.]), array([-0.5, 0.5, 0.]),
                    opt=2)  # prends pas alpha en compte
    stip = transformation(stip, largmax, Lmax, 1., 0, 0, 0, 0, 0, 0)
    if Lmax >= largmax:
        right = transformation(stip, 1, 1, 1, -3.14 / 180 * alpha, -gamma, 0, 0, 0, 0)
        left = transformation(stip, 1, 1, 1, 3.14 / 180 * alpha, gamma, 0, 0, 0, 0)
    else:
        right = transformation(stip, 1, 1, 1, -3.14 / 180 * alpha, 0, gamma, 0, 0, 0)
        left = transformation(stip, 1, 1, 1, 3.14 / 180 * alpha, 0, gamma, 0, 0, 0)
    return Group([right, left])
    # par convention choisi longueur dans direction du petiole?
    # faire porter gamma sur la plus grande direction -> marche normalement pour pois
    # reprendre gammaFeuil pour le gamma (pas IncPet petiole)


def leg_grass(Lmax, largmax, gamma=0., angfol=10., nfol=8, anginit=45., geom=True):
    anginit = anginit * 3.14 / 180
    angfol = angfol * 3.14 / 180
    ecfol = Lmax / nfol  # longueur de segment de feuille

    ls_pts = []

    leaf = quadform(array([-0.5, 0., 0.]), array([-0.5, 1., 0.]), array([0.5, 1., 0.]), array([0.5, 0., 0.]),
                    opt=2)  # prends pas alpha en compte
    leaf = transformation(leaf, largmax, ecfol, 1., 0, 0, 0, 0, 0, 0)
    bottom = transformation(leaf, 1, 1, 1, 0, 0, 3.14 / 2 - anginit, 0, 0, 0)

    # ajout points 2 premiers points debuts + 1er segments
    ls_pts.append(array([0, 0, 0]))
    ls_pts.append(array([0, sin(anginit) * ecfol, cos(anginit) * ecfol]))
    listfol = [bottom]

    for i in range(1, nfol):  # nombre de segment restants
        ang = anginit - (angfol * -i)  # (angfol * -i) - 3.14 - anginit
        z = ls_pts[-1][2] + cos(ang) * ecfol
        y = ls_pts[-1][1] + sin(ang) * ecfol
        listfol.append(transformation(leaf, 1, 1, 1, 0, 0, 3.14 / 2 - ang, 0, ls_pts[-1][1], ls_pts[-1][2]))
        ls_pts.append(array([0, y, z]))
        #print i, ecfol, angfol, ang, distance(ls_pts[-1], ls_pts[-2])

    if geom == True:
        return Group(listfol)  # groupe les differents geom
    else:
        return ls_pts[1:]  # retire premier point -> surface attribuee aux extremite de segment