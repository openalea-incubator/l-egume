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
