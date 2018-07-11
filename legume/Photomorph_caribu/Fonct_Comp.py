from numpy import *
from openalea.lpy import *
from openalea.plantgl.all import *
from generateScene import run

def Init_Lpy(lpy_filename, fichier_pois):

    lsys = Lsystem(lpy_filename)

    # Definir carto et Nombre plantes

    nbplantes = len(fichier_pois)  # Doit etre moins que 64 plantes !!!!!
    lsys.nbcote = int(sqrt(nbplantes))

    # Pattern & sol
    lsys.cote = 55
    lsys.nbplantes = nbplantes


    a = AxialTree()
    a.append(lsys.attente(1))

    lsys.cartographie = []

    for i in range(0, nbplantes):
        x = fichier_pois[i][0]
        y = fichier_pois[i][1]

        lsys.cartographie.append(array([x, y, 0.]))
        a.append(lsys.Sd(i))

    lsys.axiom = a
    axiom = lsys.axiom
    axiom = lsys.derive(axiom, 25)

    runL = run(lsys, axiom=axiom, nbstep=1)


    return runL,lsys, axiom, nbplantes






def PrepareScene(scene,runL):

    sc = runL[0]
    s_leg2 = Scene()
    s_leg = scene

    xy = runL[1].cote

    #points = [(0, 0, 0),(xy,0,0), (xy,xy,0),(0,xy,0)]
    points = [(0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]
    sol = QuadSet(points, indices, normals, indices)

    Dico_conv = {}

    nb_tri = 0

    # Pour le sol
    esp_opt = 000000000000
    n_plt = 0
    tp_org = 000000
    count_sh = esp_opt + n_plt + tp_org + nb_tri
    nb_tri += 1
    count_ID = 0


#Commente pour supprimer sol
    s_leg2.add(Shape(geometry=sol, id=count_ID, appearance=s_leg[0].appearance))
    s_leg2[count_ID].setName(str(count_sh))
    Dico_conv[count_ID] = 0
    count_ID += 1

    nb_tri = 0
    m = 0
    nump = 0

    Dico_Pet = {}
    Dico_In = {}
    Dico_Stp = {}
    Dico_Apex_ID = {}

    for x in s_leg:
        if x.id != s_leg[0].id and x.id != s_leg[1].id:

            nump = runL[0][x.id][0]

            if (sc[x.id].name == 'Pet'):
                liste = []
                liste.append(nump)
                liste.append(x.id)

                xa = runL[0][x.id + 1][0]
                ya = runL[0][x.id + 1][1]
                za = runL[0][x.id + 1][2]

                liste.append(xa)
                liste.append(ya)
                liste.append(za)

                if (xa < 0 or ya < 0 or x > 55 or ya > 55):
                    liste.append(-1)
                else:
                    liste.append(0)

                Dico_Pet[count_ID] = liste

            if (sc[x.id].name == 'In'):
                liste = []
                liste.append(nump)
                liste.append(x.id)

                liste.append(runL[0][x.id - 1][0])
                liste.append(runL[0][x.id - 1][1])
                liste.append(runL[0][x.id - 1][2])

                Dico_In[count_ID] = liste

            if (sc[x.id].name == 'A'):
                name_A = str(sc[x.id][0])+'_'+str(sc[x.id][1])+str(sc[x.id][5])
                Dico_Apex_ID[name_A] = x.id

            if (sc[x.id].name != 'solxy' and sc[x.id].name != 'p' and sc[x.id].name != 'RLAP' and sc[
                x.id].name != 'RS' and sc[x.id].name != 'RLB' and sc[x.id].name != 'RA' ):

                if (sc[x.id].name == 'A' or sc[x.id].name == 'Stp' or sc[x.id].name == 'Pet'):  # Feuille
                    esp_opt = 100000000000
                    n_plt = 100000 + 100000 * nump
                    tp_org = 00000
                    count_sh = esp_opt + n_plt + tp_org + nb_tri
                    nb_tri += 1

                    if (sc[x.id].name == 'Stp'):
                        liste = []
                        liste.append(nump)
                        liste.append(x.id)
                        liste.append(count_sh)

                        Dico_Stp[count_ID] = liste


                if (sc[x.id].name == 'S' or sc[x.id].name == 'Lf' or sc[x.id].name == 'In'):  # Tige
                    esp_opt = 100000000000
                    n_plt = 00000 + 100000 * nump
                    tp_org = 10000
                    count_sh = esp_opt + n_plt + tp_org + nb_tri
                    nb_tri += 1

                s_leg2.add(Shape(geometry=x.geometry, id=x.id, appearance=x.appearance))
                s_leg2[count_ID].setName(str(count_sh))

                Dico_conv[count_ID] = x.id

                count_ID += 1


    return Dico_Apex_ID, Dico_conv, Dico_In, Dico_Pet, Dico_Stp, s_leg2





def create_Sensors(runL):

    xyz = runL[1].dxyz
    pattern8 = runL[1].pattern8

    nb_vox = int(pattern8[1][0] / xyz[0])
    zmax_sc = int(max(runL[1].invar['Hplante'])) + 2

    s_capt=Scene()

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    points = [(0, 0, 0),(x,0,0), (x,y,0),(0,y,0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    carre = QuadSet(points, indices, normals, indices)

    ID_capt = 1
    dico_translat = {}
    dico_IDtoVox = {}
    dico_VoxtoID = {}

    for i in range(0,nb_vox):
        tx = i * x
        dico_VoxtoID[i] = {}

        for j in range(0, nb_vox):
            ty = j * y
            dico_VoxtoID[i][j] = {}

            for k in range(0, zmax_sc):
                tz = k * z
                dico_VoxtoID[i][j][k] = []
                dico_VoxtoID[i][j][k].append(ID_capt)

                Vox=Translated(geometry=carre, translation=(tx, ty, tz))
                s_capt.add(Shape(geometry=Vox, id=ID_capt))

                liste = []
                liste.append(tx)
                liste.append(ty)
                liste.append(tz)
                dico_translat[ID_capt] = liste

                liste = []
                liste.append(i)
                liste.append(j)
                liste.append(k)
                dico_IDtoVox[ID_capt] = liste

                ID_capt += 1

    #Viewer.display(s_capt)



    #Ajout des capteurs verticaux
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    points = [[(0, 0, 0), (x, 0, 0), (x, 0, z), (0, 0, z)], [(0, y, 0), (0, y, z), (x, y, z), (x, y, 0)], [(0, 0, 0), (0, y, 0), (0, y, z),(0, 0, z)], [(x, 0, 0),(x, y, 0), (x, y, z),(x, 0, z)]]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]
    for i in range(0,nb_vox):
        tx = i * x
        for j in range(0, nb_vox):
            ty = j * y
            for k in range(0, zmax_sc):
                tz = k * z

                for l in range(0, 4):
                    dico_VoxtoID[i][j][k].append(ID_capt)
                    carre = QuadSet(points[l], indices, normals, indices)

                    Vox = Translated(geometry=carre, translation=(tx, ty, tz))
                    s_capt.add(Shape(geometry=Vox, id=ID_capt))

                    liste=[]
                    liste.append(tx)
                    liste.append(ty)
                    liste.append(tz)
                    liste.append(l)
                    dico_translat[ID_capt] = liste

                    liste=[]
                    liste.append(i)
                    liste.append(j)
                    liste.append(k)
                    liste.append(l)

                    dico_IDtoVox[ID_capt] = liste

                    ID_capt += 1




    #Preparation capteurs
    pt_lst = s_capt[1].geometry.geometry.pointList
    idx_lst = s_capt[1].geometry.geometry.indexList

    Dico_Sensors = {}

    for x in s_capt:

        for i in range(0, len(idx_lst)):

            x11 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y11 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z11 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            x12 = pt_lst[idx_lst[i][1]][0] + dico_translat[x.id][0]
            y12 = pt_lst[idx_lst[i][1]][1] + dico_translat[x.id][1]
            z12 = pt_lst[idx_lst[i][1]][2] + dico_translat[x.id][2]

            x13 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y13 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z13 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            tple1 = []
            triangle = []
            tple1.append((x11, y11, z11))
            tple1.append((x12, y12, z12))
            tple1.append((x13, y13, z13))
            triangle.append(tple1)

            x21 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y21 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z21 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            x22 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y22 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z22 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            x23 = pt_lst[idx_lst[i][3]][0] + dico_translat[x.id][0]
            y23 = pt_lst[idx_lst[i][3]][1] + dico_translat[x.id][1]
            z23 = pt_lst[idx_lst[i][3]][2] + dico_translat[x.id][2]

            tple2 = []
            tple2.append((x21, y21, z21))
            tple2.append((x22, y22, z22))
            tple2.append((x23, y23, z23))
            triangle.append(tple2)

            Dico_Sensors[x.id]= triangle


    return Dico_Sensors, dico_VoxtoID, dico_IDtoVox, s_capt






def create_LstVox_V1(runL,Dico_val,nbplantes):

    xyz = runL[1].dxyz
    pattern8 = runL[1].pattern8
    origin_grid = runL[1].origin_grid
    na = runL[1].na

    Dico_ParTresh = runL[1].Dico_ParTresh

    lstParTresh = {}
    count = 0
    for x in range(0,nbplantes):

        for z in Dico_ParTresh[x]['A']:
            xvox = Dico_ParTresh[x]['A'][z][0][0]
            yvox = Dico_ParTresh[x]['A'][z][0][1]
            zvox = Dico_ParTresh[x]['A'][z][0][2]
            numVox = [xvox, yvox, zvox]
            lstParTresh[count] = numVox
            count += 1

        for z in Dico_ParTresh[x]['B']:
            z = z.split()
            xB = float(z[1])
            yB = float(z[2])
            zB = float(z[3])
            placement = array([xB, yB, zB])
            vox = WhichVoxel(placement, origin_grid, na, xyz)
            lstParTresh[count] = vox
            count += 1

        for z in Dico_ParTresh[x]['A2']:
            xvox = Dico_ParTresh[x]['A2'][z][0][0]
            yvox = Dico_ParTresh[x]['A2'][z][0][1]
            zvox = Dico_ParTresh[x]['A2'][z][0][2]



            numVox = [xvox,yvox,zvox]
            lstParTresh[count] = numVox
            count += 1


    Lst_Vox = {}
    for k in Dico_val:
        xa = Dico_val[k][2]
        ya = Dico_val[k][3]
        za = Dico_val[k][4]
        placement = array([xa, ya, za])
        vox = WhichVoxel(placement, origin_grid, na, xyz)
        Lst_Vox[k] = vox

    for count in lstParTresh:
        Lst_Vox[k+1+count] = lstParTresh[count]

    return  Lst_Vox






def create_Sensors_V2(runL,Dico_val,nbplantes):


    xyz = runL[1].dxyz
    pattern8 = runL[1].pattern8
    origin_grid = runL[1].origin_grid
    na = runL[1].na
    Lst_Vox = create_LstVox_V1(runL,Dico_val,nbplantes)


    s_capt = Scene()

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    points = [(0, 0, 0),(x,0,0), (x,y,0),(0,y,0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    carre = QuadSet(points, indices, normals, indices)

    ID_capt = 1
    dico_translat = {}
    dico_IDtoVox = {}
    dico_VoxtoID = {}

    for l in Lst_Vox:

        i = Lst_Vox[l][0]
        j = Lst_Vox[l][1]

        if (Lst_Vox[l][2] != -1):
            k = len(runL[1].res_rfr) -1 - Lst_Vox[l][2]
        else :
            k = k+1#k = Boite la plus haute de la plante ????


        tx = i * x
        try :
            dico_VoxtoID[i] == {}
        except:
            dico_VoxtoID[i] = {}

        ty = j * y
        try :
            dico_VoxtoID[i][j] == {}
        except:
            dico_VoxtoID[i][j] = {}

        tz = k * z
        try:
            dico_VoxtoID[i][j][k] == {}
        except:
            dico_VoxtoID[i][j][k] = []

        dico_VoxtoID[i][j][k].append(ID_capt)

        Vox = Translated(geometry=carre, translation=(tx, ty, tz))
        s_capt.add(Shape(geometry=Vox, id=ID_capt))

        liste = []
        liste.append(tx)
        liste.append(ty)
        liste.append(tz)
        dico_translat[ID_capt] = liste

        liste = []
        liste.append(i)
        liste.append(j)
        liste.append(k)
        dico_IDtoVox[ID_capt] = liste

        ID_capt += 1

        #ajout face sup
        tz = (k+1) * z
        try:
            dico_VoxtoID[i][j][k+1] == {}
        except:
            dico_VoxtoID[i][j][k+1] = []

        dico_VoxtoID[i][j][k+1].append(ID_capt)

        Vox = Translated(geometry=carre, translation=(tx, ty, tz))
        s_capt.add(Shape(geometry=Vox, id=ID_capt))

        liste = []
        liste.append(tx)
        liste.append(ty)
        liste.append(tz)
        dico_translat[ID_capt] = liste

        liste = []
        liste.append(i)
        liste.append(j)
        liste.append(k)
        dico_IDtoVox[ID_capt] = liste

        ID_capt += 1


    #Viewer.display(s_capt)



    #Ajout des capteurs verticaux
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    points = [[(0, 0, 0), (x, 0, 0), (x, 0, z), (0, 0, z)], [(0, y, 0), (0, y, z), (x, y, z), (x, y, 0)], [(0, 0, 0), (0, y, 0), (0, y, z),(0, 0, z)], [(x, 0, 0),(x, y, 0), (x, y, z),(x, 0, z)]]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    for m in Lst_Vox:

        i = Lst_Vox[m][0]
        j = Lst_Vox[m][1]
        k = len(runL[1].res_rfr) - Lst_Vox[m][2] -1

        tx = i * x
        ty = j * y
        tz = k * z
        for n in range(0, 4):
            dico_VoxtoID[i][j][k].append(ID_capt)
            carre = QuadSet(points[n], indices, normals, indices)

            Vox = Translated(geometry=carre, translation=(tx, ty, tz))
            s_capt.add(Shape(geometry=Vox, id=ID_capt))

            liste = []
            liste.append(tx)
            liste.append(ty)
            liste.append(tz)
            liste.append(n)
            dico_translat[ID_capt] = liste

            liste=[]
            liste.append(i)
            liste.append(j)
            liste.append(k)
            liste.append(n)

            dico_IDtoVox[ID_capt] = liste

            ID_capt += 1




    #Preparation capteurs
    pt_lst = s_capt[1].geometry.geometry.pointList
    idx_lst = s_capt[1].geometry.geometry.indexList

    Dico_Sensors = {}

    for x in s_capt:

        for i in range(0, len(idx_lst)):

            x11 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y11 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z11 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            x12 = pt_lst[idx_lst[i][1]][0] + dico_translat[x.id][0]
            y12 = pt_lst[idx_lst[i][1]][1] + dico_translat[x.id][1]
            z12 = pt_lst[idx_lst[i][1]][2] + dico_translat[x.id][2]

            x13 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y13 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z13 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            tple1 = []
            triangle = []
            tple1.append((x11, y11, z11))
            tple1.append((x12, y12, z12))
            tple1.append((x13, y13, z13))
            triangle.append(tple1)

            x21 = pt_lst[idx_lst[i][0]][0] + dico_translat[x.id][0]
            y21 = pt_lst[idx_lst[i][0]][1] + dico_translat[x.id][1]
            z21 = pt_lst[idx_lst[i][0]][2] + dico_translat[x.id][2]

            x22 = pt_lst[idx_lst[i][2]][0] + dico_translat[x.id][0]
            y22 = pt_lst[idx_lst[i][2]][1] + dico_translat[x.id][1]
            z22 = pt_lst[idx_lst[i][2]][2] + dico_translat[x.id][2]

            x23 = pt_lst[idx_lst[i][3]][0] + dico_translat[x.id][0]
            y23 = pt_lst[idx_lst[i][3]][1] + dico_translat[x.id][1]
            z23 = pt_lst[idx_lst[i][3]][2] + dico_translat[x.id][2]

            tple2 = []
            tple2.append((x21, y21, z21))
            tple2.append((x22, y22, z22))
            tple2.append((x23, y23, z23))
            triangle.append(tple2)

            Dico_Sensors[x.id]= triangle


    return Dico_Sensors, dico_VoxtoID, dico_IDtoVox, s_capt




















def Opt_And_Pattern(scene,xmax,ymax):

    # Definition des propietes optiques des feuilles et tiges
    opt = {'rs': {}, 'rc': {}}

    for obj in scene:
        if obj.name[0] != '0':
            if (obj.name[7] == '1'):
                opt['rc'][obj.id] = (0.10, 0.07, 0.10, 0.07)
                opt['rs'][obj.id] = (0.41, 0.43, 0.41, 0.43)

            if (obj.name[7] == '0'):
                opt['rc'][obj.id] = (0.10, 0.07)
                opt['rs'][obj.id] = (0.41, 0.43)

        # Propriete sol :
        if (obj.name[0] == '0'):
            opt['rc'][obj.id] = (0.1, )
            opt['rs'][obj.id] = (0.1, )

    # Definition du patern
    pat = (0, 0, xmax, ymax)


    return pat, opt



def Sky_turtle6():
    # Calcul poids ciel 5 direction
    sky_list = []

    ls_poids = [1. / 6.] + [5. / (6. * 4)] * 4
    alfa_turtle6 = 0.4637
    effet_sin = [sin(pi / 2)] + [sin(alfa_turtle6)] * 4
    x = array(ls_poids) * effet_sin
    ls_poids = x / sum(x)

    Azym = [0, 0, 90, 180, 270]

    for i in range(0, 5):
        azimuth = Azym[i]
        Energie = ls_poids[i]

        if i != 0:
            x_dir = cos(radians(azimuth)) * sin(alfa_turtle6)
            y_dir = sin(radians(azimuth)) * sin(alfa_turtle6)
            z_dir = cos(alfa_turtle6)
        else:
            x_dir = 0
            y_dir = 0
            z_dir = 1

        t_sky = tuple((float(Energie), tuple((float(x_dir), float(y_dir), float(-z_dir)))))

        sky_list.append(t_sky)

    return sky_list






def WhichVoxel(p, origin_grid, na, dxyz):
    """ en z, id=0 = haut du couvert """
    p1_rel = p - origin_grid
    vox = [int(p1_rel[0]//dxyz[0]), int(p1_rel[1]//dxyz[1]), int(-p1_rel[2]//dxyz[2])] #// division entiere

    #test si dans grille et sinon retouve indice correspondant
    test = [0 <= vox[0] <na[0], 0 <= vox[1] <na[1], 0 <= vox[2] <na[2]]
    if test != [True,True,True]:#si hors grille
        #recup les indices du fautif
        matches = [i for i in range(3) if test[i]==False]
        #faire le calcul du bon id equivalent
        for i in matches:
            if i!=2: #x,y -> couvert infini
                vox[i] = vox[i]%na[i] # marche pour > et <0 et na[i]=1
            else: #z ->
                if vox[i]<0:#au dessus: met dans la derniere strate
                    vox[i] = 0
                else: #en dessous: met dans la strate au dessus du sol
                    vox[i] = na[i]-1
    return vox



def schnute(x,a=3.09,b=1.59,c=0,d=1.12,x1=0.,x2=2.):
    """fonction schnute generale appliquee a x"""
    Y=pow((pow(c,b)+(pow(d,b)-pow(c,b)))*((1-exp(-a*(x-x1)))/(1-exp(-a*((x2-x1))))),(1/b))
    return Y










#FeedBack donnee caribu sur Lpy
def Dico_Reponse(Dico_val, runL, Dico_ValCapt_Direct_Rc, Dico_ValCapt_Direct_Rs, dico_VoxtoID, nbplantes):


    origin_grid = runL[1].origin_grid
    xyz = runL[1].dxyz
    na = runL[1].na

    Dico_rep = {}
    Dico_PARif = {}

    for nump in range(0,nbplantes+1):
        Dico_rep[nump] = {}
        Dico_PARif[nump] = {}


    lstId = []
    for k in Dico_val:

        xa = Dico_val[k][2]
        ya = Dico_val[k][3]
        za = Dico_val[k][4]
        placement = array([xa, ya, za])
        vox = WhichVoxel(placement, origin_grid, na, xyz)

        try:
            LstId_Capt = dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr) - vox[2]]
            LstId_Capt = unique(LstId_Capt)
        except:
            print('erreur avec le voxel {} plante {}').format(vox,Dico_val[k][0])

        # Ajout face du dessus
        try:
            LstId_Capt.append(dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr)-1  - vox[2]][0])
            LstId_Capt = unique(LstId_Capt)
        except:
            pass

        try:
            face = 0
            sommeZeta = 0
            for Id_Capt in LstId_Capt:
                Rc = Dico_ValCapt_Direct_Rc[Id_Capt]
                Rs = Dico_ValCapt_Direct_Rs[Id_Capt]
                Zeta = Rc / Rs
                sommeZeta += Zeta
                face += 1

            #Zeta = schnute(Rc)

            Zeta = sommeZeta / face
            # Ecriture des resultats
            nump = Dico_val[k][0]

            Dico_rep[nump][str(xa)] = {}
            Dico_rep[nump][str(xa)][str(ya)] = {}
            Dico_rep[nump][str(xa)][str(ya)][str(za)] = Zeta

            Dico_PARif[nump][str(xa)] = {}
            Dico_PARif[nump][str(xa)][str(ya)] = {}
            Dico_PARif[nump][str(xa)][str(ya)][str(za)] = Rc
        except:
            pass

    return Dico_rep, Dico_PARif











def PrePa_Reponse_PAR_Tresh(aggregated_direct, runL, dico_VoxtoID, Dico_Apex_ID):

    Dico_ParTresh = runL[1].Dico_ParTresh
    xyz = runL[1].dxyz
    surfaceVox = xyz[0] * xyz[1]
    Dico_ValCapt_Direct_Rc = aggregated_direct['rc']['sensors']['Ei']

    Dico_Direct_Rc = aggregated_direct['rc']['Ei_sup']

    Normalis = runL[1].meteo_j['I0']


    Dico_Rep = {}
    for nump in Dico_ParTresh :
        Dico_Rep[nump] = {}
        Dico_Rep[nump]['A'] = {}
        Dico_Rep[nump]['A2'] = {}
        Dico_Rep[nump]['B'] = {}

    for nump in Dico_ParTresh:
        try:

            LstID = Dico_ParTresh[nump]['A']
            for ID in LstID:
                try:
                    PAR = Dico_Direct_Rc[Dico_Apex_ID[ID]]
                    Dico_Rep[nump]['A'][ID] = PAR * Normalis * 2 # Probleme cari
                except:
                    vox = [Dico_ParTresh[nump]['A'][ID][0][0], Dico_ParTresh[nump]['A'][ID][0][1],Dico_ParTresh[nump]['A'][ID][0][2]]
                    LstId_Capt = dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr) - vox[2]]
                    # Ajout face du dessus
                    LstId_Capt.append(dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr)- vox[2]][0])
                    LstId_Capt = unique(LstId_Capt)
                    PAR = 0
                    for Id_Capt in LstId_Capt:
                        Rc = Dico_ValCapt_Direct_Rc[Id_Capt]
                        PAR += Rc

                    PAR = PAR /len(LstId_Capt)
                    Dico_Rep[nump]['A'][ID] = PAR * Normalis

        except:
           pass


        LstID = Dico_ParTresh[nump]['A2']
        for ID in LstID:

            try :

                vox = [Dico_ParTresh[nump]['A2'][ID][0][0], Dico_ParTresh[nump]['A2'][ID][0][1], Dico_ParTresh[nump]['A2'][ID][0][2]]
                LstId_Capt = dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr)-1 - vox[2]]
                # Ajout face du dessus
                LstId_Capt.append(dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr) - vox[2]][0])
                LstId_Capt = unique(LstId_Capt)
                PAR = 0
                for Id_Capt in LstId_Capt:
                    Rc = Dico_ValCapt_Direct_Rc[Id_Capt]
                    PAR += Rc

                PAR = PAR /len(LstId_Capt)

                Dico_Rep[nump]['A2'][ID] = PAR * Normalis

            except:
                pass

        LstID = Dico_ParTresh[nump]['B']
        for ID in LstID:
            try:
                vox = [Dico_ParTresh[nump]['B'][ID][0][0], Dico_ParTresh[nump]['B'][ID][0][1],0 ]
                LstId_Capt = dico_VoxtoID[vox[0]][vox[1]][0]
                # Ajout face du dessus
                LstId_Capt.append(dico_VoxtoID[vox[0]][vox[1]][0][0])
                LstId_Capt = unique(LstId_Capt)
                PAR = 0
                for Id_Capt in LstId_Capt:
                    Rc = Dico_ValCapt_Direct_Rc[Id_Capt]
                    PAR += Rc

                PAR = PAR /(len(LstId_Capt))

                Dico_Rep[nump]['B'][ID] = PAR * Normalis

            except:
                pass

    return Dico_Rep












def CompareVoxCaribu(Dico_val, runL, aggregated_direct, Dico_conv):

    origin_grid = runL[1].origin_grid
    xyz = runL[1].dxyz
    na = runL[1].na
    sc = runL[0]

    Dico_Rc = aggregated_direct['rc']['Ei_sup']
    Dico_Rs = aggregated_direct['rs']['Ei_sup']

    Dico_return = {}

    for k in Dico_val:

        xa = Dico_val[k][2]
        ya = Dico_val[k][3]
        za = Dico_val[k][4]
        placement = array([xa, ya, za])
        vox = WhichVoxel(placement, origin_grid, na, xyz)

        #Zeta Caribu Organe
        Id_Org = Dico_conv[k]
        Rc = Dico_Rc[Id_Org]
        Rs = Dico_Rs[Id_Org]

        Zeta = Rc / Rs

        # Zeta Lpy
        RFR = runL[1].res_rfr[vox[2]][vox[1]][vox[0]]
        ResTrans = runL[1].res_trans[vox[2]][vox[1]][vox[0]]


        liste = []
        liste.append(Id_Org)
        liste.append(sc[Id_Org].name)
        liste.append(RFR)
        liste.append(ResTrans)
        liste.append(Rc)
        liste.append(Zeta)

        Dico_return[k] = liste




    return Dico_return