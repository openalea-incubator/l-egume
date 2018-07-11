import legume
from numpy import *
from generateScene import *
import Fonct_Comp as Fc
from openalea.lpy import *
from openalea.plantgl.all import *



path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
lpy_filename = path_+'\l-egume.lpy'

path_pois = path_+'\Photomorph_caribu\input\Position_Pois.csv'
fichier_pois = loadtxt(path_pois, skiprows=1, delimiter=';')

runL, lsys, axiom, nbplantes = Fc.Init_Lpy(lpy_filename = lpy_filename, fichier_pois = fichier_pois)



nb_iter = 180
for i in range(nb_iter):
    runL = run(runL[1], axiom=runL[0], nbstep=1)




s_leg = runL[1].sceneInterpretation(runL[0]).deepcopy()
#Viewer.display(s_leg)

Dico_Apex_ID, Dico_conv, Dico_In, Dico_Pet, Dico_Stp, s_leg2 = Fc.PrepareScene(scene=s_leg, runL=runL)
#Viewer.display(s_leg2)


Dico_val = {}

for x in Dico_Pet:
    Dico_val[x] = Dico_Pet[x]
for x in Dico_In:
    Dico_val[x] = Dico_In[x]

#Ecriture du fichier de capteur virtuels :
    #Definir les voxels comme definis dans model Grassland
Dico_Sensors, dico_VoxtoID, dico_IDtoVox, s_capt =Fc.create_Sensors_V2(runL,Dico_val,nbplantes)



#Imports Caribu

from params_Ciel import*
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun
from alinea.caribu.caribu_shell import *


#Definition du pattern et des proprietes optique
s = s_leg2
pat, opt = Fc.Opt_And_Pattern(scene=s, xmax=lsys.cote, ymax=lsys.cote)



#Caribu Direct

#heure = 12
#DOY = runL[1].DOY

#soleil = Gensun.Gensun()(1, DOY, heure, 43)
#soleil = GetLightsSun.GetLightsSun(soleil)
#soleil = soleil.split(' ')
#soleil = [tuple((float(soleil[0]), tuple((float(soleil[1]), float(soleil[2]), float(soleil[3])))))]

#Soleil vertical
soleil = [(1.0, (0., -0.0, -1))]

cc_scene = CaribuScene(scene = s, opt = opt, light = soleil, pattern = pat)
_, aggregated_direct = cc_scene.run(direct = False, infinite = False, split_face = True, sensors = Dico_Sensors)

#Resultats direct
Dico_ValCapt_Direct_Rc = aggregated_direct['rc']['sensors']['Ei']
Dico_ValCapt_Direct_Rs = aggregated_direct['rs']['sensors']['Ei']



#Tests surface de feuille
lst_surf = []
for x in Dico_Stp:
    lst_surf.append(aggregated_direct['rc']['area'][Dico_Stp[x][1]])

moy_surf = mean(lst_surf)
min_surf = min(lst_surf)
max_surf = max(lst_surf)
xyz = runL[1].dxyz
surfaceVox = xyz[0] *xyz[1]

PourcentMoy = 100 * 2 * moy_surf / surfaceVox
PourcentMin = 100 * 2 * min_surf / surfaceVox
PourcentMax = 100 * 2 * max_surf / surfaceVox



origin_grid = runL[1].origin_grid
xyz = runL[1].dxyz
na = runL[1].na

lst_ResTrans = []
lst_vox =[]
for k in Dico_Pet:

    xa = Dico_Pet[k][2]
    ya = Dico_Pet[k][3]
    za = Dico_Pet[k][4]
    placement = array([xa,ya,za])
    vox= Fc.WhichVoxel(placement,origin_grid,na,xyz)
    lst_vox.append(str(vox))

    ResTrans = runL[1].res_trans[vox[2]][vox[1]][vox[0]]

    lst_ResTrans.append(ResTrans)


moy_ResTrans = mean(lst_ResTrans)
min_ResTrans = min(lst_ResTrans)
max_ResTrans = max(lst_ResTrans)


NbVox = len(lst_vox)
NbVox_Unique = len(unique(lst_vox))


LAIS = runL[1].m_lais[0]
lstLAI = []

for z in range(0,len(LAIS)):

    for y in range(0,len(LAIS[z])):

        for x in range(0,len(LAIS[z][y])):
            lstLAI.append(sum(LAIS[z][y][x]))





#Caribu Diffus

#sky_list = determine_sky_Caribu()
sky_list = Fc.Sky_turtle6()

cc_scene = CaribuScene(scene = s, opt = opt, light = sky_list, pattern = pat)
_, aggregated_diffus = cc_scene.run(direct = False, infinite = False, split_face = True, sensors=Dico_Sensors)

#Resultats diffus

Dico_ValCapt_Diffus_Rc = aggregated_diffus['rc']['sensors']['Ei']
Dico_ValCapt_Diffus_Rs = aggregated_diffus['rs']['sensors']['Ei']







#Calculs et sorties

PourcDifRs=0.225
PourcDifRc= 0.27

lstId = []

origin_grid = runL[1].origin_grid
xyz = runL[1].dxyz
na = runL[1].na


NormalResTrans = runL[1].meteo_j['I0']* runL[1].surf_refVOX

for k in Dico_Pet:

    xa = Dico_Pet[k][2]
    ya = Dico_Pet[k][3]
    za = Dico_Pet[k][4]
    placement = array([xa,ya,za])
    vox= Fc.WhichVoxel(placement,origin_grid,na,xyz)

    LstId_Capt= dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr)-1-vox[2]]
    #Ajout face du dessus
    LstId_Capt.append(dico_VoxtoID[vox[0]][vox[1]][len(runL[1].res_rfr)-vox[2]][0])
    LstId_Capt = unique(LstId_Capt)

    face = 0
    sommeZeta = 0
    sommeZeta2 = 0
    sommeZetaDirectDiffus = 0
    sommeRc = 0
    for Id_Capt in LstId_Capt:

        Rc = Dico_ValCapt_Direct_Rc[Id_Capt]
        Rs = Dico_ValCapt_Direct_Rs[Id_Capt]
        Zeta = Rc / Rs

        Zeta2 = Fc.schnute(Dico_ValCapt_Direct_Rc[Id_Capt])

        Rc_dir = (1 - PourcDifRc) * Dico_ValCapt_Direct_Rc[Id_Capt]
        Rs_dir = (1 - PourcDifRs) * Dico_ValCapt_Direct_Rs[Id_Capt]
        Rc_dif = PourcDifRc * Dico_ValCapt_Diffus_Rc[Id_Capt]
        Rs_dif = PourcDifRs * Dico_ValCapt_Diffus_Rs[Id_Capt]

        RcDirectDiffus  = Rc_dir + Rc_dif
        RsDirectDiffus  = Rs_dir + Rs_dif
        ZetaDirectDiffus = RcDirectDiffus / RsDirectDiffus

        #Zeta Lpy
        RFR = runL[1].res_rfr[vox[2]][vox[1]][vox[0]]
        ResTrans = runL[1].res_trans[vox[2]][vox[1]][vox[0]]
        surface = LAIS[vox[2]][vox[1]][vox[0]]
        lstId.append(Id_Capt)
        #print Id_Capt, '   ', Zeta, '  ', RFR, '  ', Zeta2,'   ',face
        face += 1
        sommeRc += Rc
        sommeZeta += Zeta
        sommeZeta2 += Zeta2
        sommeZetaDirectDiffus += ZetaDirectDiffus
    #print Id_Capt, '   ', sommeZeta/face, '   ', sommeZeta/face,'   ',sommeZeta2/face, '   ', RFR, '  ',RES_TRANS
    #print sommeRc/face, '   ',sommeZeta2/face, '   ', sommeZeta/face,'   ',sommeZetaDirectDiffus/face, '   ', RFR, '  ',ResTrans/NormalResTrans, surface










#Dico de feedback
Dico_rep, Dico_PARif = Fc.Dico_Reponse(Dico_val, runL, Dico_ValCapt_Direct_Rc, Dico_ValCapt_Direct_Rs, dico_VoxtoID, nbplantes)
Dico_rep_PAR = Fc.PrePa_Reponse_PAR_Tresh(aggregated_direct, runL, dico_VoxtoID, Dico_Apex_ID)

#injection
runL[1].Dico_rep = Dico_rep
runL[1].Dico_rep_PAR = Dico_rep_PAR



#Creation dico comparaison et ecriture fichier
Dico_Compare = Fc.CompareVoxCaribu(Dico_val, runL, aggregated_direct, Dico_conv)


filename = path_+'\Photomorph_caribu\output\Compare_Voxel_CaribuOrgane.dat'
fout = open(filename, 'w')
for k in Dico_Compare:
    fout.write("%s;%s;%s;%s;%s;%s;%s \n" % (
    k, Dico_Compare[k][0], Dico_Compare[k][1], Dico_Compare[k][2], Dico_Compare[k][3], Dico_Compare[k][4],
    Dico_Compare[k][5]))




#Remarques generale sur le code
# Changement dans LPY : Precede de #Ajout Quentin
# Ligne 1416 de Lpy la surface de cotyledon a ete commentee pour avoir le meme calcul de surface sur les deux test caribu et voxel

