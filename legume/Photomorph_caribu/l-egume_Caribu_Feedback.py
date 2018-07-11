import legume
from numpy import *
from generateScene import *
import Fonct_Comp as Fc
from openalea.lpy import *
from openalea.plantgl.all import *
#Imports Caribu
from params_Ciel import*
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun
from alinea.caribu.caribu_shell import *

path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
lpy_filename = path_+'\l-egume.lpy'

path_pois = path_+'\Photomorph_caribu\input\Position_Pois.csv'
fichier_pois = loadtxt(path_pois, skiprows=1, delimiter=';')

runL,lsys, axiom, nbplantes = Fc.Init_Lpy(lpy_filename=lpy_filename, fichier_pois=fichier_pois)


nb_iter = 180

for i in range(nb_iter):
    runL = run(runL[1], axiom=runL[0], nbstep=1)


    s_leg = runL[1].sceneInterpretation(runL[0]).deepcopy()

    Dico_Apex_ID, Dico_conv, Dico_In, Dico_Pet, Dico_Stp, s_leg2 = Fc.PrepareScene(scene=s_leg, runL=runL)


    Dico_val = {}

    for x in Dico_Pet:
        Dico_val[x] = Dico_Pet[x]
    for x in Dico_In:
        Dico_val[x] = Dico_In[x]

    if (Dico_val!= {} ):
        Dico_Sensors, dico_VoxtoID, dico_IDtoVox, s_capt = Fc.create_Sensors_V2(runL, Dico_val, nbplantes)

        #Definition du pattern et des proprietes optique
        s = s_leg2
        pat, opt = Fc.Opt_And_Pattern(scene=s, xmax=lsys.cote, ymax=lsys.cote)

        #Calculs Caribu Direct
        #Soleil vertical
        soleil = [(1.0, (0., -0.0, -1))]
        d_sphere = 0.05 # simu 1 = 0.1  simu 2 = 0.05m


        cc_scene = CaribuScene(scene = s, opt = opt, light = soleil, pattern = pat)
        _, aggregated_direct = cc_scene.run(direct = False, infinite = False, split_face = True, sensors= Dico_Sensors, d_sphere= d_sphere)


        #Resultats
        Dico_ValCapt_Direct_Rc = aggregated_direct['rc']['sensors']['Ei']
        Dico_ValCapt_Direct_Rs = aggregated_direct['rs']['sensors']['Ei']



        #Prepa Dico de reponse
        try :
            Dico_rep, Dico_PARif = Fc.Dico_Reponse(Dico_val, runL, Dico_ValCapt_Direct_Rc, Dico_ValCapt_Direct_Rs, dico_VoxtoID, nbplantes)
        except:
            pass
        try :
            Dico_rep_PAR = Fc.PrePa_Reponse_PAR_Tresh(aggregated_direct, runL, dico_VoxtoID, Dico_Apex_ID)
        except:
            pass

        Dico_Compare = Fc.CompareVoxCaribu(Dico_val, runL, aggregated_direct, Dico_conv)
        fout = open('Compare_Voxel_CaribuOrgane.dat', 'a')
        for k in Dico_Compare:
            fout.write("%s;%s;%s;%s;%s;%s;%s \n" % (k, Dico_Compare[k][0], Dico_Compare[k][1], Dico_Compare[k][2], Dico_Compare[k][3], Dico_Compare[k][4], Dico_Compare[k][5]  ))

        #injection
        runL[1].Dico_rep = Dico_rep
        runL[1].Dico_rep_PAR = Dico_rep_PAR
        runL[1].Dico_PARif = Dico_PARif




Viewer.display(s_leg)


print("Simulation finie")
#Recuperer les sorties !!

s_leg.save('Pois_Caribu_dsphere=0.05_DOY 185.bgeom')
Viewer.display(s_leg)