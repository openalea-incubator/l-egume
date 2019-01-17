#batch pour L-egume 23/12/17
#version GL utilisee pour melanges binaires (8*8)

#import des fonctions de Quentin
from generateScene import run
import Fonct_Comp as Fc

#Imports Caribu comme Quentin
from params_Ciel import*
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun
from alinea.caribu.caribu_shell import *
from openalea.plantgl.all import *

#import the modules necessary to initiate the L-systems
from openalea.lpy import *
import multiprocessing

import os
import sys

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'

print ('path', path_)

sys.path.insert(0, path_)
import IOxls
import IOtable


#lecture de la liste des usm
#path_ = r'H:\devel\grassland\grassland\L-gume'
mn_path = os.path.join(path_,'multisim','liste_usms_mix.xls')#(path_,'liste_usms_mix_these lucas.xls')#
ongletBatch = 'SimTest'#'complement'#'Feuil1'#'Sensi'#
usms = IOxls.xlrd.open_workbook(mn_path)
ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))


#lecture scenario de changement des parametres (passer dans le script -> pourrait ne le faire que dans le batch)
mn_sc = os.path.join(path_,'liste_scenarios.xls')#(path_,'liste_scenarios_these_lucas.xls')#
#ongletBatch = 'Feuil1'
#usc = IOxls.xlrd.open_workbook(mn_sc)
#ls_sc = IOtable.conv_dataframe(IOxls.get_xls_col(usc.sheet_by_name(ongletBatch)))





#cree la liste de L-systems et liste des noms
testsim={}
names = []
for i in range(1):#len(ls_usms['ID_usm'])):
    if int(ls_usms['torun'][i]) == 1:#si 1 dans la colonne 'torun' l'ajoute a la liste
        name = str(int(ls_usms['ID_usm'][i]))+'_'+str(ls_usms['l_system'][i])[0:-4]
        seednb = int(ls_usms['seed'][i])
        names.append(name)
        path_lsys = os.path.join(path_, str(ls_usms['l_system'][i]))
        testsim[name]=Lsystem(path_lsys)  #objet l-system
        
        
        #testsim[name].ongletM = str(ls_usms['ongletM'][i])
        meteo_path_ =  os.path.join(path_, 'input',str(ls_usms['meteo'][i]))
        ongletM_ = str(ls_usms['ongletM'][i])
        testsim[name].meteo = IOxls.read_met_file(meteo_path_, ongletM_)

        #testsim[name].ongletMn = str(ls_usms['ongletMn'][i])
        mn_path_ = os.path.join(path_, 'input',str(ls_usms['mng'][i]))
        ongletMn_ = str(ls_usms['ongletMn'][i])
        testsim[name].mng = IOxls.read_met_file(mn_path_, ongletMn_)

        ini_path_ = os.path.join(path_, 'input',str(ls_usms['inis'][i]))
        ongletIni_ = str(ls_usms['ongletIn'][i])
        testsim[name].inis = IOxls.read_plant_param(ini_path_, ongletIni_)

        #testsim[name].ongletP = str(ls_usms['ongletP'][i])
        path_plante = os.path.join(path_, 'input',str(ls_usms['plante'][i]))
        ongletP = str(ls_usms['ongletP'][i])
        ongletPvois = str(ls_usms['ongletVoisin'][i])
        ### g4 = IOxls.read_plant_param(path_plante, ongletP)
        ### g5 = IOxls.read_plant_param(path_plante, ongletvois)
        
        #la, lire scenario et changer parametres
        idscenar1 = int(ls_usms['scenario1'][i])
        idscenar2 = int(ls_usms['scenario2'][i])
        ongletScenar2 = ongletPvois #fait porter les changements sur fichier parametre voisin
        ongletScenar1 = ongletP
        
        #nbcote=7 # a passer ext
        #testsim[name].ParamP = [g4]*int(ls_usms['nbplt'][i])#nbcote*nbcote
        optdamier = int(ls_usms['damier'][i])
        nbcote = int(ls_usms['nbcote'][i])
        ### testsim[name].ParamP = damier8(g4,g5,opt=optdamier)
        nommix = '_'+ongletP+'-'+ongletPvois+'_'+'damier'+str(optdamier)+'_scenario'+str(idscenar2)+'-'+str(idscenar1)
        
        testsim[name].ongletS = str(ls_usms['ongletS'][i])
        testsim[name].ongletP = ongletP
        testsim[name].ongletPvois = ongletPvois
        testsim[name].nbcote = nbcote
        testsim[name].cote = float(ls_usms['cote'][i])
        testsim[name].deltalevmoy = float(ls_usms['retard'][i])
        testsim[name].deltalevsd = float(ls_usms['sd_retard'][i])
        testsim[name].typearrangement = str(ls_usms['arrangement'][i])
        testsim[name].optdamier = optdamier
        testsim[name].idscenar1 = idscenar1
        testsim[name].idscenar2 = idscenar2
        testsim[name].ongletScenar2 = ongletScenar2
        testsim[name].ongletScenar1 = ongletScenar1
        testsim[name].Rseed = seednb
        testsim[name].DOYdeb = int(ls_usms['DOYdeb'][i])
        testsim[name].DOYend = int(ls_usms['DOYend'][i])
        
        #mise a jour derivartionLength & axiom
        testsim[name].derivationLength = 1#int(ls_usms['DOYend'][i]) - int(ls_usms['DOYdeb'][i])#derivationLength variable predefinie dans L-py
        nbplantes = nbcote*nbcote
        a=AxialTree()
        a.append(testsim[name].attente(1))
        for j in range(0,nbplantes):
           a.append(testsim[name].Sd(j))

        testsim[name].axiom = a#passe un axial tree, pas de chaine de caractere


        #path fichiers de sortie
        testsim[name].path_out = os.path.join(path_, str(ls_usms['folder_out'][i]))
        testsim[name].outvarfile = 'toto_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.csv'
        testsim[name].lsorgfile = 'lsAxes_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.csv'
        testsim[name].outHRfile = 'outHR_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.csv'
        testsim[name].resrootfile = 'resroot_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.csv'
        testsim[name].outBilanNfile = 'BilanN_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.csv'
        testsim[name].outimagefile = 'scene_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.bmp'#'scene.bmp'

        #plante si dossier out pas cree
        #pourrait faire la lecture les ls_usm directement dans le l-system pour faciliter...+

        print 'Fin initialisation '

nb_usms =len(names)#len(ls_usms['ID_usm'])#len(names)#




#function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    lsys =  testsim[names[n]]
    lstring = lsys.axiom
    nb_iter = 45

    for i in range(nb_iter):
        print 'iter ',i,n
        lstring = lsys.derive(lstring, i, 1)
        #runL = run(testsim[names[n]], axiom=axiom, nbstep=1)
        print 'ici',lsys.cote
        #axiom = runL[0]

        s_leg = lsys.sceneInterpretation(lstring)

        s_leg.save("s_leg.bgeom")

        Dico_Apex_ID, Dico_conv, Dico_In, Dico_Pet, Dico_Stp, s_leg2 = Fc.PrepareScene(scene=s_leg, runL=lsys)


        Dico_val = {}

        for x in Dico_Pet:
            Dico_val[x] = Dico_Pet[x]
        for x in Dico_In:
            Dico_val[x] = Dico_In[x]

        if (Dico_val!= {} ):
            Dico_Sensors, dico_VoxtoID, dico_IDtoVox, s_capt = Fc.create_Sensors_V2(lsys, Dico_val, nbplantes)

            #Definition du pattern et des proprietes optique
            s = s_leg2
            pat, opt = Fc.Opt_And_Pattern(scene=s, xmax=lsys.cote, ymax=lsys.cote)

            #Calculs Caribu Direct
            #Soleil vertical
            soleil = [(1.0, (0., -0.0, -1))]#Fc.Sky_turtle6()
            d_sphere = 0.05 # simu 1 = 0.1  simu 2 = 0.05m

            print 'Calcul Caribu'
            cc_scene = CaribuScene(scene = s, opt = opt, light = soleil, pattern = pat)
            _, aggregated_direct = cc_scene.run(direct = False, infinite = False, split_face = True, sensors= Dico_Sensors, d_sphere= d_sphere)


            #Resultats
            Dico_ValCapt_Direct_Rc = aggregated_direct['rc']['sensors']['Ei']
            Dico_ValCapt_Direct_Rs = aggregated_direct['rs']['sensors']['Ei']



            #Prepa Dico de reponse
            try :
                Dico_rep, Dico_PARif = Fc.Dico_Reponse(Dico_val, lsys, Dico_ValCapt_Direct_Rc, Dico_ValCapt_Direct_Rs, dico_VoxtoID, nbplantes)
            except:
                pass
            try :
                Dico_rep_PAR = Fc.PrePa_Reponse_PAR_Tresh(aggregated_direct, lsys, dico_VoxtoID, Dico_Apex_ID)
            except:
                pass

            Dico_Compare = Fc.CompareVoxCaribu(Dico_val, lsys, aggregated_direct, Dico_conv)
            fout = open('Compare_Voxel_CaribuOrgane.dat', 'a')
            for k in Dico_Compare:
                fout.write("%s;%s;%s;%s;%s;%s;%s \n" % (k, Dico_Compare[k][0], Dico_Compare[k][1], Dico_Compare[k][2], Dico_Compare[k][3], Dico_Compare[k][4], Dico_Compare[k][5]  ))

            #injection
            lsys.Dico_rep = Dico_rep
            lsys.Dico_rep_PAR = Dico_rep_PAR
            lsys.Dico_PARif = Dico_PARif

            print 'Pas de temps ',lsys.DOY
    #testsim[names[n]].derive()
    testsim[names[n]].clear()
    print(''.join((names[n]," - done")))

#run the L-systems

### Ancienne methode ###
#if __name__ == '__main__':
#    multiprocessing.freeze_support()
#    CPUnb=multiprocessing.cpu_count()-1###nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
#    print 'nb CPU: '+str(CPUnb)
#    #for i in range(nb_usms):#sans parallellisme
#    #    runlsystem(i)
#    for i in range(int(((nb_usms-1)/CPUnb)+1)):#avec parrallellisme
#        pool = multiprocessing.Pool(processes=CPUnb)#doit etre superieur ou egal au range... PC lucas = 8 cpus.
#        if i<int(((nb_usms-1)/CPUnb)+1)-1:
#            pool.map(runlsystem, range(CPUnb*(i),CPUnb*(i+1)))
#            pool.close()
#            pool.join()
#        else:
#            pool.map(runlsystem, range(CPUnb*(i),nb_usms))
#            pool.close()
#            pool.join()
#    #Parallel(n_jobs=multiprocessing.cpu_count())(delayed(montrucaparalleliser)(sp) for sp in range(6))
#    #pool=multiprocessing.Pool(6)
#    #pool.map(parallel, range(6))


### Nouvelle methode ###
if __name__ == '__main__':
    #multiprocessing.freeze_support()
    #CPUnb=multiprocessing.cpu_count()-1 #nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
    print 'debiut'#nb CPU: '+str(CPUnb)
    #pool = multiprocessing.Pool(processes=CPUnb)
    for i in range(1):#int(nb_usms)):
        runlsystem(i)
        #pool.apply_async(runlsystem, args=(i,)) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
    pool.close()
    pool.join()