#### initialisation functions for l-egume

import os
import sys

try:
    import legume

    path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'  # r'C:\devel\grassland'


sys.path.insert(0, path_)

import numpy as np
from scipy import *
from copy import deepcopy
import string
import time
import pandas as pd

try:
    from .soil3ds import soil_moduleN as solN #import de la version develop si module soil3ds est installe
    #import soil_moduleN3 as solN
except:
    import soil_moduleN3 as solN #soil_moduleN2_bis as solN #! renommer car dans nouvelle version Lpy, mot module est reserve et fait planter!

try:
    from .riri5 import RIRI5 as riri #import de la version develop si module soil3ds est installe
except:
    import RIRI5 as riri

import RootDistrib as rtd
import RootMorpho2 as rt
import ShootMorpho as sh
import IOtable
import IOxls



def init_glob_variables_simVGL(meteo, mng, DOYdeb, path_station, ongletSta):
    """ Initialise global variables used within the L-egume L-system """

    ## station
    station = IOxls.read_plant_param(path_station, ongletSta)  # ou lire ds fichier inis

    ## meteo du jour
    DOY = DOYdeb
    meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay', 'RG', 'Et0', 'Precip', 'Tmin', 'Tmax', 'Tsol'], 'DOY', val=DOY)
    meteo_j['I0'] = [0.48 * meteo_j['RG'][0] * 10000 / (3600 * 24)]  # flux PAR journalier moyen en W.m-2 / RG en j.cm-2
    mng_j = IOxls.extract_dataframe(mng, ['Coupe', 'Irrig', 'FertNO3', 'FertNH4', 'Hcut'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k] = meteo_j[k][0]
    for k in list(mng_j.keys()): mng_j[k] = mng_j[k][0]
    meteo_j['durjour'] = sh.DayLength(station['latitude'], sh.DecliSun(DOY % 365))

    TT = 0
    TTsol = 0
    STEPS_ = meteo_j['TmoyDay'] - 5.  # dTT(meteo_j['TmoyDay'], [ParamP[0]['Tdev']])# #variable remise  ajour chaque jour
    STEPSsol_ = meteo_j['Tsol'] - 5.
    ls_epsi = [0.]  # [0.4, 0.4]##!! correspondance avec les nb de root systems!

    ##coupe
    TT_repousse = 0  # TT de la derniere coupev #resoud pb: TT utilise pour LAI pour NI revenait jamais a zero
    isTTcut = False
    wasTTcut = False
    # TTcutFreq = 18.*32#15.*32 #phyllochron
    isRegrowth = False  # indicateur: est on a la pousse initiale ou bien plus tard?
    Hcut = 1.  # 3.#simple initialisation : est passe en lecture fichier management
    cutNB = 0

    ## divers
    start_time, past_time = time.time(), 0.  # pour recuperer temps de calcul

    # distribution des retard a levee
    # deltalevmoy = 30 #degre.jours
    # deltalevsd = 15
    # pour gerer un decalage a la levee/reprise
    # tir1, tir2 = [], []
    # for i in range(1000):
    #  tir1.append(max(0,random.gauss(deltalevmoy,deltalevsd)))#test_retard.append(random.uniform(0,60))
    #  tir2.append(max(0,random.gauss(deltalevmoy,deltalevsd)))#test_retard.append(random.uniform(0,60))

    # test_retard = [tir1, tir2] #pour gerer deux especes

    return DOY, TT, TTsol, meteo_j, mng_j, STEPS_, STEPSsol_, ls_epsi, TT_repousse, isTTcut, wasTTcut, isRegrowth, cutNB, Hcut, start_time, past_time, station




def init_sol_fromLpy(inis, meteo_j, par_sol, par_SN, Lsol, largsol, discret_solXY, dz_sol, pattern8, opt_residu, obstarac=None):
    """ 3DS soil initialisation from L-py - 3 couches de sol avec propirete differentes maxi"""

    #Tsol lu dans meteo
    Tsol = meteo_j['Tsol']  # 15. #degresC

    # vecteurs d'initialisation du sol (pour 3 couches maxi)
    num_nb = list(map(int, inis['num_nb']))  # [6,6,18] #nbr de couche de chaque num de sol
    ncouches_sol = int(inis['ncouches_sol'])#num_nb[0] + num_nb[1] + num_nb[2]
    vsoilnumbers = [1] * num_nb[0] + [2] * num_nb[1] + [3] * num_nb[2]  # convention autorise 3 types d'horizon max
    # vDA = [par_SN['DA'][0]]*num_nb[0] + [par_SN['DA'][1]]*num_nb[1] + [par_SN['DA'][2]]*num_nb[2] #densite apparente de sol
    vCN = [par_SN['CN0_30']] * num_nb[0] + [par_SN['CN30_60']] * num_nb[1] + [par_SN['CN60_90']] * num_nb[2]  # maxi 3 horizons
    vMO = [par_SN['MO0_30']] * num_nb[0] + [par_SN['MO30_60']] * num_nb[1] + [par_SN['MO60_90']] * num_nb[2]  # maxi 3 horizons
    vARGIs = [par_SN['ARGIs0_30']] * num_nb[0] + [par_SN['ARGIs30_60']] * num_nb[1] + [par_SN['ARGIs60_90']] * num_nb[2]
    vCALCs = [par_SN['CALCs']] * ncouches_sol
    vNH4 = inis['NH4']  # [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    vNO3 = inis['NO3']  # [0.]*ncouches_sol
    HRpinit = inis['HRp']  # []
    if min(HRpinit) < 0:  # code -1 pour pas d'initialisation
        HRpinit = []

    vDA = []
    for i in vsoilnumbers:
        vDA.append(par_sol[str(i)]['DA'])

    # vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11
    # vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
    # vCN = [par_SN['CN0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vMO = [par_SN['MO0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vARGIs = [par_SN['ARGIs']]*ncouches_sol #maxi 90cm
    # vCALCs = [par_SN['CALCs']]*ncouches_sol
    # vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    # coeff = 0.#0.09#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
    # vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
    # vNO3 = array([16.96, 16.07, 15.17, 33.92, 33.92, 33.92, 33.92, 62.49, 82.13, 89.27, 76.77, 107.13, 124.98, 142.84, 124.98, 142.84, 160.69, 151.76, 151.76, 142.84, 178.55, 133.91, 98.20, 89.27, 83.92, 89.27, 73.20, 89.27, 87.45, 62.49])*coeff #issu du profil en sol nu
    # HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

    ## soil initialisation
    S = solN.SoilN(par_sol, par_SN, soil_number=vsoilnumbers,
                   dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1],
                         [dz_sol / 100.] * ncouches_sol], vDA=vDA, vCN=vCN, vMO=vMO, vARGIs=vARGIs, vNO3=vNO3,
                   vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol, pH=par_SN['pH'], ZESX=par_SN['ZESX'], CFES=par_SN['CFES'],
                   obstarac=obstarac, pattern8=pattern8)

    ## initialise humidite si fournit
    if HRpinit != []:  # initialise humidite si un vecteur est fourni
        S.init_asw(HRp_init=HRpinit)

    # Uval = 0.9*2.61#(epaisseur de sol* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
    # Uval = par_SN['q0']*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
    Uval = par_SN['q0']
    stateEV = [0., 0., 0.]  # pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
    HXs = par_sol[str(vsoilnumbers[0])]['teta_fc']  # humidite a la capacite au champ de l'horizon de surface
    b_ = solN.bEV(par_SN['ACLIMc'], par_SN['ARGIs'], HXs)  # HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63
    # print (par_SN['q0'],'Uval', ' b ', Uval, b_, sum(S.m_QH20fc[0]), surfsolref , (S.dxyz[2][0]*100.))

    return S, Tsol, Uval, stateEV, b_




def init_scene_fromLpy(ParamP, inis, cote, nbcote, station, lsidP, type='damier8'):
    # initialoise la scene L-egume: arrangement des plantes (carto), discretisation souterraine, discretisation aerienne
    # 1) CARTO
    distplantes = cote / nbcote  # 1. #cm

    ## pour ilot
    # carto = [array([0.,0.,0.]), array([distplantes,0.,0.]),array([-distplantes,0.,0.]),array([0.5*distplantes,0.866*distplantes,0.]),array([-0.5*distplantes,0.866*distplantes,0.]), array([0.5*distplantes,-0.866*distplantes,0.]), array([-0.5*distplantes,-0.866*distplantes,0.])]#, array([-10.,0.,0.]), array([0.,7.,0.])] #liste des localisations (1pt par plante) -> a lire en fichier #LF - cos (pi/3) = 0.5   sin (pi/3) = 0.866

    # carto = [array([0.,0.,0.]), array([distplantes,0.,0.]),array([-distplantes,0.,0.]),array([0.5*distplantes,0.866*distplantes,0.]),array([-0.5*distplantes,0.866*distplantes,0.]), array([0.5*distplantes,-0.866*distplantes,0.]), array([-0.5*distplantes,-0.866*distplantes,0.]),array([2*distplantes,0.,0.]),array([-2*distplantes,0.,0.]),array([2*0.5*distplantes,2*0.866*distplantes,0.]),array([-2*0.5*distplantes,2*0.866*distplantes,0.]), array([2*0.5*distplantes,-2*0.866*distplantes,0.]), array([-2*0.5*distplantes,-2*0.866*distplantes,0.])] #carto13

    ## pour champ en rangs
    # yyy = [-15]*9+[0]*9+[15]*9
    # xxx = range(-20,25,5)*3
    # yyy = [-15]*9+[-7.5]*9+[0]*9+[7.5]*9+[15]*9+[22.49]*9
    # xxx = range(-20,25,5)*6
    # yyy = [-15]*5+[0]*5+[15]*5
    # xxx = range(-20,25,10)*3
    # yyy = [-15]*2+[0]*2+[15]*2
    # xxx = range(-20,25,25)*3
    # yyy = [-15]*1+[0]*1+[15]*1
    # xxx = range(-20,25,60)*3
    # yyy = [-15]*23+[-7.5]*23+[0]*23+[7.5]*23+[15]*23+[22.49]*23
    # xxx = range(-20,25,2)*6

    # pour grand rhizotron
    # yyy = [-12.25, 4.4, 21.05]
    # xxx = [-10.75, -8.35, -5.95, -3.55, -1.15, 1.25, 3.65, 6.05, 8.45, 10.85]

    if type == 'row4':  # pour carre 4 rangs heterogenes
        Param_, carto = sh.row4(1, 2, Lrow=cote, nbprow=nbcote)
    elif type == 'random8' or type == 'random8':
        # pour tirage random
        carto = sh.random_planter(nbcote * nbcote, cote, cote)
    else:
        # pour carre distance homogene
        carto = sh.regular_square_planter(nbcote, distplantes)

    # reduit a une espece si veut simul separee
    if type == 'damier8_sp1' or type == 'damier8_sp2' or type == 'damier16_sp1' or type == 'damier16_sp2' or 'row4_sp1' or 'row4_sp2':
        carto = sh.reduce_carto(carto, lsidP)

    # 2) definition du pattern et discretisation sol
    pattern8 = [[0, 0], [cote, cote]]
    # pattern8 =[[min(xxx)-dp,min(yyy)-dp], [max(xxx)+dp,max(yyy)+dp]]
    # [[-2.5,-2.5], [2.5,2.5]]#[[-5.,-5.],[5.,5.]]#[[-2.5,-2.5], [5,5]]#[[-12.5,-12.5],[12.5,12.5]]
    # pattern8 = [[-22.3/2.,-49.5/2.], [22.3/2.,49.5/2.]]#pattern.8 rhizotron equivalent (cm)

    Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
    largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
    surfsolref = Lsol * largsol  # m2
    dz_sol = inis['dz_sol']  # 4.#5. #cm
    ncouches_sol = int(inis['ncouches_sol'])  # 4#10#30
    prof_sol_max = ncouches_sol * dz_sol  # 80.

    discret_solXY = list(map(int, inis['discret_solXY']))  # [10,10]# nb de discretisation du sol en X et en Y
    lims_sol = rtd.lims_soil(pattern8, dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1], [dz_sol / 100.] * ncouches_sol])

    # 3) discretisation au niveau aerien
    # creation des grid3D pour calcul de rayonnement
    ls_gammagroup = list(map(int, riri.get_lsparami(ParamP, 'gammagroup')))
    setp = list(set(ls_gammagroup))  # set equivalent fonction r.unique!
    n_gamagroup = len(setp)

    for nump in range(len(ParamP)):  # ajout a chaque plante de son id grille dans ParamP
        ParamP[nump]['id_grid'] = setp.index(int(ParamP[nump]['gammagroup']))  # pour savoir id de grille de la plante

    na, dxyz, lims_aer, origin_grid, surf_refVOX = riri.def_na_lims(pattern8, station['dz_aerien'], station['Hmaxcouv'], opt=station['opt1D3D'])  # '1D'
    # na, dxyz, lims_aer, origin_grid, surf_refVOX = riri.def_na_lims(pattern8, dz_aerien, Hmaxcouv,opt='1D') #version 1D, comme avant

    m_lais = np.zeros([n_gamagroup, na[2], na[1], na[0]])  # ngamma, Z,Y,,X
    m_lais_construct = deepcopy(m_lais)  # pour contsruction du m_lais a t+1
    m_laiPlt = np.zeros([len(ParamP), na[2], na[1], na[0]])  # Distrib lai des plantes indivs: nump, Z,Y,,X
    triplets = riri.get_ls_triplets(m_lais[0], opt=station['sky'])  # 'VXpXmYpYm')#opt='V')#

    # liste des k_teta (coeff extinction directionnels) par entite
    ls_dif = []
    # prepa des k_teta par entite
    for i in setp:
        nump = ls_gammagroup.index(i)  # retrouve numero de plante de premiere occurence de gammagroup
        ls_dif.append(ParamP[nump]['k_teta_distf'])

    # print('triplets', len(triplets))
    # print('ls_dif', len(ls_dif), ls_dif)

    res_trans, res_abs_i = [], []
    res_rfr = []

    return carto, distplantes, pattern8, Lsol, largsol, surfsolref, dz_sol, ncouches_sol, prof_sol_max, discret_solXY, lims_sol, ls_gammagroup, setp, n_gamagroup, na, dxyz, lims_aer, origin_grid, surf_refVOX, m_lais, m_lais_construct, m_laiPlt, triplets, ls_dif, res_trans, res_abs_i, res_rfr
    # par logique, passer cote et nbcote dans ini?



def init_plant_residues_fromParamP(S, opt_residu, ParamP):
    """ separate initialisation of plant residue from plant files / lystem"""
    # initialise soil residues from a ParamP of plants + update ParamP[nump]['CNRES']

    if opt_residu == 1:  # initialisatio de residus

        # number of groupes?
        ls_groupres = list(map(int, riri.get_lsparami(ParamP, 'groupe_resid')))
        setg = list(set(ls_groupres))  # set equivalent fonction r.unique!
        n_groupres = len(setg)
        nbplantes = len(ls_groupres)
        # recup 1 jeu de param pour chaque groupe
        CNRES, CC, WC, Nmires = [], [], [], []
        for i in range(len(setg)):
            for nump in range(nbplantes):
                if ParamP[nump]['groupe_resid'] == setg[i]:
                    ParamP[nump]['CNRES'] = [ParamP[nump]['CNRESlf'], ParamP[nump]['CNRESst'], ParamP[nump]['CNRESr'], ParamP[nump]['CNRESpiv']]
                    CNRES = CNRES + ParamP[nump]['CNRES']  # ajout des 4 classes pardefaut
                    CC = CC + ParamP[nump]['CC']
                    WC = WC + ParamP[nump]['WC']
                    Nmires = Nmires + ParamP[nump]['Nmires']
                    # print ('par',CNRES, CC, WC, Nmires)
                    break  # s'arrete a premiere plante de ce groupe

        if len(setg) == 1:  # si 1 seul grope, met qd meme un deuxieme residu de meme type pour pas planter
            for nump in range(nbplantes):
                if ParamP[nump]['groupe_resid'] == setg[0]:
                    ParamP[nump]['CNRES'] = [ParamP[nump]['CNRESlf'], ParamP[nump]['CNRESst'], ParamP[nump]['CNRESr'],
                                             ParamP[nump]['CNRESpiv']]
                    CNRES = CNRES + ParamP[nump]['CNRES']
                    CC = CC + ParamP[nump]['CC']
                    WC = WC + ParamP[nump]['WC']
                    Nmires = Nmires + ParamP[nump]['Nmires']
                    # print ('par',CNRES, CC, WC, Nmires)
                    break
        # Nmires not used???

        # distrib dans le sol en dur!
        nb_res = len(CNRES)  # 4 types de residus par espece (4 compatiment du papier) * 2 especes #pas utilise jusque la? (force donc cycles boucle pas? ou  ajuster a 1 moment?)
        vAmount = [0.1] * nb_res  # [20.]# T Fresh Weight.ha-1 (equivalent QRES)
        Vprop1 = [1. / 3., 1. / 3., 1. / 3.] + 50 * [0.]  # distribution dans les horizons #-> change 27 en 50 pour etre sur d'avoir le nb d'horizon-> a adapter selon le vrai nbr d'horizons!!!
        vProps = [Vprop1] * nb_res  # [Vprop1]#[Vprop1, Vprop1, Vprop1]

        # S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

        print('soil init', CNRES, vAmount, vProps, WC, CC)
        S.init_residues(CNRES, vAmount, vProps, WC, CC)

        print('ls_CRES', np.shape(S.ls_CRES))
        print('ls_CBio', np.shape(S.ls_CBio))
        print('parResi', S.parResi)

        return CC



#def init_ParamP(path_plante, ongletP, ongletPvois, nbcote, deltalevmoy, deltalevsd, seed_=0, type='homogeneous', opt=4, ongletScenar1='default', ongletScenar2='default', idscenar1=1, idscenar2=1, mn_sc=None):
def init_ParamP_VGL(path_plante, ongletP, ongletPvois, nbcote, deltalevmoy, deltalevsd, Plt_seed, seed_=0, type='homogeneous', opt=4, ongletScenar1='default', ongletScenar2='default', idscenar1=1, idscenar2=1, mn_sc=None, opt_sd=0, opt_covar=0, path_variance_geno=None, path_variance_matrix=None, idscenar1_sd=None, idscenar2_sd=None):
    """ """

    # 1) cree liste des paramtres plante (1 dico par plante)
    # nbcote = nombre de plante sur un cote en supposant repartition homogene
    g4 = IOxls.read_plant_param(path_plante, ongletP)
    g4 = IOxls.modif_param(g4, ongletP, ongletScenar1, idscenar1, mn_sc=mn_sc)
    g5 = IOxls.read_plant_param(path_plante, ongletPvois)
    g5 = IOxls.modif_param(g5, ongletPvois, ongletScenar2, idscenar2, mn_sc=mn_sc)
    if type == 'homogeneous':  # cas d'un couvert monospe homogene
        ParamP = [g4] * nbcote * nbcote
    elif type == 'damier8' or type == 'random8' or type == 'damier8_sp1' or type == 'damier8_sp2':  # damier binaire 64 plantes
        if nbcote == 8:
            ParamP = sh.damier8(g4, g5, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 64 plant design')
    elif type == 'damier16' or type == 'random16' or type == 'damier16_sp1' or type == 'damier16_sp2':  # damier binaire 256 plantes
        if nbcote == 16:
            ParamP = sh.damier16(g4, g5, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 256 plant design')
    elif type == 'row4' or 'row4_sp1' or 'row4_sp2':  # 4 rangs - 500pl.m-2
        ParamP, cart_ = sh.row4(g4, g5, nbprow=nbcote, opt=opt)
    else:  # defautl= force nb plante comme nbcote
        ParamP = [g4] * nbcote
        # ParamP = [g4]*10#*2#*30#[g6, g4, g6, g4, g6, g4, g6]#[g6,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4]#[g4]
        # [g1, g2, g3]#[g4, g5, g5, g5, g5, g5, g5]

    # 2) ajout variabilite sd si opt_sd==1 (variabilite intra) ; possible seulement si pas analyse de sensibilite (onglet scenar=default)
    # test pour esp 1, Len avec sd=0.5
    # ls_sdpar = [0.5] #ecart type parametre - a passer via un fichier d'entree comme scenar? autrement (multivarie ou directement dans fichier parametre plante?)
    # ls_parname = ['Len'] #liste a recuperer via un fichier d'entree

    if opt_sd == 2:  # new version: lecture des CV et calcul des SD a partir des parametres moyens

        sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
        ls_parname_g4 = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        CVs4 = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar1_sd][ls_parname_g4][0:]
        CVs4 = CVs4.iloc[0]
        # calcul des sigma a partir des cv et moy
        moyP4 = []
        for p in ls_parname_g4:
            moyP4.append(g4[p])

        ls_sdpar_g4 = abs(CVs4 * np.array(moyP4))

        sd_fichier_g5 = pd.read_excel(path_variance_geno, sheet_name=ongletPvois)
        ls_parname_g5 = list(sd_fichier_g5.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        CVs5 = sd_fichier_g5.loc[sd_fichier_g5["id_scenario"] == idscenar2_sd][ls_parname_g5][0:]
        CVs5 = CVs5.iloc[0]
        # calcul des sigma a partir des cv et moy
        moyP5 = []
        for p in ls_parname_g5:
            moyP5.append(g5[p])

        ls_sdpar_g5 = abs(CVs5 * np.array(moyP5))

        # print(CVs5, CVs4)
        # print(ls_sdpar_g5, ls_sdpar_g4)
        # print(ls_parname_g5, ls_parname_g4)

        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname=  ls_parname_g4, ls_sdpar=ls_sdpar_g4)
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname=  ls_parname_g5, ls_sdpar=ls_sdpar_g5)
        # print(df1, df2)

        if opt_covar == 1:
            # lecture de deux matrices de correlation -> bon onlet / scenario (id colonne dedie) dans bon fichier "corr_Matrix" -> a faire!
            # idscenar1_sd et idscenar2_sd: meme id que pour opt_sd!
            covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
            cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar1_sd]["correlation"])
            size_mat4 = int(sqrt(len(cor_g4)))
            cor_mat4 = cor_g4.reshape((size_mat4, size_mat4))

            covar_fichier_g5 = pd.read_excel(path_variance_matrix, sheet_name=ongletPvois)
            cor_g5 = np.array(covar_fichier_g5[covar_fichier_g5['id_scenario'] == idscenar2_sd]["correlation"])
            size_mat5 = int(sqrt(len(cor_g5)))
            cor_mat5 = cor_g5.reshape((size_mat5, size_mat5))

            # print('covar g4', cor_mat4, cor_mat5)
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=cor_mat4)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=cor_mat5)

        else:
            # defaut = no matrix of covariance = independant
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=None)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=None)

    elif opt_sd == 1:  # remet old file = lecture directe des valeurs de SD independament des valeurs moyennes

        sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
        ls_parname_g4 = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        ls_sdpar_g4 = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar1_sd][ls_parname_g4][0:]
        ls_sdpar_g4 = ls_sdpar_g4.iloc[0]

        sd_fichier_g5 = pd.read_excel(path_variance_geno, sheet_name=ongletPvois)
        ls_parname_g5 = list(sd_fichier_g5.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        ls_sdpar_g5 = sd_fichier_g5.loc[sd_fichier_g5["id_scenario"] == idscenar2_sd][ls_parname_g5][0:]
        ls_sdpar_g5 = ls_sdpar_g5.iloc[0]

        if opt_covar == 1:
            covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
            cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar1_sd]["correlation"])
            size_mat4 = int(sqrt(len(cor_g4)))
            cor_mat4 = cor_g4.reshape((size_mat4, size_mat4))

            covar_fichier_g5 = pd.read_excel(path_variance_matrix, sheet_name=ongletPvois)
            cor_g5 = np.array(covar_fichier_g5[covar_fichier_g5['id_scenario'] == idscenar2_sd]["correlation"])
            size_mat5 = int(sqrt(len(cor_g5)))
            cor_mat5 = cor_g5.reshape((size_mat5, size_mat5))

            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=cor_mat4)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=cor_mat5)

        else:
            # defaut = no matrix of covariance = independant
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=None)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=None)

    # 3) ajout de parametre 'recalcule'
    # roots
    for nump in range(len(ParamP)):
        # update des parametre racinaire
        rt.update_root_params(ParamP[nump])  # 'lsDrac', 'nb_ordre_rac', 'lsVrac', 'lsDemanDRac', 'LDs'
        # print 'LDs', ParamP[nump]['LDs2'], ParamP[nump]['LDs3'], ParamP[nump]['lsDrac'], ParamP[nump]['GDs2'], ParamP[nump]['GDs3']

        ParamP[nump]['profilRoot'] = rt.rootTropism(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.)
        # ParamP[nump]['profilRoot'] = rt.rootTropism_df(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.) #plus lent!

    # shoot profiles
    ParamP = sh.update_shoot_params(ParamP)

    # 4) random number generators per plant and test_retard vector
    # seed_=1
    # creation d'une liste d'objet random pour chaque plante (tirage random par plante et plus pour toutes les plantes)
    ls_seeds = []
    for nump in range(len(ParamP)):
        ls_seeds.append(Plt_seed(nump + seed_))
        # ls_seeds.append(Plt_seed(seed_))

    test_retard = []
    for nump in range(len(ParamP)):
        test_retard.append(max(0, ls_seeds[nump].randgen.normal(deltalevmoy, deltalevsd)))
        # max(0,random.gauss(deltalevmoy,deltalevsd))

    # !! ls_seeds pas passe pour tirage multivarie opt_sd? -> non, tirage multivarie avec Rssed pour ttes les plantes ensemble

    # 5) reduit a une espece si veut simul separee
    lsidP = range(len(ParamP))  # default value = all
    if type == 'damier8_sp1' or type == 'damier16_sp1' or type == 'row4_sp1':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ongletP)
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)
    elif type == 'damier8_sp2' or type == 'damier16_sp2' or type == 'row4_sp2':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ongletPvois)
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)

    # print('test_retard', test_retard)

    nbplantes = len(ParamP)
    return ParamP, nbplantes, ls_seeds, lsidP, test_retard
