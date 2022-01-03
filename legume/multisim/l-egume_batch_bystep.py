# batch pour L-egume pour simul step-by-step and external coupling with environment 19/01/19


# import the modules necessary to initiate the L-systems
from openalea.lpy import *
import multiprocessing

import os
import sys

try:
    import legume

    path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'  # r'C:\devel\grassland'

print(('path', path_))

sys.path.insert(0, path_)
import IOxls
import IOtable
import run_legume_usm as runl
import ShootMorpho as sh
import daily_loop as loop
import numpy as np

try:
    from .riri5 import RIRI5 as riri #import de la version develop si module soil3ds est installe
except:
    import RIRI5 as riri

global foldin, fxls, ongletBatch
# to define if used for multisimulation or non-regression tests
opttest = 'exemple'#'mayssa'#'sdBea'#'OATbea'#1#2#'Histor'#1#4 ##2#1#  5#4#2#'autre'#0#13#'exemple_BA'#
if opttest == 1 or opttest == 2 or opttest == 3 or opttest == 4 or opttest == 5:  # si multisim des test de non regression (1 or 2)
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'test\inputs'
    fxls = 'liste_usms_nonregression.xls'
    if opttest == 1:  # non regression
        ongletBatch = 'test'
        foldout = os.path.join(path_, 'test\lastcheck')
    elif opttest == 2:  # obssim
        ongletBatch = 'valid'
        foldout = os.path.join(path_, 'test\lastvalidBis')
    elif opttest == 3:  # solnu
        ongletBatch = 'solnu'  #
        foldout = os.path.join(path_, 'output')
    elif opttest == 4:  # champ
        ongletBatch = 'test_champ'  #
        foldout = os.path.join(path_, 'test', 'test_champ')
    elif opttest == 5:  # pour bea
        ongletBatch = 'test_beajul'  #
        foldout = os.path.join(path_, 'test', 'test2')
elif opttest == 'exemple':
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'multisim'
    foldout = os.path.join(path_, 'output')
    fxls = 'liste_usms_exemple.xls'
    ongletBatch = 'exemple'
    #fscenar = 'liste_scenarios_exemple.xls'
elif opttest == 'exemple_BA':
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'multisim'
    fxls = 'liste_usms_exemple_BA.xls'
    ongletBatch = 'exemple'
    foldout = os.path.join(path_, 'test', 'test2')
elif opttest == 'mayssa':
    # global foldin, fxls, ongletBatch, fscenar
    # to be manually updated
    foldin =  'C:\inputs\inputs mayssa\DIGITLUZ'#r'C:\inputs\inputs test variance BLW'#os.path.join(path_, 'input')#'input'  # 'multisim'
    fxls = 'liste_usms_eval.xls'#'liste_usms_exemple.xls'#'liste_usms_essais.xls'  # 'liste_usms_mix.xls'
    ongletBatch = 'valid'#'Param1GL'#'OATbea'#'Histor'#'Champs'  # 'SimTest'#
    foldout =  'C:\inputs\inputs mayssa\output'#os.path.join(path_, 'output')
else:  # to personalize - other multisimulation to be defined (0)
    # global foldin, fxls, ongletBatch, fscenar
    # to be manually updated
    foldin =  r'C:\inputs\inputs test variance BLW'#os.path.join(path_, 'input')#'input'  # 'multisim'
    fxls = 'liste_usms_exemple.xls'#'liste_usms_essais.xls'  # 'liste_usms_mix.xls'
    ongletBatch = 'Param1GL'#'OATbea'#'Histor'#'Champs'  # 'SimTest'#
    foldout = os.path.join(path_, 'output')

#scenar et fsd definis en dur dans la fonction lsystemInputOutput_usm
#fscenar = 'liste_scenarios.xls'
#fsd = 'exemple_sd.xls'



#lecture de la liste des usm
#path_ = r'H:\devel\grassland\grassland\L-gume'
mn_path = os.path.join(path_, foldin, fxls)#(path_,'test\inputs','liste_usms_nonregression.xls')#(path_,'multisim','liste_usms_mix.xls')#(path_,'liste_usms_mix_these lucas.xls')#
#ongletBatch = 'test'#'SimTest'#'complement'#'Feuil1'#'Sensi'#
usms = IOxls.xlrd.open_workbook(mn_path)
ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))



# cree la liste de L-systems et liste des noms
testsim = {}
names = []
for i in range(len(ls_usms['ID_usm'])):
    if int(ls_usms['torun'][i]) == 1:  # si 1 dans la colonne 'torun' l'ajoute a la liste
        #mylsys = runl.lsystemInputOutput_usm(path_, fxls, i, foldin=foldin, ongletBatch=ongletBatch)
        mylsys = runl.lsystemInputOutput_usm(fxls, foldin=foldin, ongletBatch=ongletBatch, i=i, path_OUT=foldout)
        name = list(mylsys)[0]
        names.append(name)
        testsim[name] = mylsys[name]


nb_usms = len(names)  # len(ls_usms['ID_usm'])#len(names)#

# print (nb_usms, names)


# print nb_usms, names


# function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    testsim[names[n]].derive()
    testsim[names[n]].clear()
    print((''.join((names[n], " - done"))))


def runlsystem_bystep(n):
    """run du niem l-system by step dans une liste (names)"""

    lsys = testsim[names[n]]
    lstring = lsys.axiom
    nb_iter = lsys.derivationLength #lire dans derivation_length #335 #30
    lsys.opt_external_coupling = 1 # met a un l'option external coupling

    #option stress a zero si besoin
    #lsys.opt_stressN = 0
    #lsys.opt_stressW = 0

    for i in range(nb_iter+1):
        print('iter ',i,n)
        lstring = lsys.derive(lstring, i, 1)

        ## daily loop
        tag_loop_inputs = lsys.tag_loop_inputs
        invar, outvar, invar_sc, ParamP, station, carto, meteo_j, mng_j, DOY, cutNB, start_time, nbplantes, surfsolref, m_lais, dicFeuilBilanR, surf_refVOX, triplets, ls_dif, S, par_SN, lims_sol, ls_roots, stateEV, Uval, b_, ls_mat_res, vCC, ls_ftswStress, ls_NNIStress, ls_TStress, lsApex, lsApexAll, dicOrgans, deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant, I_I0profilInPlant, NlClasses, NaClasses, NlinClasses, opt_stressW, opt_stressN, opt_stressGel, opt_residu = tag_loop_inputs

        ############
        # step light transfer coupling
        ############

        # PAR / Blue voxel
        tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, meteo_j['I0'] * surf_refVOX]  # input tag

        # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
        local_res_trans, local_res_abs_i = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])

        res_trans, res_abs_i = local_res_trans, local_res_abs_i  # mise a jour variables globales

        # R_FR voxel (calcul de zeta)
        tag_light_inputs2 = [res_trans / (meteo_j['I0'] * surf_refVOX)]  # input tag
        local_res_rfr = riri.rfr_calc_relatif(*tag_light_inputs2)  # (res_trans/(meteo_j['I0']*surf_refVOX))

        res_rfr = local_res_rfr  # mise a jour variables globales

        # calul des interception feuille et ls_epsi plante
        dicFeuilBilanR = sh.calc_paraF(dicFeuilBilanR, m_lais, res_abs_i)
        ls_epsi, invar = loop.step_epsi(invar, res_trans, dicFeuilBilanR, meteo_j, surfsolref)

        ##########
        # Step Potential plant growth
        ##########

        invar, outvar, ls_demandeN_bis, temps = loop.daily_growth_loop(ParamP, invar, outvar, ls_epsi, meteo_j, mng_j,
                                                                       nbplantes, surfsolref, ls_ftswStress,
                                                                       ls_NNIStress, ls_TStress, lsApex, lsApexAll,
                                                                       opt_stressW, opt_stressN, opt_stressGel)

        ##########
        # step soil
        ##########

        tag_inputs_soil_step = [S, par_SN, surfsolref, stateEV, Uval, b_, meteo_j, mng_j, ParamP, ls_epsi, ls_roots, ls_demandeN_bis, opt_residu] # input tag

        res_soil_step = loop.step_bilanWN_sol(*tag_inputs_soil_step)
        S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = res_soil_step  # unpacks results from a list and updates global variables

        ##########
        # setp update plant stress variables
        ##########

        tag_inputs_stress = [ParamP, invar, invar_sc, temps, DOY, nbplantes, surfsolref, ls_epsi, ls_ftsw, ls_transp,
                             ls_Act_Nuptake_plt, ls_demandeN_bis, ls_ftswStress, ls_TStress, dicOrgans, dicFeuilBilanR, lsApex,
                             start_time, cutNB, deltaI_I0, nbI_I0, I_I0profilLfPlant, I_I0profilPetPlant,
                             I_I0profilInPlant, NlClasses, NaClasses, NlinClasses, outvar]

        invar, invar_sc, outvar, I_I0profilInPlant, ls_ftswStress, ls_NNIStress, ls_TStress = loop.Update_stress_loop(*tag_inputs_stress)

        ##########
        # step update soil residues senescence
        ##########
        tag_inputs_residue_updt = [ls_mat_res, vCC, S, ls_roots, par_SN['PROFHUMs'], ParamP, invar, opt_residu,opt_stressGel] # input tag

        res_residue_step = loop.update_residue_mat(*tag_inputs_residue_updt)
        ls_mat_res, S = res_residue_step  # unpacks results from a list and updates global variables


        #########
        # reinjecte les sorties midiee dans le lsystem
        #########
        lsys.invar = invar
        lsys.outvar = outvar
        lsys.invar_sc = invar_sc

        lsys.S = S
        lsys.stateEV = stateEV

        lsys.res_trans = res_trans
        lsys.res_abs_i = res_abs_i
        lsys.res_rfr = res_rfr

        lsys.ls_ftswStress = ls_ftswStress
        lsys.ls_NNIStress = ls_NNIStress
        lsys.ls_TStress = ls_TStress
        lsys.I_I0profilInPlant = I_I0profilInPlant

        #ls_mat_res bien pris en compte??? -> yes

        #yes!! + produit bien la meme sortie!! que fichier sans by_pass
        #par contre temps de calcul plus long!!!  525s au lieu de 323!!!

        #pour sauvegarder geometrie
        #s_leg = lsys.sceneInterpretation(lstring)
        #s_leg.save("s_leg.bgeom")

    testsim[names[n]].clear()
    print((''.join((names[n], " - done"))))

#check / debug partie N / NNI a la coupe??




def runl2system_bystep(n, m):
    """run de deux l-system by step dans une liste (names) - couplage lumiere et sol"""
    # pour le moment: run en // pas couple! A tester avec usm 17111 et 17112 dans exemple

    #lecture des deux l-systems
    lsys1 = testsim[names[n]]
    lstring1 = lsys1.axiom
    nb_iter1 = lsys1.derivationLength #lire dans derivation_length #335 #30
    lsys1.opt_external_coupling = 1 # met a un l'option external coupling

    lsys2 = testsim[names[m]]
    lstring2 = lsys2.axiom
    nb_iter2 = lsys2.derivationLength  # lire dans derivation_length #335 #30
    lsys2.opt_external_coupling = 1  # met a un l'option external coupling

    for i in range(nb_iter1):
        print('iter ',i,n,m)
        lstring1 = lsys1.derive(lstring1, i, 1)
        lstring2 = lsys2.derive(lstring2, i, 1)

        ## daily loop
        # bien recuperer les ttes les variables des 2 l-systems
        # mettre en commun les variables qui reunissent les 2
        # pb sur les invar pour reunir??! -> OK pour lumiere, mais pb pour le sol!!! -> faire passer des sorties et MAJ invar en dehors de loop.step_bilanWN_sol et loop.update_residue_mat et loop.step_epsi

        tag_loop_inputs1 = lsys1.tag_loop_inputs
        invar1, outvar1, invar_sc1, ParamP1, station1, carto1, meteo_j1, mng_j1, DOY1, cutNB1, start_time1, nbplantes1, surfsolref1, m_lais1, dicFeuilBilanR1, surf_refVOX1, triplets1, ls_dif1, S1, par_SN1, lims_sol1, ls_roots1, stateEV1, Uval1, b_1, ls_mat_res1, vCC1, ls_ftswStress1, ls_NNIStress1, ls_TStress1, lsApex1, lsApexAll1, dicOrgans1, deltaI_I01, nbI_I01, I_I0profilLfPlant1, I_I0profilPetPlant1, I_I0profilInPlant1, NlClasses1, NaClasses1, NlinClasses1, opt_stressW1, opt_stressN1, opt_stressGel1, opt_residu1 = tag_loop_inputs1

        tag_loop_inputs2 = lsys2.tag_loop_inputs
        invar2, outvar2, invar_sc2, ParamP2, station2, carto2, meteo_j2, mng_j2, DOY2, cutNB2, start_time2, nbplantes2, surfsolref2, m_lais2, dicFeuilBilanR2, surf_refVOX2, triplets2, ls_dif2, S2, par_SN2, lims_sol2, ls_roots2, stateEV2, Uval2, b_2, ls_mat_res2, vCC2, ls_ftswStress2, ls_NNIStress2, ls_TStress2, lsApex2, lsApexAll2, dicOrgans2, deltaI_I02, nbI_I02, I_I0profilLfPlant2, I_I0profilPetPlant2, I_I0profilInPlant2, NlClasses2, NaClasses2, NlinClasses2, opt_stressW2, opt_stressN2, opt_stressGel2, opt_residu2 = tag_loop_inputs2

        #TO DO: definir variables communes dans run couple: sol S de S1 et variables reunissant facilement les 2
        # sortir les invar plantes du sol! (simplifiant sans doute les entree/sorties)
        # pour le moment run //

        ############
        # step light transfer coupling
        ############

        # PAR / Blue voxel
        tag_light_inputs1 = [m_lais1 / surf_refVOX1, triplets1, ls_dif1, meteo_j1['I0'] * surf_refVOX1]  # input tag
        tag_light_inputs2 = [m_lais2 / surf_refVOX2, triplets2, ls_dif2, meteo_j2['I0'] * surf_refVOX2]  # input tag

        # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
        local_res_trans1, local_res_abs_i1 = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs1, optsky=station1['optsky'], opt=station1['sky'])
        local_res_trans2, local_res_abs_i2 = riri.calc_extinc_allray_multi_reduced(*tag_light_inputs2, optsky=station2['optsky'], opt=station2['sky'])

        res_trans1, res_abs_i1 = local_res_trans1, local_res_abs_i1  # mise a jour variables globales
        res_trans2, res_abs_i2 = local_res_trans2, local_res_abs_i2  # mise a jour variables globales

        # R_FR voxel (calcul de zeta)
        tag_light_inputs2_1 = [res_trans1 / (meteo_j1['I0'] * surf_refVOX1)]  # input tag
        tag_light_inputs2_2 = [res_trans2 / (meteo_j2['I0'] * surf_refVOX2)]  # input tag
        local_res_rfr1 = riri.rfr_calc_relatif(*tag_light_inputs2_1)
        local_res_rfr2 = riri.rfr_calc_relatif(*tag_light_inputs2_2)

        res_rfr1 = local_res_rfr1  # mise a jour variables globales
        res_rfr2 = local_res_rfr2  # mise a jour variables globales

        # ls_epsi
        # calul des interception feuille et ls_epsi plante
        dicFeuilBilanR1 = sh.calc_paraF(dicFeuilBilanR1, m_lais1, res_abs_i1)
        dicFeuilBilanR2 = sh.calc_paraF(dicFeuilBilanR2, m_lais2, res_abs_i2)

        ls_epsi1, invar1 = loop.step_epsi(invar1, res_trans1, dicFeuilBilanR1, meteo_j1, surfsolref1)
        ls_epsi2, invar2 = loop.step_epsi(invar2, res_trans2, dicFeuilBilanR2, meteo_j2, surfsolref2)


        ##########
        # Step Potential plant growth
        ##########

        invar1, outvar1, ls_demandeN_bis1, temps1 = loop.daily_growth_loop(ParamP1, invar1, outvar1, ls_epsi1, meteo_j1, mng_j1, nbplantes1, surfsolref1, ls_ftswStress1, ls_NNIStress1, ls_TStress1, lsApex1, lsApexAll1, opt_stressW1, opt_stressN1, opt_stressGel1)
        invar2, outvar2, ls_demandeN_bis2, temps2 = loop.daily_growth_loop(ParamP2, invar2, outvar2, ls_epsi2, meteo_j2, mng_j2, nbplantes2, surfsolref2, ls_ftswStress2, ls_NNIStress2, ls_TStress2, lsApex2, lsApexAll2, opt_stressW2, opt_stressN2, opt_stressGel2)

        ##########
        # step soil
        ##########

        tag_inputs_soil_step1 = [S1, par_SN1, surfsolref1, stateEV1, Uval1, b_1, meteo_j1, mng_j1, ParamP1, ls_epsi1, ls_roots1, ls_demandeN_bis1, opt_residu1]  # input tag
        tag_inputs_soil_step2 = [S2, par_SN2, surfsolref2, stateEV2, Uval2, b_2, meteo_j2, mng_j2, ParamP2, ls_epsi2, ls_roots2, ls_demandeN_bis2, opt_residu2]  # input tag


        res_soil_step1 = loop.step_bilanWN_sol(*tag_inputs_soil_step1)
        res_soil_step2 = loop.step_bilanWN_sol(*tag_inputs_soil_step2)
        S1, stateEV1, ls_ftsw1, ls_transp1, ls_Act_Nuptake_plt1, temps_sol1 = res_soil_step1  # unpacks results from a list and updates global variables
        S2, stateEV2, ls_ftsw2, ls_transp2, ls_Act_Nuptake_plt2, temps_sol2 = res_soil_step2  # unpacks results from a list and updates global variables

        ##########
        # setp update plant stress variables
        ##########

        tag_inputs_stress1 = [ParamP1, invar1, invar_sc1, temps1, DOY1, nbplantes1, surfsolref1, ls_epsi1, ls_ftsw1, ls_transp1,
                             ls_Act_Nuptake_plt1, ls_demandeN_bis1, ls_ftswStress1, ls_TStress1, dicOrgans1, dicFeuilBilanR1, lsApex1,
                             start_time1, cutNB1, deltaI_I01, nbI_I01, I_I0profilLfPlant1, I_I0profilPetPlant1,
                             I_I0profilInPlant1, NlClasses1, NaClasses1, NlinClasses1, outvar1]

        tag_inputs_stress2 = [ParamP2, invar2, invar_sc2, temps2, DOY2, nbplantes2, surfsolref2, ls_epsi2, ls_ftsw2, ls_transp2,
                              ls_Act_Nuptake_plt2, ls_demandeN_bis2, ls_ftswStress2, ls_TStress2, dicOrgans2, dicFeuilBilanR2, lsApex2,
                              start_time2, cutNB2, deltaI_I02, nbI_I02, I_I0profilLfPlant2, I_I0profilPetPlant2,
                              I_I0profilInPlant2, NlClasses2, NaClasses2, NlinClasses2, outvar2]

        invar1, invar_sc1, outvar1, I_I0profilInPlant1, ls_ftswStress1, ls_NNIStress1, ls_TStress1 = loop.Update_stress_loop(*tag_inputs_stress1)
        invar2, invar_sc2, outvar2, I_I0profilInPlant2, ls_ftswStress2, ls_NNIStress2, ls_TStress2 = loop.Update_stress_loop(*tag_inputs_stress2)

        ##########
        # step update soil residues senescence
        ##########
        tag_inputs_residue_updt1 = [ls_mat_res1, vCC1, S1, ls_roots1, par_SN1['PROFHUMs'], ParamP1, invar1, opt_residu1, opt_stressGel1]  # input tag
        tag_inputs_residue_updt2 = [ls_mat_res2, vCC2, S2, ls_roots2, par_SN2['PROFHUMs'], ParamP2, invar2, opt_residu2, opt_stressGel2]  # input tag

        res_residue_step1 = loop.update_residue_mat(*tag_inputs_residue_updt1)
        res_residue_step2 = loop.update_residue_mat(*tag_inputs_residue_updt2)
        ls_mat_res1, S1 = res_residue_step1  # unpacks results from a list and updates global variables
        ls_mat_res2, S2 = res_residue_step2  # unpacks results from a list and updates global variables

        #########
        # reinjecte les sorties dans le lsystem
        #########
        #lsystem 1
        lsys1.invar = invar1
        lsys1.outvar = outvar1
        lsys1.invar_sc = invar_sc1

        lsys1.S = S1
        lsys1.stateEV = stateEV1

        lsys1.res_trans = res_trans1
        lsys1.res_abs_i = res_abs_i1
        lsys1.res_rfr = res_rfr1

        lsys1.ls_ftswStress = ls_ftswStress1
        lsys1.ls_NNIStress = ls_NNIStress1
        lsys1.ls_TStress = ls_TStress1
        lsys1.I_I0profilInPlant = I_I0profilInPlant1

        # lsystem 2
        lsys2.invar = invar2
        lsys2.outvar = outvar2
        lsys2.invar_sc = invar_sc2

        lsys2.S = S2
        lsys2.stateEV = stateEV2

        lsys2.res_trans = res_trans2
        lsys2.res_abs_i = res_abs_i2
        lsys2.res_rfr = res_rfr2

        lsys2.ls_ftswStress = ls_ftswStress2
        lsys2.ls_NNIStress = ls_NNIStress2
        lsys2.ls_TStress = ls_TStress2
        lsys2.I_I0profilInPlant = I_I0profilInPlant2

    testsim[names[n]].clear()
    print((''.join((names[n], " - done"))))
    testsim[names[m]].clear()
    print((''.join((names[m], " - done"))))



def runl2systemLight_bystep(n, m):
    """run de deux l-system by step dans une liste (names) - couplage lumiere et sol"""
    # pour le moment: couplage light seulement avec RIRI5! A tester avec usm 17111 et 17112 dans exemple

    #lecture des deux l-systems
    lsys1 = testsim[names[n]]
    lstring1 = lsys1.axiom
    nb_iter1 = lsys1.derivationLength #lire dans derivation_length #335 #30
    lsys1.opt_external_coupling = 1 # met a un l'option external coupling

    lsys2 = testsim[names[m]]
    lstring2 = lsys2.axiom
    nb_iter2 = lsys2.derivationLength  # lire dans derivation_length #335 #30
    lsys2.opt_external_coupling = 1  # met a un l'option external coupling

    #force les option stress a zero pour n'intergir que pour lumiere
    lsys1.opt_stressN = 0
    lsys1.opt_stressW = 0
    lsys2.opt_stressN = 0
    lsys2.opt_stressW = 0

    for i in range(nb_iter1):
        print('iter ',i,n,m)
        lstring1 = lsys1.derive(lstring1, i, 1)
        lstring2 = lsys2.derive(lstring2, i, 1)

        ## daily loop
        # bien recuperer les ttes les variables des 2 l-systems
        # mettre en commun les variables qui reunissent les 2
        # pb sur les invar pour reunir??! -> OK pour lumiere, mais pb pour le sol!!! -> faire passer des sorties et MAJ invar en dehors de loop.step_bilanWN_sol et loop.update_residue_mat et loop.step_epsi

        tag_loop_inputs1 = lsys1.tag_loop_inputs
        invar1, outvar1, invar_sc1, ParamP1, station1, carto1, meteo_j1, mng_j1, DOY1, cutNB1, start_time1, nbplantes1, surfsolref1, m_lais1, dicFeuilBilanR1, surf_refVOX1, triplets1, ls_dif1, S1, par_SN1, lims_sol1, ls_roots1, stateEV1, Uval1, b_1, ls_mat_res1, vCC1, ls_ftswStress1, ls_NNIStress1, ls_TStress1, lsApex1, lsApexAll1, dicOrgans1, deltaI_I01, nbI_I01, I_I0profilLfPlant1, I_I0profilPetPlant1, I_I0profilInPlant1, NlClasses1, NaClasses1, NlinClasses1, opt_stressW1, opt_stressN1, opt_stressGel1, opt_residu1 = tag_loop_inputs1
        tag_loop_inputs2 = lsys2.tag_loop_inputs
        invar2, outvar2, invar_sc2, ParamP2, station2, carto2, meteo_j2, mng_j2, DOY2, cutNB2, start_time2, nbplantes2, surfsolref2, m_lais2, dicFeuilBilanR2, surf_refVOX2, triplets2, ls_dif2, S2, par_SN2, lims_sol2, ls_roots2, stateEV2, Uval2, b_2, ls_mat_res2, vCC2, ls_ftswStress2, ls_NNIStress2, ls_TStress2, lsApex2, lsApexAll2, dicOrgans2, deltaI_I02, nbI_I02, I_I0profilLfPlant2, I_I0profilPetPlant2, I_I0profilInPlant2, NlClasses2, NaClasses2, NlinClasses2, opt_stressW2, opt_stressN2, opt_stressGel2, opt_residu2 = tag_loop_inputs2

        #def variables communes
        meteo_j, station, surf_refVOX, triplets, surfsolref = meteo_j1, station1, surf_refVOX1, triplets1, surfsolref1

        ls_dif = ls_dif1
        m_lais = m_lais1 + m_lais2
        ## ??vraiment bon?? -> comme si une seule dist!!

        #gere difference de dsitib par especes
        #if (ls_dif1[0] == ls_dif2[0]).all(): # cas ou ce sont les memes distrib d'angles
        #    ls_dif = ls_dif1
        #    m_lais = m_lais1 + m_lais2
        #else: #distrib differents
        #    ls_dif = ls_dif1 + ls_dif2
        #    m_lais = np.array([m_lais1[0], m_lais2[0]])  # matrice 4D , 1ere dim = 1 matrice 3D par distrib angles espece


        #TO DO: definir variables communes dans run couple: sol S de S1 et variables reunissant facilement les 2
        # sortir les invar plantes du sol! (simplifiant sans doute les entree/sorties)

        ############
        # step light transfer coupling
        ############

        # PAR / Blue voxel
        tag_light_inputs = [m_lais / surf_refVOX, triplets, ls_dif, meteo_j['I0'] * surf_refVOX]  # input tag

        # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
        res_trans, res_abs_i =  riri.calc_extinc_allray_multi_reduced(*tag_light_inputs, optsky=station['optsky'], opt=station['sky'])


        # R_FR voxel (calcul de zeta)
        tag_light_inputs2 = [res_trans / (meteo_j['I0'] * surf_refVOX)]  # input tag
        res_rfr = riri.rfr_calc_relatif(*tag_light_inputs2)


        # ls_epsi et invar mis a jour par espece ls_epsi1, invar1 / ls_epsi2, invar2
        transmi_sol = np.sum(res_trans[-1][:][:]) / (meteo_j['I0'] * surfsolref)  # bon
        epsi = 1. - transmi_sol  # bon

        dicFeuilBilanR1 = sh.calc_paraF(dicFeuilBilanR1, m_lais, res_abs_i)
        dicFeuilBilanR2 = sh.calc_paraF(dicFeuilBilanR2, m_lais, res_abs_i)
        sh.calc_para_Plt(invar1, dicFeuilBilanR1)
        sh.calc_para_Plt(invar2, dicFeuilBilanR2)

        ls_epsi1 = epsi * invar1['parip'] / (np.sum(invar1['parip']) + np.sum(invar2['parip']) + 10e-15)
        ls_epsi2 = epsi * invar2['parip'] / (np.sum(invar1['parip']) + np.sum(invar2['parip']) + 10e-15)
        # mettre a jour epsi dans invar?

        print('espi', epsi, np.sum(ls_epsi1), np.sum(ls_epsi2))
        # espi ok - espi1 et 2 oscille bizarement?? option d'inclinaison des tiges ombrees??  bug?
        # sera a comparer avec des sim activant l'absence de stress edaphiques (compet lumiere uniquement)...

        ##########
        # Step Potential plant growth
        ##########

        invar1, outvar1, ls_demandeN_bis1, temps1 = loop.daily_growth_loop(ParamP1, invar1, outvar1, ls_epsi1, meteo_j1, mng_j1, nbplantes1, surfsolref1, ls_ftswStress1, ls_NNIStress1, ls_TStress1, lsApex1, lsApexAll1, opt_stressW1, opt_stressN1, opt_stressGel1)
        invar2, outvar2, ls_demandeN_bis2, temps2 = loop.daily_growth_loop(ParamP2, invar2, outvar2, ls_epsi2, meteo_j2, mng_j2, nbplantes2, surfsolref2, ls_ftswStress2, ls_NNIStress2, ls_TStress2, lsApex2, lsApexAll2, opt_stressW2, opt_stressN2, opt_stressGel2)

        ##########
        # step soil
        ##########

        tag_inputs_soil_step1 = [S1, par_SN1, surfsolref1, stateEV1, Uval1, b_1, meteo_j1, mng_j1, ParamP1, ls_epsi1, ls_roots1, ls_demandeN_bis1, opt_residu1]  # input tag
        tag_inputs_soil_step2 = [S2, par_SN2, surfsolref2, stateEV2, Uval2, b_2, meteo_j2, mng_j2, ParamP2, ls_epsi2, ls_roots2, ls_demandeN_bis2, opt_residu2]  # input tag

        res_soil_step1 = loop.step_bilanWN_sol(*tag_inputs_soil_step1)
        res_soil_step2 = loop.step_bilanWN_sol(*tag_inputs_soil_step2)
        S1, stateEV1, ls_ftsw1, ls_transp1, ls_Act_Nuptake_plt1, temps_sol1 = res_soil_step1  # unpacks results from a list and updates global variables
        S2, stateEV2, ls_ftsw2, ls_transp2, ls_Act_Nuptake_plt2, temps_sol2 = res_soil_step2  # unpacks results from a list and updates global variables

        ##########
        # setp update plant stress variables
        ##########

        tag_inputs_stress1 = [ParamP1, invar1, invar_sc1, temps1, DOY1, nbplantes1, surfsolref1, ls_epsi1, ls_ftsw1, ls_transp1,
                             ls_Act_Nuptake_plt1, ls_demandeN_bis1, ls_ftswStress1, ls_TStress1,dicOrgans1, dicFeuilBilanR1, lsApex1,
                             start_time1, cutNB1, deltaI_I01, nbI_I01, I_I0profilLfPlant1, I_I0profilPetPlant1,
                             I_I0profilInPlant1, NlClasses1, NaClasses1, NlinClasses1, outvar1]

        tag_inputs_stress2 = [ParamP2, invar2, invar_sc2, temps2, DOY2, nbplantes2, surfsolref2, ls_epsi2, ls_ftsw2, ls_transp2,
                              ls_Act_Nuptake_plt2, ls_demandeN_bis2, ls_ftswStress2, ls_TStress2, dicOrgans2, dicFeuilBilanR2, lsApex2,
                              start_time2, cutNB2, deltaI_I02, nbI_I02, I_I0profilLfPlant2, I_I0profilPetPlant2,
                              I_I0profilInPlant2, NlClasses2, NaClasses2, NlinClasses2, outvar2]

        invar1, invar_sc1, outvar1, I_I0profilInPlant1, ls_ftswStress1, ls_NNIStress1, ls_TStress1 = loop.Update_stress_loop(*tag_inputs_stress1)
        invar2, invar_sc2, outvar2, I_I0profilInPlant2, ls_ftswStress2, ls_NNIStress2, ls_TStress2 = loop.Update_stress_loop(*tag_inputs_stress2)

        ##########
        # step update soil residues senescence
        ##########
        tag_inputs_residue_updt1 = [ls_mat_res1, vCC1, S1, ls_roots1, par_SN1['PROFHUMs'], ParamP1, invar1, opt_residu1, opt_stressGel1]  # input tag
        tag_inputs_residue_updt2 = [ls_mat_res2, vCC2, S2, ls_roots2, par_SN2['PROFHUMs'], ParamP2, invar2, opt_residu2, opt_stressGel2]  # input tag

        res_residue_step1 = loop.update_residue_mat(*tag_inputs_residue_updt1)
        res_residue_step2 = loop.update_residue_mat(*tag_inputs_residue_updt2)
        ls_mat_res1, S1 = res_residue_step1  # unpacks results from a list and updates global variables
        ls_mat_res2, S2 = res_residue_step2  # unpacks results from a list and updates global variables

        #########
        # reinjecte les sorties dans le lsystem
        #########
        #lsystem 1
        lsys1.invar = invar1
        lsys1.outvar = outvar1
        lsys1.invar_sc = invar_sc1

        lsys1.S = S1
        lsys1.stateEV = stateEV1

        lsys1.res_trans = res_trans
        lsys1.res_abs_i = res_abs_i
        lsys1.res_rfr = res_rfr

        lsys1.ls_ftswStress = ls_ftswStress1
        lsys1.ls_NNIStress = ls_NNIStress1
        lsys1.ls_TStress = ls_TStress1
        lsys1.I_I0profilInPlant = I_I0profilInPlant1

        # lsystem 2
        lsys2.invar = invar2
        lsys2.outvar = outvar2
        lsys2.invar_sc = invar_sc2

        lsys2.S = S2
        lsys2.stateEV = stateEV2

        lsys2.res_trans = res_trans
        lsys2.res_abs_i = res_abs_i
        lsys2.res_rfr = res_rfr

        lsys2.ls_ftswStress = ls_ftswStress2
        lsys2.ls_NNIStress = ls_NNIStress2
        lsys2.ls_TStress = ls_TStress2
        lsys2.I_I0profilInPlant = I_I0profilInPlant2

    testsim[names[n]].clear()
    print((''.join((names[n], " - done"))))
    testsim[names[m]].clear()
    print((''.join((names[m], " - done"))))

#yes! mais gere pas les distrib angles differentes entre especes!
# sorties sont differents: !! gere pas les nump et tirages aleatoires / delaiinin -> a completer avec MAJ test_retard
# faire 1 test_retard par esp!!? -> fait pour tirage; OK pour les timing de depart
# mais tres vite pour la suite autres processus aleatoires donnent qd meme resultats differents pour la meme simuls
# introduit generateur de nb aleatoire par plante -> OK pour le premier mais apres derive toujours??
# utilise le generateur de numpy pour aussi capter tirage de loi binomiale ; tout avec meme seed -> ko
#desactive tirages loi binomiale racine? -> pas ; ko tjrs different


# to do:
# 1) modif light pour gerer esp avec differentes distributions! (fonction pour join ls_dif et m_lais)
# 2) test d'un couplage light et sol!!
# - sept sol n'a plus invar en entree! bien!!
# - mise en commun des ParamP, ls_epsi, ls_roots, ls_demandeN_bis
# - puis split apres le step des ls_ftsw1, ls_transp1, ls_Act_Nuptake_plt1, temps_sol1 pour continuer vers stress
# - gestion de une seule ls_mat_res avec differentes esp




#tuto: https://lpy-fb.readthedocs.io/en/latest/user/integration.html



# run the L-systems

### Ancienne methode ###
# if __name__ == '__main__':
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
    multiprocessing.freeze_support()
    CPUnb = multiprocessing.cpu_count() - 1  # nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
    print('nb CPU: ' + str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb)
    for i in range(1):#(int(nb_usms)):
        #pool.apply_async(runlsystem, args=(i,))   # Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
        runlsystem(i) #pour debug hors multisim (messages d'ereur visible)
        #runlsystem_bystep(i)
        #runl2system_bystep(i, i+1)
        #runl2systemLight_bystep(i, i+1)

    pool.close()
    pool.join()
