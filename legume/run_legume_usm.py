
#import the modules necessary to initiate the L-systems
from openalea.lpy import *

import os
import sys

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'

print(('path', path_))

sys.path.insert(0, path_)
import IOxls
import IOtable

import getopt





def lsystemInputOutput_usm(path_, fxls_usm, i=0, foldin = 'input', ongletBatch = 'exemple'):
    """" cree et update l-system en fonction du fichier usm """

    # lecture de la liste des usm
    # path_ = r'H:\devel\grassland\grassland\L-gume'
    usm_path = os.path.join(path_, foldin, fxls_usm)
    usms = IOxls.xlrd.open_workbook(usm_path)
    ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))
    #foldin = pour cas ou fichier d'usm dans un sous dossier different de input / tous les autres sont dans input

    #nom fichier en dur (pas en entree de la fonction) + onglet determine par geno
    fscenar = 'liste_scenarios.xls' #'liste_scenarios_exemple.xls'
    fsd = 'exemple_sd.xls' #nom mis a jour mais pas table variance_geno


    testsim = {} #dico sortie avec un nom d'usm
    name = str(int(ls_usms['ID_usm'][i])) + '_' + str(ls_usms['l_system'][i])[0:-4]
    seednb = int(ls_usms['seed'][i])
    #names.append(name)
    path_lsys = os.path.join(path_, str(ls_usms['l_system'][i]))
    testsim[name] = Lsystem(path_lsys)  # objet l-system

    # testsim[name].ongletM = str(ls_usms['ongletM'][i])
    meteo_path_ = os.path.join(path_, 'input', str(ls_usms['meteo'][i]))
    ongletM_ = str(ls_usms['ongletM'][i])
    testsim[name].meteo = IOxls.read_met_file(meteo_path_, ongletM_)

    # testsim[name].ongletMn = str(ls_usms['ongletMn'][i])
    mn_path_ = os.path.join(path_, 'input', str(ls_usms['mng'][i]))
    ongletMn_ = str(ls_usms['ongletMn'][i])
    testsim[name].mng = IOxls.read_met_file(mn_path_, ongletMn_)

    ini_path_ = os.path.join(path_, 'input', str(ls_usms['inis'][i]))
    ongletIni_ = str(ls_usms['ongletIn'][i])
    testsim[name].inis = IOxls.read_plant_param(ini_path_, ongletIni_)

    # testsim[name].ongletP = str(ls_usms['ongletP'][i])
    path_plante = os.path.join(path_, 'input', str(ls_usms['plante'][i]))
    ongletP = str(ls_usms['ongletP'][i])
    ongletPvois = str(ls_usms['ongletVoisin'][i])
    testsim[name].path_plante = path_plante

    path_scenar = os.path.join(path_, 'input', fscenar)
    testsim[name].mn_sc = path_scenar

    path_variance_geno = os.path.join(path_, 'input', fsd)
    testsim[name].path_variance_geno = path_variance_geno
    # la, lire scenario et changer parametres
    idscenar1 = int(ls_usms['scenario1'][i])
    idscenar2 = int(ls_usms['scenario2'][i])
    idscenar1_sd = int(ls_usms['scenario1_sd'][i])
    idscenar2_sd = int(ls_usms['scenario2_sd'][i])
    ongletScenar2 = ongletPvois  # fait porter les changements sur fichier parametre voisin
    ongletScenar1 = ongletP

    # sol
    path_sol = os.path.join(path_, 'input', str(ls_usms['sol'][i]))
    ongletS = str(ls_usms['ongletS'][i])
    par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)
    par_SN['concrr'] = 0.  # force eau de pluie dans ls test (a retirer)
    # testsim[name].ongletS = str(ls_usms['ongletS'][i])
    testsim[name].par_SN = par_SN
    testsim[name].par_sol = par_sol

    # nbcote=7 # a passer ext
    # testsim[name].ParamP = [g4]*int(ls_usms['nbplt'][i])#nbcote*nbcote
    optdamier = int(ls_usms['damier'][i])
    nbcote = int(ls_usms['nbcote'][i])
    ### testsim[name].ParamP = damier8(g4,g5,opt=optdamier)
    if str(ls_usms['arrangement'][i]) == 'damier8':
        arrang = 'damier' + str(optdamier)
    elif str(ls_usms['arrangement'][i]) == 'row4':
        arrang = 'row' + str(optdamier)
    else:
        arrang = str(ls_usms['arrangement'][i]) + str(optdamier)

    nommix = '_' + ongletP + '-' + ongletPvois + '_' + arrang + '_scenario' + str(idscenar2) + '-' + str(idscenar1)

    testsim[name].ongletP = ongletP
    testsim[name].ongletPvois = ongletPvois
    testsim[name].nbcote = nbcote
    testsim[name].opt_sd = int(ls_usms['opt_sd'][i])  # 1
    testsim[name].cote = float(ls_usms['cote'][i])
    testsim[name].deltalevmoy = float(ls_usms['retard'][i])
    testsim[name].deltalevsd = float(ls_usms['sd_retard'][i])
    testsim[name].typearrangement = str(ls_usms['arrangement'][i])
    testsim[name].optdamier = optdamier
    testsim[name].idscenar1 = idscenar1
    testsim[name].idscenar2 = idscenar2
    testsim[name].ongletScenar2 = ongletScenar2
    testsim[name].ongletScenar1 = ongletScenar1
    testsim[name].idscenar1_sd = idscenar1_sd
    testsim[name].idscenar2_sd = idscenar2_sd
    testsim[name].Rseed = seednb
    testsim[name].DOYdeb = int(ls_usms['DOYdeb'][i])
    testsim[name].DOYend = int(ls_usms['DOYend'][i])

    # mise a jour derivartionLength & axiom
    testsim[name].derivationLength = int(ls_usms['DOYend'][i]) - int(ls_usms['DOYdeb'][i])  # derivationLength variable predefinie dans L-py
    if str(ls_usms['arrangement'][i]) == 'row4':  # carre rang heterogene
        nbplantes = nbcote * 4
    else:  # carre homogene
        nbplantes = nbcote * nbcote

    a = AxialTree()
    a.append(testsim[name].attente(1))
    for j in range(0, nbplantes):
        a.append(testsim[name].Sd(j))

    testsim[name].axiom = a  # passe un axial tree, pas de chaine de caractere

    if int(ls_usms['opt_sd'][i]) == 1:
        sdname = '_SD' + str(idscenar2_sd) + '-' + str(idscenar1_sd)
    else:
        sdname = '_-'

    # path fichiers de sortie
    testsim[name].path_out = os.path.join(path_, str(ls_usms['folder_out'][i]))
    testsim[name].outvarfile = 'toto_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].lsorgfile = 'lsAxes_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outHRfile = 'outHR_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].resrootfile = 'resroot_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outBilanNfile = 'BilanN_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outimagefile = 'scene_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.bmp'  # 'scene.bmp'
    testsim[name].outsdfile = 'paramSD_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + '_' + sdname + '_' + '.csv'

    # plante si dossier out pas cree
    # pourrait faire la lecture les ls_usm directement dans le l-system pour faciliter...+
    return testsim #dico avec {nom:lsystem}


def runlsystem(lsys, name):
    """ run the lsystem from a dict with its name"""
    try:
        lsys[name].derive()
        lsys[name].clear()
        print((''.join((name, " - done"))))
    except Exception as e:
        print(e)


def animatelsystem(lsys, name):
    lsys[name].animate()
    lsys[name].clear()
    print((''.join((name," - done"))))



if __name__ == '__main__':

    path_input = path_
    usm_file = 'liste_usms_exemple.xls' #ex fxls
    IDusm = 0

    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:s:d", ["usm_file=", "usm_scenario="])#cf l-grass
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-f", "--file"):
            usm_file = arg
        elif opt in ("-u", "--usm"):
            IDusm = int(arg)
        #elif opt in ("-o", "--outputs"):
        #    outputs = arg
        # pour le moment output folde fourni dans usm -> a changer

    mylsys = lsystemInputOutput_usm(path_input, usm_file, IDusm, foldin='multisim', ongletBatch='exemple')
    keyname = list(mylsys)[0]
    runlsystem(mylsys, keyname)
    #animatelsystem(mylsys, keyname)


#finir rendre accessible en externe le fichier de gestion des sorties!
# -> input and output folders
#rendre output accessible hors usm?
#peut appeler en ligne de commande "python run_l-egume_usm.py -s 1"


#a tester en remplacement dans le batch!
#prevoir de tout mettre dans le meme dossier en entree


