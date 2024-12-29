import xlrd
from copy import deepcopy
import random
import numpy as np
import pandas as pd
#from rpy_options import set_options
#set_options(RHOME='c:/progra~1/R/R-2.12.1')
#from rpy import r
import os
import xml.etree.ElementTree as ET # for XML files


def get_xls_col(sheet):
    """ recupere dans une feuille excel donnees par colone  """
    res=[]
    for i in range(sheet.ncols):
        res.append(sheet.col(i))

    # retient seulement les valeurs du dictionnaire
    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j]=res[i][j].value

    return res

def get_xls_row(sheet):
    """ recupere dans une feuille excel donnees par colone  """
    res=[]
    for i in range(sheet.nrows):
        res.append(sheet.row(i))

    # retient seulement les valeurs du dictionnaire
    for i in range(len(res)):
        for j in range(len(res[0])):
            res[i][j]=res[i][j].value

    return res



def t_list(tab):
    """transpose tab"""
    res = []
    for j in range(len(tab[0])):
        v = []
        for i in range(len(tab)):
            v.append(tab[i][j])
        
        res.append(v)

    return res

#def as_matrix(tab):
#    """ converts a list of list or a python array into an R matrix Robj """
#    r.rbind.local_mode(0)
#    r.c.local_mode(0)
#    x = r.c(tab[0])
#    for i in range (1,len(tab)):
#        x = r.rbind(x, r.c(tab[i]))
#
#    return x

def conv_dataframe(tab):
    """ converti liste de liste en dictionnaire; prend la cle comme le pemier element de la liste"""
    """ format a priori compatible pour conversion en data.frame R"""
    dat = {}
    for i in range(len(tab)):
        dat[str(tab[i][0])] = tab[i][1:]

    return dat #r.as_data_frame(dat)

def conv_list(tab):
    """ converti dictionnaireen liste de liste en ;  cle comme pemier element de la liste"""
    """ format compatible pour mes_csv"""
    dat = []
    for i in list(tab.keys()):
        v = [i]
        dat.append(v)

    count = 0
    for i in list(tab.keys()):
        for j in range(len(tab[i])):
            dat[count].append(tab[i][j])

        count = count+1

    return dat 

def extract_dataframe(dat, ls_cles, cle, val=None):
    """ extrait dans listes de cles ls_cles les lignes pour lesquelles cle=val; toutes si val=None """
    #cree liste d'index ou cle = val

    id = []
    for i in range(len(dat[cle])):
        if val == None:
            id.append(i)
        else:
            if dat[cle][i] == val:
                id.append(i)

    x = {}
    for k in ls_cles: # recupere les paires interessantes
        v = []
        for i in id: # les id respectant cle=val
            v.append(dat[k][i])

        x[k] = v

    return x
    #extract_dataframe(dat, cles, 'geno', geno)
    #extract_dataframe(dat, cles, 'geno')



def extract_list(dat, ls_id, ls_vals, L1=1):
    """ extrait avec un ET les lignes pour les quelles les colonnes numerotees ls_id prennent les valeurs ls_vals"""

    res = []
    for i in range(L1, len(dat)):
        bol = 1
        for j in range(len(ls_id)):
            if dat[i][ls_id[j]] == ls_vals[j]:
                bol = bol*1
            else :
                bol = bol*0

        if bol == 1:
            res.append(dat[i])

    return res


def read_plant_param(xls_path, onglet):
    """ lit l'onglet d'un fichier xls pour cree un dico de parametre plante (L-egume / VGL)- se base sur 3 colones 'name', 'nb_par' (0 si sclaire, n si vecteur), 'id_par'
    ; presuppose que valeurs a mettre dans un vecteur sont deja ordonnees"""
    book=xlrd.open_workbook(xls_path)
    shc = get_xls_col(book.sheet_by_name(onglet))
    dico_shc = conv_dataframe(shc)

    g = {}
    g['name']=onglet
    for i in range(len(dico_shc['name'])):
        if dico_shc['nb_par'][i] == 1: #si scalaire
            g[str(dico_shc['name'][i])]=dico_shc['value'][i]
        elif dico_shc['nb_par'][i] > 1 and dico_shc['id_par'][i] == 0:#si vecteur et sur la premier valeur
            paramlist=[]
            for j in range(int(dico_shc['nb_par'][i])):
                paramlist.append(dico_shc['value'][i+j])

            g[str(dico_shc['name'][i])]=paramlist   

    return g
    #path_source=r'H:\devel\grassland\grassland\L-gume\Parametres_plante.xls' 
    #onglet = 'geno_test'
    #g4 = read_plant_param(path_source, onglet)


def read_sol_param(xls_path, onglet):
    """ lecture du fichier sol VGL dans par_SN et mise en forme du dictionnaire par_sol"""
    par_SN = read_plant_param(xls_path, onglet)
    par_sol = {}
    for i in range(len(par_SN['soil_number'])):
        id = str(int(par_SN['soil_number'][i]))
        par_sol[id] = {'soil number': id, 'teta_sat': par_SN['teta_sat'][i], 'teta_fc': par_SN['teta_fc'][i], 'teta_wp': par_SN['teta_wp'][i], 'teta_ad': par_SN['teta_ad'][i],'DA':par_SN['DA'][i]}

    return par_SN, par_sol
    #idealement a faire: enlever par_sol et lire dnas le module sol directement dans le par_SN


def read_met_file(meteo_path, ongletM):
    """ lecture de fichier meteo / Mng type VGL """
    #meteo_path = os.path.join(path_leg,'meteo_exemple2.xls')#r'H:\devel\grassland\grassland\L-gume\meteo_exemple2.xls'
    #ongletM = 'Avignon30'#'exemple'#'morpholeg15'#'testJLD'#'competiluz'#
    met = xlrd.open_workbook(meteo_path)
    meteo = conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
    for k in ['month', 'day', 'DOY']: meteo[k] = list(map(int, meteo[k]))

    return meteo


def get_lsparami(ParamP, param):
    """ recupere une liste des parametre param de chaque plante de L-egume """
    v = []
    nbplt = len(ParamP)
    for i in range(nbplt):
        v.append(ParamP[i][param])
    
    return v

def modif_param(gx, ongletP, ongletScenar, idscenar, idlist=1, mn_sc=None):
    """ met a jour ParamP d'un genotype gx pour les variables et valeurs indiques dans ongletScnar """
    """ fait rien si ongletScenar='default' ou iscenar<0 ou onglet correspond pas a celui a modifier; sinon va modifier selon fichier scenario"""
    if ongletP == ongletScenar and idscenar > 0 and mn_sc != None:  # si onglet correspond a ongletscenat a modifier
        usc = xlrd.open_workbook(mn_sc)
        ls_sc = conv_dataframe(get_xls_col(usc.sheet_by_name(ongletScenar)))
        nb_modif = len(list(ls_sc.keys())) - 1  # nb de param a modifier
        ls_sc['id_scenario'] = list(map(int, ls_sc['id_scenario']))

        if nb_modif > 0:  # s'il y a des parametre a modifier
            idok = ls_sc['id_scenario'].index(idscenar)
            keys_modif = list(ls_sc.keys())
            keys_modif.remove('id_scenario')
            for k in keys_modif:
                if str(ls_sc[k][idok]) != '' or str(ls_sc[k][idok]) != 'NA':  # y a un valeur specifiee
                    if type(gx[k]) == type(0.):  # gere uniquement les parametres avec 1 seule valeur float (pas liste)
                        gx[k] = ls_sc[k][idok]
                    elif type(gx[k]) == list:  # pour les liste, gere uniquement un id de la liste (1 par defaut)
                        gx[k][idlist] == ls_sc[k][idok]
    return gx



def modif_ParamP_sd(ParamP, g4, ls_parname, ls_sdpar):
    """ modif paramP tirages 1 a 1 independents : loi normale monovariee """
    #ls_sdpar = [0.5]  # ecart type parametre - a passer via un fichier d'entree comme scenar? autrement (multivarie ou directement dans fichier parametre plante?)
    #ls_parname = ['Len']  # liste a recuperer via un fichier d'entree
    name1 = g4['name']
    for nump in range(len(ParamP)):
        if ParamP[nump]['name'] == name1:  # issu du bon onglet
            g = deepcopy(g4)
            for j in range(len(ls_parname)):
                parname = ls_parname[j]
                if parname in ['Largfeuille']: #liste param qui peuvent etre negatifs -> pas de contrainte de positivite
                    g[parname] = random.gauss(ParamP[nump][parname], ls_sdpar[j])
                else:
                    g[parname] = max(0.0000000001, random.gauss(ParamP[nump][parname], ls_sdpar[j]))  # seulement pour paramtere scalaie (pas liste) et positif (valeur >0)
                #!! ici a revoir car certain parametre pruvent etre negatifs + peut vouloir rirer dans differentes lois de distrib!

            ParamP[nump] = g

    return ParamP
    # modif_ParamP_sd(ParamP, g4, ls_parname= ['Len'], ls_sdpar= [0.5])

def modif_ParamP_sdMulti(ParamP, g4, ls_parname, ls_sdpar, corrmatrix=None):
    """ modif paramP loi normale multivariee """
    name1 = g4['name'] #nom de la bonne pop/sp
    nbp = len(ParamP) #nb de plantes max a resimule
    nbpar = len(ls_parname)

    # calul d'une matrice produits sigmax,sigmay
    res = []
    for i in ls_sdpar:
        res.append(i * np.array(ls_sdpar))

    sigmaproduct = np.array(res)

    #matrice de correlation
    if corrmatrix is None:
        correl_matrix = np.ones((nbpar, nbpar)) * 0. #covariances seront a zero par defaut
        np.fill_diagonal(correl_matrix, 1.)
    else:#matrice fournie en entree
        correl_matrix = np.array(corrmatrix) #sinon, matrice de correlation a fournir en entree
        if correl_matrix.shape[0] != nbpar or correl_matrix.shape[1] != nbpar:
            print('matrice de correlation pas a la bonne dimension!')

    #definition de la matrice de covariance
    matcov = sigmaproduct * correl_matrix

    # construction du vecteur des valeurs moyennes
    idpOK = 0 #id premiere plante bonne pop
    for nump in range(len(ParamP)):
        if ParamP[nump]['name'] == name1:  # issu du bon onglet
            idpOK = nump
            break

    mean_ = []
    for param in ls_parname:
        mean_.append(ParamP[idpOK][param])

    #tirage multivarie
    x = np.random.multivariate_normal(np.array(mean_), matcov, nbp)

    df = pd.DataFrame(x, columns=ls_parname)

    #boucle pour reaffecter les valeurs tirees dans g, puis ParamP
    for nump in range(len(ParamP)):
        if ParamP[nump]['name'] == name1:  # issu du bon onglet
            g = deepcopy(g4)
            for j in range(len(ls_parname)):
                parname = ls_parname[j]
                if parname in ['Largfeuille']:  # liste param qui peuvent etre negatifs -> pas de contrainte de positivite
                    g[parname] = df[parname][nump]
                else:
                    g[parname] = max(0.0000000001, df[parname][nump]) # seulement pour paramtere scalaie (pas liste) et positif (valeur >0)

            ParamP[nump] = g

    return ParamP, df


def dic2vec(nbplantes, dic):
    """ mise en liste 'nump'par plante un dico deja par plante """
    res3 = []
    for nump in range(nbplantes):
        try:
            key_ = str(nump)
            res3.append(dic[key_])
        except:
            res3.append(0.)

    return res3


# plus utilise ds l-egume
def dic_sum(ls_dict):
    "somme par cle les element de dico d'un meme format ; e.g [{0: 0, 1: 1, 2: 2}, {0: 3, 1: 4, 2: 5}] "
    # prepa d'un dico nul avec meme cles
    res = {}
    for k in list(ls_dict[0].keys()): res[k] = 0.
    # somme des dico
    for k in list(ls_dict[0].keys()):
        for i in range(len(ls_dict)):
            res[k] += ls_dict[i][k]

    return res


def append_dic(dic, key, element):
    """ add an element to a list in a dictionnary or create the list if the key is not present"""
    try:
        dic[key].append(element)
    except:
        dic[key] = [element]


def add_dic(dadd, dini):
    """ add the values of the k keys of a dictionary dadd to an existing dictionnary dini with the same keys - creates keys in dini if not already existing"""
    for k in list(dadd.keys()):
        try:
            dini[k] += dadd[k]
        except:
            dini[k] = dadd[k]
    return dini


def sum_ls_dic(dic):
    """ sum of the element by keys in a dictionnary of lists """
    for k in list(dic.keys()):
        dic[k] = sum(dic[k])



#### to read STICS files


def fill_metAN_fromSTICS(input_folder, tabSTICS):
    """ Converts STICS meteo file table into VGL meteo dictionnary

    :param input_folder: input folder containing the 'meteo_exemple.xls' file
    :param tabSTICS: panda data.frame of the STICS meteo file of a given year
    :return: VGL meteo dictionnary
    """

    # ! necessite 'meteo_exemple.xls' dans les inputs!!
    path_m = os.path.join(input_folder, 'meteo_exemple.xls')  #

    tab = tabSTICS
    lentab = tab.shape[0]
    if lentab==365:
        basemet = read_met_file(path_m, 'default365') #deepcopy(met365)
    elif lentab==366:
        basemet = read_met_file(path_m, 'default366') #deepcopy(met366)
    else:
        print("annee meteo incomplete")

    #an
    basemet['year'] = tab['AN'].tolist()
    #Tmin, Tmax, Tmoy
    basemet['Tmin'] = tab['TN'].tolist() #degresC
    basemet['Tmax'] = tab['TX'].tolist() #degresC
    Tmoy = (tab['TN'] + tab['TX'])/2.
    basemet['TmoyDay'] = Tmoy.tolist() #degresC
    # Tsol
    basemet['Tsol'] = Tmoy.tolist() #degresC #-> sans mesure, supose egal a air...
    #Et0, PP
    basemet['Et0'] = tab['ETPP'].tolist() #mm
    basemet['Precip'] = tab['RR'].tolist() #mm
    #RG, I0
    RG = tab['RG']*100 # conversion J.cm-2 #!! unite a ce niveau
    basemet['RG'] = RG.tolist() #J.cm-2
    I0 = 0.48*RG/(3600*24)*10000 # conversion W PAR.m-2
    basemet['I0'] = I0.tolist() #W PAR.m-2
    #V et CO2 (not used)
    basemet['V'] = tab['V'].tolist() #m.s-1
    basemet['CO2'] = tab['CO2'].tolist() #ppm

    #to to: prevoir filling quand ETPP manquant! + test unite RG

    return basemet


def fill_mng_fromSTICS(input_folder, mngSTICS, met):
    """ Converts STICS tec file into VGL management dictionnary

    :param input_folder: input folder containing the mngSTICS file
    :param mngSTICS: name of the STICS tec file
    :param met: VGL meteo dictionnary
    :return: VGL management dictionnary, list of information extracted from XML file
    """

    #marche pour usm avec dates fixes d'intervention (pas management automatique)

    # construction du fichier management base a partir de met et xml stics
    # default base
    basemng = {'Date': met['Date'], 'year': met['year'], 'month': met['month'], 'day': met['day'], 'DOY': met['DOY']}
    basemng['Coupe'] = [0] * len(basemng['DOY']) #yes/no : 1/0
    basemng['Irrig'] = [0] * len(basemng['DOY']) #mm
    basemng['FertNO3'] = [0] * len(basemng['DOY']) #kg N/ha
    basemng['FertNH4'] = [0] * len(basemng['DOY']) #kg N/ha
    basemng['Hcut'] = [3.] * len(basemng['DOY']) #cm
    basemng['WidthCut'] = [1000.] * len(basemng['DOY']) #cm # not used

    #mise a jour avec coupe et ferti/irrig du managment
    path_mng = os.path.join(input_folder, mngSTICS)
    treem = ET.parse(path_mng)
    rootm = treem.getroot()

    ## dates de coupe en dur depuis xml
    mngcut = rootm.find("./formalisme[@nom='special techniques']")
    mngcutta = mngcut.find("./option/choix/option/choix/ta[@nom='cutting management']")

    #nb_cut
    nb_cut = int(mngcutta.attrib['nb_interventions'])
    # recup des DOY et Hcut
    ls_DOYcut = []
    ls_Hcut = []
    if nb_cut>0:
        cuts_ = mngcutta.findall("./intervention/colonne[@nom='julfauche']")
        for cut_ in cuts_:
            ls_DOYcut.append(int(cut_.text))


        cuts_ = mngcutta.findall("./intervention/colonne[@nom='hautcoupe']")
        for cut_ in cuts_:
            hcut = float(cut_.text)*100. -2. #cm
            ls_Hcut.append(hcut)

        #remplacement des DOYcut et Hcut dans basemng
        for i in range(len(ls_DOYcut)):
            idcut = ls_DOYcut[i]-1
            if idcut < basemng['DOY'][-1]: # si bien dans la periode de simul
                basemng['Coupe'][idcut] = 1
                basemng['Hcut'][idcut] = ls_Hcut[i]

    # to do: autres options de coupe: auto / TT ... a faire a partir de management auto VGL


    ## irrigation
    mngirr = rootm.find("./formalisme[@nom='irrigation']")
    mngirrta = mngirr.find("./option/choix/ta[@nom='water inputs']")
    #nb_irrig
    nb_irrig = int(mngirrta.attrib['nb_interventions'])
    ls_DOYirr = []
    ls_amountirr = []
    if nb_irrig>0:

        irrs_ = mngirrta.findall("./intervention/colonne[@nom='julapI_or_sum_upvt']")
        for irr_ in irrs_:
            ls_DOYirr.append(int(irr_.text))

        irrs_ = mngirrta.findall("./intervention/colonne[@nom='amount']")
        for irr_ in irrs_:
            ls_amountirr.append(float(irr_.text))

        # remplacement des DOYirr et amount dans basemng
        for i in range(len(ls_DOYirr)):
            idcut = ls_DOYirr[i]-1
            if idcut < basemng['DOY'][-1]: # si bien dans la periode de simul
                basemng['Irrig'][idcut] = ls_amountirr[i]



    ## fertilisation
    mngfert = rootm.find("./formalisme[@nom='fertilisation']")
    mngfertta = mngfert.find("./ta[@nom='mineral nitrogen inputs']")
    #nb_fert
    nb_fert = int(mngfertta.attrib['nb_interventions'])
    ls_DOYfert = []
    ls_Amountfert = []
    ls_typefert = []
    if nb_fert>0:

        ferts_ = mngfertta.findall("./intervention/colonne[@nom='julapN_or_sum_upvt']")
        for fert_ in ferts_:
            ls_DOYfert.append(int(fert_.text))

        ferts_ = mngfertta.findall("./intervention/colonne[@nom='absolute_value/%']")
        for fert_ in ferts_:
            ls_Amountfert.append(float(fert_.text))

        ferts_ = mngfertta.findall("./intervention/colonne[@nom='engrais']")
        for fert_ in ferts_:
            ls_typefert.append(int(fert_.text))

        # remplacement des DOYfert et amount dans basemng
        for i in range(len(ls_DOYfert)):
            idcut = ls_DOYfert[i] - 1
            if idcut < basemng['DOY'][-1]:  # si bien dans la periode de simul
                # id '1': ammonium nitrate
                # id '2': UAN solution
                # id '3': urea
                # id '4': anydride ammonium
                # id '5': ammonium sulfate
                # id '6': calcium nitrate
                # id '8': fixed efficiency fertilizer
                # suppose doses toutes en kg N.ha-1 -> a verifier (pas poids masse absolue?)
                if ls_typefert[i] in [4,5]: #100% ammonium
                    basemng['FertNH4'][idcut] = ls_amountirr[i]
                elif ls_typefert[i] in [6]: #100% nitrate
                    basemng['FertNO3'][idcut] = ls_amountirr[i]
                elif ls_typefert[i] in [1]: #50% nitrate / 50% ammonium
                    basemng['FertNO3'][idcut] = 0.5*ls_amountirr[i]
                    basemng['FertNH4'][idcut] = 0.5*ls_amountirr[i]
                else:
                    basemng['FertNH4'][idcut] = ls_amountirr[i] #! autres fertilisant pour le moment mis en ammonium

                # check unite apport bien kg N.ha-1!

    info_mng = [[nb_cut, ls_DOYcut, ls_Hcut], [nb_irrig, ls_DOYirr, ls_amountirr], [nb_fert, ls_DOYfert, ls_Amountfert, ls_typefert]]

    return basemng, info_mng
    #mng, info_mng = fill_mng_fromSTICS(path_leg, nom_mng, met)


def read_usmXML_fromSTICS(input_folder, usm_file, usm_name):
    """ Converts STICS usm file into input dictionnaries for VGL model

    :param input_folder: input folder containing the STICS and 'meteo_exemple.xls' files
    :param usm_file: name of the STICS usm file (.xml)
    :param usm_name: name of STICS usm
    :return: VGL meteo dictionnary, VGL management dictionnary, list of information extracted from XML files
    """

    # recumperer nom_m1 et nom_m2 dans fichier XML
    # usm_file = 'usms.xml'
    # usm_name = 'DivLegLuz15'

    path_usm = os.path.join(input_folder, usm_file)
    tree = ET.parse(path_usm)
    root = tree.getroot()

    # lecture noms fichiers dans usm xml
    usm = root.find("./usm[@nom='" + usm_name + "']")
    nom_m1 = usm.find("./fclim1").text
    nom_m2 = usm.find("./fclim2").text
    nom_sol = usm.find("./nomsol").text
    # prend mng plante1 (stics cuture pure)
    usmplt1 = usm.find("./plante[@dominance='1']")
    nom_mng = usmplt1.find("./ftec").text
    usmdoydeb = int(usm.find("./datedebut").text)
    usmdoyend = int(usm.find("./datefin").text)

    # lecture fichiers meteo STICS
    # faire par anne puis fonction concat des dico
    tab1 = pd.read_csv(os.path.join(input_folder, nom_m1), sep=' ', decimal='.', header=None)
    tab1 = tab1.iloc[:, 0:13]  # garde les 13 colonnes
    tab1.columns = ["NUM_POSTE", "AN", "MOIS", "JOUR", "DOY", "TN", "TX", "RG", "ETPP", "RR", "V", "VPD", "CO2"]
    lentab1 = tab1.shape[0]

    tab2 = pd.read_csv(os.path.join(input_folder, nom_m2), sep=' ', decimal='.', header=None)
    tab2 = tab2.iloc[:, 0:13]  # garde les 13 colonnes
    tab2.columns = ["NUM_POSTE", "AN", "MOIS", "JOUR", "DOY", "TN", "TX", "RG", "ETPP", "RR", "V", "VPD", "CO2"]
    lentab2 = tab2.shape[0]

    # conversion fichiers meteo pour VGL
    met1 = fill_metAN_fromSTICS(input_folder, tab1)
    met2 = fill_metAN_fromSTICS(input_folder, tab2)

    # si 2 ans -> concatene les 2 annees
    if nom_m1 != nom_m2:
        met = add_dic(met2, met1)
        met['DOY'] = list(range(1, len(met['DOY']) + 1))  # update DOY
    else:
        met = met1

    # definition management
    mng, info = fill_mng_fromSTICS(input_folder, nom_mng, met)
    info.append([nom_sol, nom_m1, nom_m2, nom_mng, usmdoydeb, usmdoyend])

    #renvoie pas le sol -> traite dans sol3ds
    return met, mng, info
    #met, mng, info = read_usmXML_fromSTICS(path_leg, 'usms.xml', 'DivLegLuz15')
    #ajouter ecriture des fichiers en option
