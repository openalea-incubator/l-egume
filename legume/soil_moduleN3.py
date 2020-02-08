from scipy import *
#from rpy import r
import sys
path_ = r'C:\devel\l-egume\legume'
path2_ = r'C:\devel\l-egume\legume'
sys.path.insert(0, path_)
sys.path.insert(0, path2_)
from soil_module5 import *


##changement version3
#depend de soil_module5
#ajout entree  obstarac + effet obstrac sur dispo Nmin pour uptke plante->OK
#debug L708 dans #stepNuptakePlt -> test des None avec is et pas ==
#retirer dependence a R -> a priori pas utilise
#debug FN_factor: ajout d'un epsilon pour eviter division par zero
#ajout d'une methode mixResMat() pour ger ajour d'un redisu deja existant
#ajout d'un bilanN['NminfromNres'] pour comptabiliser Nmin produit par chaque type de residu (liste de liste), temporel et son equivalent cumule bilanN['NminfromNresCum']
#modif des bilanN['cumNRes123'] pour les avoir en journalier (et pas en vrac de tous les residus
#debug bilanN['NminfromNresCum'] dans closeNbalance qd pas de residus (try:)

#Afaire:


#revoir test uni




## changements version 2bis
#depend de soil_module3 ->OK, en fait soil_module4 (pour couplage avec plante)
#pour init, besoin de : par_sol, soil_number = [13]*10, dxyz = [[1.], [1.], [0.2]*10], vDA=[1.2]*10, ZESX=0.3, CFES=1.
# -> rajouter vDA, ZESX, CFES ; a initialisation du sol
# -> retire DA en dur + passer pHeau et CALCs en parametres pour initialisation de sol N
# -> passer vecteur de profil vMO[z], vCN[z], parSN, vNH4[z], vNO3[z] = utilise a l'initialisation, mais pas passer en parametre (variables gloable
# modifier partout ou DA, HR() qui etaient utilise -> OK, sequence test fonctionne
# cree methode mask_PROFUM(par_SN) pour appliquer humification sur profhum uniquement
# filtre de PROFHUM-> OK
# applique sur Corg/Norg, plutot que sur KH2HUM
# mettre intitalisation des m_NO3 / m_NH4 aux bonnes unite!! (init en kg d'N.ha-1)-> a passer en kg d'N sur la bonne surface ->OK
# cree methode updateTsol pour mettre a jour m_Tsol
# ajoute InertCorg, InertNorg pour les retirer des pools actifs de MO
# debug stepNitrif -> retirer de NH4 convertie etait oubli! -> OK
# debug flux d'N ds le sol: flux lixiviation etaient pas bon -> pb dans distrib_NO3 -> correction OK, bilan equilibre
# ajout de fonctions pour ouvrir / clore / afficher bilan C et N -> base sur 2 dictionnaires qui stockent la dynamqiue des variables
# bilan C N avec et sans residus -> OK
# ajout des methode ls_NRES() et ls_NBio() qui etaient necessaire au bilan Norg
#ajout d'entree differencie d'azote mineral pour rain, riirigation, fertlizers dans stepNINFILT


#ajout de prelevement d'N plante :stepNuptakePlt() et toutes les fonctions assoiees

#Rq! !! PAs de garde fou pour empecher m_N4 et m_NO3 de passer sous zero!
# peut arriver avec residus (stepMicrobioMin: besoin de Neccessi) -> verif

#23/02/18
#ajoute la sortie de idmin dans stepNuptakePlt (pour faciliter deboggage)
#change Distrib_Potential_Nuptake_Plt pour fair la partition en fonction des activites d'absorption d'N mineral, pas uniquement des longueurs de racines


#decrire les differentes matrices!!



class SoilN(Soil):
    def __init__(self, par_sol, parSN, soil_number , dxyz, vDA,  vCN, vMO, vARGIs, vNO3, vNH4,  vCALCs, Tsol, pH, ZESX, CFES, obstarac=None, pattern8=[[0,0],[100.,100.]]):
        """
        par_sol SN contient en plus pour l'azote:
            'FMIN1G'        #(day-1) (p145)
            'FMIN2G'        #(%ARGIS-1) para pot rate min a ARGIs (p145)
            'FMIN3G'        #(% CALC-1)para pot rate min a CALCss (p145)
            'FINERTG'       #0.65 = default fraction of N pool inactive for minearlisation (p145) (could be smaller in grassland & forest)
            'PROFHUMs'      # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220
            'HMinMg'        #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
            'HoptMg'        #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)
            'TRefg'         #reference temperature (degreC)
            'FTEMHAg'       #asymptotic value of FTH (seen as a logistic response)
            'FTEMHg'        #(K-1) 
            'FTEMHB'
            'FNXg'          #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
            'PHMinNITg'     #pH min de nitrification (prop of Field Capacity) #value p149
            'PHMaxNITg'     #pH max de nitrification (prop of Field Capacity) #value p149
            'HMinNg'        #Humidite min de nitrification #value p149
            'HoptNg'        #Humidite opt de nitrification  #value p149
            'TNITMINg'      #Temperature min de nitrification #  (degreC)value p151
            'TNITOPTg'      #Temperature opt de nitrification #  (degreC)value p151
            'TNITMAXg'      #Temperature max de nitrification #  (degreC)value p151
            'RATIONITs'     #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
            'DIFNg'         #N diffusion coefficient at field capacity (cm2.day-1, p 161)
        
        CALCs
        pHeau
        ARGIs
        m_Tsol

        Corg: (kg C dans le voxel)
        Norg: (kg N dans le voxel)
        InertCorg: (kg C dans le voxel) 
        InertNorg: (kg N dans le voxel)
        m_NH4: (kg N NH4 dans le voxel)
        m_NO3: (kg N NO3 dans le voxel)
        m_MO
        m_CNHUM

        K2HUM
        N2ONitrif
        N2ODenitrif
        lixiNO3

        bilanC et bilanN: dictionnaires contenant les variables dynamiques et les cumuls necessaire a l'etablissement des bilans C et N (kg C/N.ha-1)
        """
        #initialisation sol et teneur en eau
        Soil.__init__(self,par_sol, soil_number , dxyz, vDA, ZESX, CFES, obstarac, pattern8)
        self.compute_teta_lim(par_sol)
        self.init_asw()

        self.CALCs = vCALCs[0] #(%) value for non calcareous soils (p220 STICS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A mesurer/mettre a jour (rq: mesure CaO anterieure parcelle C2: 1.58 g.kg-1 dans 0-30)
        self.pHeau = pH
        self.ARG = vARGIs[0]
        #a faire: NORGs matrices de NORGs, CORGs, K2HUM...
        
        #initialisation matrice Temp et et compartiments azote sol
        self.m_MO, self.m_CNHUM, self.K2HUM, self.m_Tsol, self.m_NH4, self.m_NO3 = [], [], [], [], [], []
        for z in range(len(dxyz[2])):
            v, v1, v2, v3, v4, v5 = [], [], [], [], [], []
            for x in range(len(dxyz[0])):
                vv, vv1, vv2, vv3, vv4, vv5 = [], [], [], [], [], []
                for y in range(len(dxyz[1])):
                    surf= self.dxyz[0][x] * self.dxyz[1][y]
                    vv.append(vMO[z])
                    vv1.append(vCN[z])
                    vv2.append(self.Pot_rate_SOMMin(vCALCs[z], vARGIs[z], parSN))
                    vv3.append(Tsol) #suppose cste et egale a Tair -> a adapter
                    vv4.append(vNH4[z]/10000.*surf) #en kg d'N par voxel
                    vv5.append(vNO3[z]/10000.*surf) #en kg d'N par voxel

                v.append(vv)
                v1.append(vv1)
                v2.append(vv2)
                v3.append(vv3)
                v4.append(vv4)
                v5.append(vv5)

            self.m_MO.append(v)
            self.m_CNHUM.append(v1)
            self.K2HUM.append(v2)
            self.m_Tsol.append(v3)
            self.m_NH4.append(v4)
            self.m_NO3.append(v5)

        self.m_MO, self.m_CNHUM,  self.m_Tsol, self.m_NH4, self.m_NO3  = array(self.m_MO), array(self.m_CNHUM),  array(self.m_Tsol), array(self.m_NH4), array(self.m_NO3)
        self.K2HUM = array(self.K2HUM) 
        self.Corg = multiply(self.m_soil_vol , self.m_MO)* self.mask_PROFUM(parSN) * self.m_DA / 1.72 #en kg de C (#aplique sur Profhum)
        self.Norg = divide(self.Corg, self.m_CNHUM) #en kg de N
        self.InertCorg = self.Corg * parSN['FINERTG']
        self.InertNorg = self.Norg * parSN['FINERTG']
        self.OpenCbalance()
        self.OpenNbalance()
        
        self.N2ONitrif, self.N2ODenitrif = 0., 0. #kg N for the whole soil volume
        self.lixiNO3 = 0. #kg N
        #reprendre les Temperature avec un modele plus elabore!! ici = Tair!


    def Pot_rate_SOMMin(self, CALCs, ARGIs, par):
        """ Potential rate of SOM mineralisation - eq. 8.5 p145 """
        K2HUMi = par['FMIN1G']*exp(-par['FMIN2G']*ARGIs)/(1+par['FMIN3G']*CALCs)
        return K2HUMi
        #!! revoir ARGIs et CALCs!!

    def SOMMin_RespT(self, par):
        """ reponse de la mineralisation (ammonification) de la SOM a la temperature - described as a sigmoid process (FTH) - Eq 8.3 p143 et (pour residus FTR p146-147) """
        if self.m_Tsol.min()<=0.:#min(min(min(self.m_Tsol)))<=0.:
            FTH=0.*self.m_1
        else:
            FTH = par['FTEMHAg'] / (1 + par['FTEMHB']*exp(-par['FTEMHg']*self.m_Tsol))
        return FTH 

    def SOMMin_RespHum(self, par):
        """ reponse de la mineralisation (ammonification) de la SOM a l'humidite relative (FH) - Eq. 8.2 p 143 - aussi utilise pour residus """
        HR = self.HRv()
        FH = (HR - par['HMinMg']*100) / (par['HoptMg']*100 - par['HMinMg']*100)
        for i in range(len(FH)):
            for j in range(len(FH[i])):
                for k in range(len(FH[i][j])):
                    if FH[i][j][k]<0.:
                        FH[i][j][k]=0.
                    elif FH[i][j][k]>1.:
                        FH[i][j][k]=1.
        return FH

    def Act_rate_SOMMin (self, par):
        """ (%??) """
        return self.K2HUM*self.SOMMin_RespHum(par)*self.SOMMin_RespT(par)

    def stepNB(self, par):
        """ """
        #Mineralisation of soil organic matter
        NHUM = self.Norg - self.InertNorg#* (1-par['FINERTG']) #!!! *PROFHUMs : mask pour retirer profondeur ou mineralisation pas significative???
        dN_NH4 = NHUM * self.Act_rate_SOMMin(par)
        dC_NH4 = dN_NH4*self.m_CNHUM
        self.Norg = self.Norg - dN_NH4
        self.Corg = self.Corg - dC_NH4 #suppose CN constant/pas affecte par SOM mineralisation!
        self.m_MO = self.Corg*1.72 # pas vraiment utile
        self.m_NH4 = self.m_NH4 + dN_NH4 #verif unite!!! ici en kg d'N dans le voxel de sol ->OK
        #Update bilanCN
        self.bilanN['cumMinN'].append(sum3(dN_NH4) / self.soilSurface() *10000)
        self.bilanC['cumMinC'].append(sum3(dC_NH4) / self.soilSurface() *10000)

    def init_residues(self, vCNRESt=[], vAmount=[], vProps=[], vWC=[], vCC=[], forced_Cres=None):
        """ initialisation des compartiments en relation avec gestion des residus """
        #dictionnaire de parametres des residus
        self.parResi = {}
        self.parResi['CNRESt'], self.parResi['CNBio'], self.parResi['KRES'], self.parResi['YRES'], self.parResi['KBio'], self.parResi['HRES'] = [],[],[],[],[],[]
        self.parResi['TRefg'] = 15. #reference temperature 
        self.parResi['FTEMHAg'] = 12. #asymptotic value of FTR (seen as a logistic response) p147
        self.parResi['FTEMHg'] = 0.103 #(K-1) p147
        self.parResi['FTEMHB'] = 52. #from eq 8.3

        for val in vCNRESt:
            self.addResPAR(self.parResi, val)

        #liste de matrices pour CRES et CBio
        self.ls_CRES, self.ls_CBio = [], []
        self.bilanN['NminfromNres'] = []  # liste de liste de deltaNmin produits par residus pour les bilans
        for i in range(len(vAmount)):
            self.addResMat(vAmount[i], vProps[i], vWC[i], vCC[i], forced_Cres) #utilise fonction de distribution verticale par default, mais peut forcer une matrice donnee

        # bilan!! -> a revoir (pour le moment gere seulement qd tous les residus sont aplliques en meme temps-> OK : MxResMat() pemet d'initiliser en meme temps (meme a zeo) puis d'ajouter qd on veut un amount
        self.bilanN['initialNres'] = sum(self.ls_NRES()) / self.soilSurface() * 10000

        #CO2 resp
        self.CO2respSoil = 0. #separer (kg de C par le volume total de sol) #rq: suppose aucune recapture

    def addResPAR(self, par, CNRes):
        """ add a new series of parammeters for a residue according to Nicolardot et al. (2001) in function of its CN ratio"""
        par['CNRESt'].append(CNRes) #CSURNRESt C/N des residus #liste pour les different residus 
        par['CNBio'].append(max(7.8, 16.1-123./CNRes)) #CN ratio of the zymogeneous biomass -> fontion du CN des residus (eq. 8.6) -> des biomasses microbiennes associee a chaque type de residu?
        par['KRES'].append(0.07+1.94/CNRes) # decomposition rate constant (day-1- normalised day at 15dC) from organic residue to microbil biomass?  (fig 8.4 p146) -> fonction de CNRESt
        par['YRES'].append(0.62) # Assimilation yield of residue-C by microbial biomass - partition parameter between CO2 and microbil biomass (fig 8.4 p146) -> constant
        par['KBio'].append(0.0110) # decomposition rate constant (day-1 - normalised day at 15dC) from microbial biomass to humus? (fig 8.4 p146) -> constant
        par['HRES'].append(1-(0.69*CNRes)/(11.2+CNRes)) #  partition parameter between CO2 and humus - Humification rate of microbial biomass -(fig 8.4 p146) -> fonction de CNRESt
        return par

    def addResMat(self, Amount, Vprop, Wcontent=0.8,Ccontent=0.42, forced_Cres=None):
        """ add matrice associee au C des residus (CRES) et de leur microbial biomass """
        self.ls_CBio.append(0*self.m_1)#matrice for the Amount of C in the microbial biomass (kg C in the voxel) -> initialise a zero
        if forced_Cres == None:
            cres = VdistribResidues(self, Amount, Vprop, Wcontent, Ccontent)
            self.ls_CRES.append(cres)
            #bilan
            self.bilanC['initialCres'] += sum(cres)/ self.soilSurface() *10000
            self.bilanN['NminfromNres'].append([])  # ajoute une liste vide pour le residu
        else: #pour gerer le cas ou passera une matrice (3D) deja faite a partir de maquettes
            self.ls_CRES.append(forced_Cres)
            #bilan
            self.bilanC['initialCres'] += sum(forced_Cres)/ self.soilSurface() *10000
            self.bilanN['NminfromNres'].append([])  # ajoute une liste vide pour le residu

    def Pot_rate_ResidueMin(self, res_id, par):
        """ changement potentiel en C du residu id (dans conditions donnes d'humidite et T) eq. 8.7 p 146"""
        #par : par_Sol, pour reponse a humidite
        pDCRES = -self.parResi['KRES'][res_id] * self.ls_CRES[res_id] * self.SOMMin_RespHum(par) * self.SOMMin_RespT(self.parResi) #sans *FN (dispo en N) ->  pot microbial growth dans ces condition
        return pDCRES
        #faudrait lire les FH plutot que de les recalculer a chaque fois
        #laisser en valeur positive?

    def Pot_Ndemand_microbialBio(self, par):
        """  demande totale en azote pour atteindre croissance optimale microbio de tous les residus (kg de N par voxel) - somme des demandes pour chaque residu """
        #par : par_Sol, pour reponse a humidite
        res = self.m_1*0
        for i in range(len(self.parResi['KRES'])):
            pDCRESi = -self.Pot_rate_ResidueMin(i, par)
            CBio_poti = pDCRESi *   self.parResi['YRES'][i]
            deltaNi = CBio_poti/self.parResi['CNBio'][i] - pDCRESi/self.parResi['CNRESt'][i] #ce qu'il manque entre N issu de biomasse decomposee et N requis pour microbio! (>0 de ce qu'il faut en N pour atteindre potentiel)
            res = res + deltaNi

        return res

    def FN_factor(self, par):
        """ calcul du facteur de reduction lie a la disponibilite en Nmin a proximite des residus """
        #!! demande des plantes pas prise en compte -> servie slmt s'il en reste apres microbio!
        MND = self.Pot_Ndemand_microbialBio(par)
        Nmin = self.m_NH4 + self.m_NO3
        delta = Nmin - MND #<0 -> manque de N 
        ratio = Nmin/(MND+1e-15)
        FN = deepcopy(delta)
        for i in range(len(FN)):
            for j in range(len(FN[i])):
                for k in range(len(FN[i][j])):
                    if FN[i][j][k]>=0.:#pas de limitation
                        FN[i][j][k] = 1.
                    else:#<0 -> manque de N -> au max prend Nmin d'ou reduction de Nmin/MND
                        FN[i][j][k] = ratio[i][j][k] #rajouter un petit% de securite?
        return FN

    def FBIO_factor(self, par):
        """ p147 - a faire et introduire dans stepMicrobioMin"""
        pass

    def stepResidueMin(self, par):
        """  """
        #par : par_Sol, pour reponse a humidite
        FN = self.FN_factor(par) #!! faudrait calculer FN une seule fois pour tous les residus!
        cumNRes1, cumNRes2= [],[]
        for i in range(len(self.ls_CRES)):
            #fux de C
            DCRESi = self.Pot_rate_ResidueMin(i, par) * FN
            self.ls_CRES[i] = self.ls_CRES[i] + DCRESi #ou - si garde DCRESi positif
            DCBioi = DCRESi*self.parResi['YRES'][i]
            self.ls_CBio[i] = self.ls_CBio[i] - DCBioi #ou + si garde DCRESi positif
            DCO2 = DCRESi-DCBioi
            self.CO2respSoil = self.CO2respSoil - sum3(DCO2)#ou + si garde DCRESi positif
            #bilan
            self.bilanC['cumCO2Res1'].append(- sum3(DCO2)/ self.soilSurface() *10000)

            #flux de N
            #ajouter a NH4 le complement N de CO2 emission?
            DNCO2 = DCO2/self.parResi['CNRESt'][i]
            self.m_NH4 = self.m_NH4 - DNCO2 #ou + si garde DCRESi positif
            Nmin = self.m_NH4 + self.m_NO3 + 0.000000000001 #pour eviter les Nan dans f_NH4 et f_NO3
            Ndemandi = -DCBioi/self.parResi['CNBio'][i] + DCBioi/self.parResi['CNRESt'][i] #verif signes!
            #preleve dans NH4 et NO3 a hauteur de leur contribution a Nmin + verif que aucun devient negatif
            f_NH4 = self.m_NH4/Nmin
            f_NO3 = self.m_NO3/Nmin
            self.m_NH4 = self.m_NH4 - f_NH4*Ndemandi
            self.m_NO3 = self.m_NO3 - f_NO3*Ndemandi
            #bilan
            cumNRes1.append(- sum3(DNCO2) / self.soilSurface() *10000)
            cumNRes2.append(- sum3(Ndemandi)/ self.soilSurface() *10000)
            self.bilanN['NminfromNres'][i].append(- sum3(DNCO2) / self.soilSurface() *10000)

        self.bilanN['cumNRes1'].append(sum(cumNRes1))#pour ajouter uniquement cumul journalier de tous les residus
        self.bilanN['cumNRes2'].append(sum(cumNRes2))


    def stepMicrobioMin(self, par):
        """ """
        #par : par_Sol, pour reponse a humidite
        cumNRes3= []
        for i in range(len(self.ls_CBio)):
            #fux de C
            DCBioi = self.parResi['KBio'][i] * self.ls_CBio[i] * self.SOMMin_RespHum(par) * self.SOMMin_RespT(self.parResi)
            self.ls_CBio[i] = self.ls_CBio[i] - DCBioi
            DCHUM = DCBioi * self.parResi['HRES'][i]
            self.Corg = self.Corg + DCHUM
            DCO2 = DCBioi - DCHUM
            self.CO2respSoil = self.CO2respSoil + sum3(DCO2)
            #bilan
            self.bilanC['cumCO2Res2'].append(sum3(DCO2)/ self.soilSurface() *10000)

            #flux de N
            #equilibre Norg 
            self.Norg = self.Norg + DCHUM/self.m_CNHUM
            #ajouter a NH4 le complement N de CO2 emission? + ajouter N lie au changement de CN de DCHUM
            Neccesi = DCHUM/self.parResi['CNBio'][i] - DCHUM/self.m_CNHUM #en ecces par rapport a CN de MOS #!!! si <0 pourrait potentiellement faire passer Nmin <0 -> a ajuster avec facteur FBIO!
            self.m_NH4 = self.m_NH4 + DCO2/self.parResi['CNBio'][i] + Neccesi 
            
            #bilan
            cumNRes3.append(sum3(DCO2/self.parResi['CNBio'][i] + Neccesi)/ self.soilSurface() *10000)

        self.bilanN['cumNRes3'].append(sum(cumNRes3))

        #facteur FBIO!! = 1 par defaut -> peut augmenter <-> baisse de CN ratio de la biomasse microbienne si Nmin exhauted (pour eviter bilan negatif)
        #a faire! avec approche similaire a FN_factor!

    def ls_NRES(self):
        """ calcul N dans les residus - liste equivalente de ls_CRES """
        lsNRES = []
        for i in range(len(self.ls_CRES)):
            lsNRES.append(self.ls_CRES[i]/self.parResi['CNRESt'][i])
        return lsNRES

    def ls_NBio(self):
        """ calcul N dans les biomasse microbienne de residus - liste equivalente de ls_CBio """
        lsNbio = []
        for i in range(len(self.ls_CBio)):
            lsNbio.append(self.ls_CBio[i]/self.parResi['CNBio'][i])
        return lsNbio

    def mixResMat(self, mat_res, idres, Ccontent=0.42):
        """ add matrice associee des residus (en g MS par voxl) au C (CRES) d' un residu avec id deja existant """
        cres = mat_res * Ccontent / 1000.  # conversion #kg of C per voxel
        self.ls_CRES[idres] += cres
        # bilan
        self.bilanC['initialCres'] += sum(cres) / self.soilSurface() * 10000
        self.bilanN['initialNres'] += sum(self.ls_CRES[idres] / self.parResi['CNRESt'][idres]) / self.soilSurface() * 10000

    # suppose CsurN comme residu existant; a faire (?) cree noueau residu si tres different? ajuster bilan C/N # faire evoluer lrd parametres des reisdus selon C/N vrai?

    def Nitrif_RespHum(self, par):
        """ reponse de la nitrification de NH4+ a l'humidite relative (FHN) - Eq. 8.14 p 150 - increasing sigmoid-like curve """
        HR = self.HRv()
        FHN = (HR - par['HMinNg']*100) / (par['HoptNg']*100 - par['HMinNg']*100)
        for i in range(len(FHN)):
            for j in range(len(FHN[i])):
                for k in range(len(FHN[i][j])):
                    if FHN[i][j][k]<0.:
                        FHN[i][j][k]=0.
                    elif FHN[i][j][k]>1.:
                        FHN[i][j][k]=1.
        return FHN

    def Nitrif_RespPH(self, par):
        """ reponse de la nitrification de NH4+ au pH (FH) - Eq. 8.13 p 150 - increasing sigmoid-like curve """
        pHs = self.pHeau
        FPHN = (pHs - par['PHMinNITg']) / (par['PHMaxNITg'] - par['PHMinNITg'])
        if FPHN<0.:
            FPHN=0.
        elif FPHN>1.:
            FPHN=1.
        return FPHN #a priori scalaire a ne calculer qu'une fois (comme pH change pas)

    def Nitrif_RespT(self, par):
        """ reponse de la nitrification de NH4+ a la temperature (FTN) - Eq. 8.15 p 151 - bilinear beta-like curve """
        FTN = deepcopy(self.m_Tsol)
        for i in range(len(FTN)):
            for j in range(len(FTN[i])):
                for k in range(len(FTN[i][j])):
                    Tsol = FTN[i][j][k]
                    if Tsol <= par['TNITOPTg']:
                        ratio = (Tsol - par['TNITMINg']) / (par['TNITOPTg'] - par['TNITMINg'])
                        FTN[i][j][k] = max(0., ratio)
                    else:#superieur a Topt
                        ratio = (Tsol - par['TNITMAXg']) / (par['TNITOPTg'] - par['TNITMAXg'])
                        FTN[i][j][k] = max(0., ratio)
        return FTN

    def stepNitrif(self, par):
        """ eq. 8.12 et 8.16 """
        TNITRIF = self.m_NH4 * par['FNXg'] * self.Nitrif_RespHum(par) * self.Nitrif_RespPH(par) * self.Nitrif_RespT(par)
        NITRIF = (1 - par['RATIONITs']) * TNITRIF
        self.m_NO3 = self.m_NO3 + NITRIF
        self.m_NH4 = self.m_NH4 - TNITRIF
        dN2ONitrif = sum3(TNITRIF-NITRIF)
        self.N2ONitrif = self.N2ONitrif + dN2ONitrif
        self.bilanN['cumN2O'].append(dN2ONitrif / self.soilSurface() *10000) #UpdateNminbalance


    def infil_layerNO3(self, in_N, out_Water , idz, opt=1):
        new = self.m_NO3[idz] + in_N
        propNO3 = out_Water / (self.m_QH20max[idz] + out_Water)#prop de nitrate qui part est fraction du volume d'eau max qui passe (jamais >1)
        #putmask(propNO3, propNO3>1. ,1.)#!!!verif pas superieur a 1 et sinon remplace par 1!!  syntaxe interessante
        out_N = new*propNO3
        new = new-out_N

        if opt==2: #ditribution pas juste verticale
            out_N2 = deepcopy(out_N)
            out_N2.fill(0.)
            for x in range(len(out_N2)):
                for y in range(len(out_N2[x])):
                    q_out = out_N[x][y]
                    ls_v = ls_1storder_vox(self.dxyz, x,y,idz, opt)#distribution entre les 1st order ; mettre opt=1 si veut forcer verticalement / 2 si
                    if len(ls_v)>1:
                        ponder = [0.0416666, 0.0416666, 0.0416666, 0.0416666, 2/3., 0.0416666, 0.0416666, 0.0416666, 0.0416666]# 2/3 en dessous 1/3 au premier ordre
                    else:
                        ponder = [1.]

                    for i in range(len(ls_v)):#distribution dans les voxels du dessous
                        nx, ny = ls_v[i][0], ls_v[i][1]
                        out_N2[nx][ny] = out_N2[nx][ny]+ponder[i]*q_out  

            out_N = out_N2

        return new, out_N


    def distrib_NO3(self, map_N, ls_outWater, opt=1):#map_N = map application nitrates en surface
        """  """
        
        in_N = map_N
        matNO3_t = deepcopy(self.m_NO3)
        #ls_out = []
        for z in range(len(matNO3_t)):
            new, out_ = self.infil_layerNO3(in_N, ls_outWater[z], z, opt)
            #ls_out.append(out_)
            in_N = out_
            matNO3_t[z] = new

        return matNO3_t, out_
        #A faire: distinguer les entree N pluie, irrigation, fertilisation


    def stepNINFILT(self, mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, ls_outWater, opt=1):#(self, map_N, ls_outWater, opt=1):
        """ """
        #ajout N NO3 mobile
        map_N = mapN_Rain + mapN_Irrig + mapN_fertNO3 #+ mapN_fertNH4
        matNO3_t, out_NO3 = self.distrib_NO3(map_N, ls_outWater, opt)
        self.m_NO3 = matNO3_t
        Lix = sum(sum(out_NO3))
        self.lixiNO3 = self.lixiNO3 + Lix

        #ajout N NH4 non mobile dans 1ere couche
        self.m_NH4[0,:,:] = self.m_NH4[0,:,:] + mapN_fertNH4

        #bilans
        self.bilanN['cumRain'].append(sum(mapN_Rain)/ self.soilSurface() *10000) 
        self.bilanN['cumIrrig'].append(sum(mapN_Irrig)/ self.soilSurface() *10000)
        self.bilanN['cumfertNO3'].append(sum(mapN_fertNO3)/ self.soilSurface() *10000) 
        self.bilanN['cumfertNH4'].append(sum(mapN_fertNH4)/ self.soilSurface() *10000)
        self.bilanN['cumLix'].append(Lix / self.soilSurface() *10000) 
     

    def mask_PROFUM(self, parSN):
        """ pour creer un mask pour profhum """
        PROFHUM = parSN['PROFHUMs']/100. #en m
        limz = [0.]
        for i in range(len(self.dxyz[2])): 
            limz.append(limz[-1]+self.dxyz[2][i])

        limz = array(limz)

        v = limz<PROFHUM 
        v=v*1.
        v = v.tolist()
        idlim = v.index(0)
        limz[idlim-1]
        v[idlim] = (PROFHUM - limz[idlim-1]) / (limz[idlim]-limz[idlim-1])
        v = v[1:] #mask 1D
        
        #applique a matrice sol
        res = self.m_1*1.
        for i in range(len(v)):
            res[i,:,:] = res[i,:,:]*v[i]

        return res

    def updateTsol(self, Tair):
        """ Tsol= Tair """
        self.m_Tsol = self.m_1*Tair
        # A ameliorer avec bilan d'E! et TCULT


    def OpenCbalance(self):
        """ Dictionnary for soil Carbon balance (kg C.ha-1)
        Keys for Daily outputs: 'cumMinC', 'cumCO2Res1', 'cumCO2Res2'
        Keys for Total Input: 'intialInertC', 'intialActiveC', 'initialCZygo', 'initialCres'
        Keys for Total outputs: 'FinalInertC', 'FinalActiveC', 'MinCtot'
        Keys for totals: 'InputCtot', 'OutputCtot'
        """
        surfsolref = self.soilSurface()
        self.bilanC = {}
        #Humus
        self.bilanC['intialInertC'] = sum3(self.InertCorg) / surfsolref *10000
        self.bilanC['intialActiveC'] = sum3(self.Corg-self.InertCorg) / surfsolref *10000
        self.bilanC['cumMinC'] = []
        #residus
        self.bilanC['cumCO2Res1'], self.bilanC['cumCO2Res2'] = [], []
        try:
            self.bilanC['initialCZygo'] = sum(self.ls_CBio)/ surfsolref *10000
            self.bilanC['initialCres'] = sum(self.ls_CRES)/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanC['initialCZygo'] = 0.
            self.bilanC['initialCres'] = 0.

        # autres termes avec microbio et residus a ajouter!

    #def UpdateCbalance(self, dCMin):
    #    surfsolref = self.soilSurface()
    #    self.bilanC['cumMinC'].append(sum3(dCMin) / self.soilSurface() *10000)
    #    # autres termes avec microbio et residus a ajouter!!

    def CloseCbalance(self, print_=1):
        surfsolref = self.soilSurface()
        #Humus Mineralisation
        #input  
        self.bilanC['InputCtot'] = self.bilanC['intialInertC'] + self.bilanC['intialActiveC'] + self.bilanC['initialCZygo'] + self.bilanC['initialCres']
        #output
        self.bilanC['FinalInertC'] = sum3(self.InertCorg) / surfsolref *10000
        self.bilanC['FinalActiveC'] = sum3(self.Corg-self.InertCorg) / surfsolref *10000

        #Residus
        try:
            self.bilanC['finalCZygo'] = sum(self.ls_CBio)/ surfsolref *10000
            self.bilanC['finalCres'] = sum(self.ls_CRES)/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanC['finalCZygo'] = 0.
            self.bilanC['finalCres'] = 0.

        self.bilanC['MinCtot'] = sum(self.bilanC['cumMinC']) + sum(self.bilanC['cumCO2Res1']) + sum(self.bilanC['cumCO2Res2'])
        self.bilanC['OutputCtot'] = self.bilanC['FinalInertC'] + self.bilanC['FinalActiveC'] + self.bilanC['MinCtot'] + self.bilanC['finalCZygo'] + self.bilanC['finalCres']

        if print_==1:
            self.PrintCbalance()
        #pourrait le diriger vers un fichier de sortie texte?

    def PrintCbalance(self):
        bilanC = self.bilanC
        print ("")
        print ("Carbon Balance Input (kg C.ha-1)\t\t\t Carbon Balance Output (kg C.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Active Humified Pool:\t {0:8.1f}\t\t\t Active Humified Pool:\t {1:8.1f}".format(bilanC['intialActiveC'], bilanC['FinalActiveC'])))
        print(("Inert Humified Pool:\t {0:8.1f}\t\t\t Inert Humified Pool:\t {1:8.1f}".format(bilanC['intialInertC'], bilanC['FinalInertC'])))
        print(("Zymogeneous Bio Pool:\t {0:8.1f}\t\t\t Zymogeneous Bio Pool:\t {1:8.1f}".format(bilanC['initialCZygo'], bilanC['finalCZygo'])))
        print(("Added organic matter:\t {0:8.1f}\t\t\t Added organic matter:\t {1:8.1f}".format(bilanC['initialCres'], bilanC['finalCres'])))
        print(("                            \t\t\t\t Mineralisation:\t\t {0:8.1f}".format(bilanC['MinCtot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t\t {1:8.1f}".format(bilanC['InputCtot'], bilanC['OutputCtot'])))
        print ("")

    def OpenNbalance(self):
        """ Dictionnary for soil Organic and Mineral balance (kg N.ha-1)
        Keys for Daily outputs: 
        Keys for Total Input: 
        Keys for Total outputs: 
        Keys for totals: 'InputNtot', 'OutputNtot', 'InputNmintot', 'OutputNmintot'
        """
        self.bilanN = {}
        self.bilanN['intialInertN'] = sum3(self.InertNorg) / self.soilSurface() *10000
        self.bilanN['intialActiveN'] = sum3(self.Norg-self.InertNorg) / self.soilSurface() *10000
        self.bilanN['cumMinN'] = []
        self.bilanN['intialNO3'] = sum3(self.m_NO3) / self.soilSurface() *10000
        self.bilanN['intialNH4'] = sum3(self.m_NH4) / self.soilSurface() *10000
        self.bilanN['cumLix'], self.bilanN['cumN2O'] = [],[]
        self.bilanN['cumRain'] = []
        self.bilanN['cumIrrig'] = []
        self.bilanN['cumfertNO3'] = []
        self.bilanN['cumfertNH4'] = []
        self.bilanN['cumUptakePlt'] = []
        self.bilanN['azomes'] = [] #equivalent a stics
        #residus
        self.bilanN['cumNRes1'], self.bilanN['cumNRes2'], self.bilanN['cumNRes3'] = [], [], []
        self.bilanN['NminfromNresCum'] = []

        try:
            self.bilanN['initialNZygo'] = sum(self.ls_NBio())/ self.soilSurface() *10000
            self.bilanN['initialNres'] = sum(self.ls_NRES())/ self.soilSurface() *10000
        except:#si pas de residus ajoute
            self.bilanN['initialNZygo'] = 0.
            self.bilanN['initialNres'] = 0.

        #calcul des Nzigo et Nres prevus nulle part: seulement C qui etait compte (a developper!)


    #def UpdateNorgbalance(self, dNMin):
    #    surfsolref = self.soilSurface()
    #    self.bilanN['cumMinN'].append(sum3(dNMin) / surfsolref *10000)
    #    # autres termes avec microbio et residus a ajouter!!

    #distribue dans differentes fonctions
    #def UpdateNminbalance(self, Lix, dN2O):
    #    surfsolref = self.soilSurface()
    #    self.bilanN['cumLix'].append(Lix / surfsolref *10000)
    #    self.bilanN['cumN2O'].append(dN2O/ surfsolref *10000)
    #    # autres termes avec plantes...

    def CloseNbalance(self, print_=1):
        surfsolref = self.soilSurface()
        #N org humus
        self.bilanN['FinalInertN'] = sum3(self.InertNorg) / surfsolref *10000
        self.bilanN['FinalActiveN'] = sum3(self.Norg-self.InertNorg) / surfsolref *10000


        #residus
        try:
            self.bilanN['finalNZygo'] = sum(self.ls_NBio())/ surfsolref *10000
            self.bilanN['finalNres'] = sum(self.ls_NRES())/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanN['finalNZygo'] = 0.
            self.bilanN['finalNres'] = 0.

        self.bilanN['ResidueMinNtot'] = sum(self.bilanN['cumNRes1']) + sum(self.bilanN['cumNRes2']) + sum(self.bilanN['cumNRes3'])
        try:
            self.bilanN['NminfromNresCum'] = list(map(sum, self.bilanN['NminfromNres']))
        except:
            self.bilanN['NminfromNresCum'] = 0.

        self.bilanN['HumusMinNtot'] = sum(self.bilanN['cumMinN'])
        self.bilanN['MinNtot'] = self.bilanN['ResidueMinNtot'] + self.bilanN['HumusMinNtot'] 
        self.bilanN['InputNtot'] = self.bilanN['intialInertN'] + self.bilanN['intialActiveN'] + self.bilanN['initialNres'] + self.bilanN['initialNZygo']
        self.bilanN['OutputNtot'] = self.bilanN['FinalInertN'] + self.bilanN['FinalActiveN'] + self.bilanN['finalNres'] + self.bilanN['finalNZygo'] + self.bilanN['MinNtot'] 


        #Input Min
        self.bilanN['TotNRain'] = sum(self.bilanN['cumRain'])
        self.bilanN['TotNIrrig'] = sum(self.bilanN['cumIrrig'])
        self.bilanN['TotFertNO3'] = sum(self.bilanN['cumfertNO3'])
        self.bilanN['TotFertNH4'] = sum(self.bilanN['cumfertNH4'])

        self.bilanN['InputNmintot'] = self.bilanN['intialNO3'] + self.bilanN['intialNH4'] + self.bilanN['MinNtot'] + self.bilanN['TotNRain'] + self.bilanN['TotNIrrig'] + self.bilanN['TotFertNO3'] + self.bilanN['TotFertNH4']
        
        #Output Min
        self.bilanN['FinalNO3'] = sum3(self.m_NO3) / surfsolref *10000
        self.bilanN['FinalNH4'] = sum3(self.m_NH4) / surfsolref *10000
        self.bilanN['Lixtot'] = sum(self.bilanN['cumLix'])
        self.bilanN['N2Otot'] = sum(self.bilanN['cumN2O']) #!! manque denitrif!
        self.bilanN['TotUptPlt'] = sum(self.bilanN['cumUptakePlt'])
        self.bilanN['OutputNmintot'] = self.bilanN['FinalNO3'] + self.bilanN['FinalNH4'] + self.bilanN['Lixtot'] + self.bilanN['N2Otot'] + self.bilanN['TotUptPlt']

        if print_==1:
            self.PrintNbalance()


    def PrintNbalance(self):
        bilanN= self.bilanN
        #Norg
        print ("")
        print ("Organic N Balance Input (kg N.ha-1)\t\t\t Organic N Balance Output (kg N.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Active Humified Pool:\t {0:8.1f}\t\t\t Active Humified Pool:\t {1:8.1f}".format(bilanN['intialActiveN'], bilanN['FinalActiveN'])))
        print(("Inert Humified Pool:\t {0:8.1f}\t\t\t Inert Humified Pool:\t {1:8.1f}".format(bilanN['intialInertN'], bilanN['FinalInertN'])))
        print(("Zymogeneous Bio Pool:\t {0:8.1f}\t\t\t Zymogeneous Bio Pool:\t {1:8.1f}".format(bilanN['initialNZygo'], bilanN['finalNZygo'])))
        print(("Added organic matter:\t {0:8.1f}\t\t\t Added organic matter:\t {1:8.1f}".format(bilanN['initialNres'], bilanN['finalNres'])))
        print(("                            \t\t\t\t Mineralisation:\t\t {0:8.1f}".format(bilanN['MinNtot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t\t {1:8.1f}".format(bilanN['InputNtot'], bilanN['OutputNtot'])))
        print ("")

        #Nmin
        print ("")
        print ("Mineral N Balance Input (kg N.ha-1)\t\t\t Mineral N Balance Output (kg N.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Initial soil NO3 :\t\t {0:8.1f}\t\t\t Final soil NO3:\t {1:8.1f}".format(bilanN['intialNO3'], bilanN['FinalNO3'])))
        print(("Initial soil NH4:\t\t {0:8.1f}\t\t\t Final soil NH4:\t {1:8.1f}".format(bilanN['intialNH4'], bilanN['FinalNH4'])))
        print(("Humus Mineralisation:\t {0:8.1f}\t\t\t Leaching:\t\t\t {1:8.1f}".format(bilanN['HumusMinNtot'], bilanN['Lixtot'])))
        print(("Resid. Mineralisation:\t {0:8.1f}\t\t\t N2O:\t\t\t\t {1:8.1f}".format(bilanN['ResidueMinNtot'], bilanN['N2Otot'])))
        print(("Rain:\t\t\t\t\t {0:8.1f}\t\t\t Uptake plant:\t\t {1:8.1f}".format(bilanN['TotNRain'], bilanN['TotUptPlt'])))
        print(("Irrigation:\t\t\t\t {0:8.1f}\t\t\t ".format(bilanN['TotNIrrig'])))
        print(("Fertilizers NO3:\t\t {0:8.1f}\t\t\t ".format(bilanN['TotFertNO3'])))
        print(("Fertilizers NH4:\t\t {0:8.1f}\t\t\t ".format(bilanN['TotFertNH4'])))


        #print ("Rain:\t\t\t\t\t {0:8.1f}\t\t\t ".format(bilanN['TotNRain']))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t {1:8.1f}".format(bilanN['InputNmintot'], bilanN['OutputNmintot'])))
        print ("")

    def ConcNO3(self):
        """ calculation of the nitrate concentration (kg N.mm-1) by voxel - 8.34 p 160 """
        #mm d'eau libre? (eau liee retiree -  pas dans les eaux de drainage)?
        #non - rq: dans devienne-baret (2000): utilise toute l'eau du sol pour calcul de concentration
        return self.m_NO3 / self.tsw_t#(S.tsw_t - S.m_QH20min + 0.00000001)

    def ConcN(self):
        """ calculation of the molar concentration of mineral nitrogen (micromole N.L-1) by voxel - 8.36 p 161 """
        #L d'eau libre (eau liee retiree - car pas dans les eaux de drainage)
        MMA = 142.85 #mole d'N.kg-1 (g.mole-1)
        moleN = (self.m_NO3 + self.m_NH4)/(MMA)*10**6 #micromole d'N
        return moleN / (self.tsw_t * self.m_vox_surf) * 10000   #remis pour conc sur 1ha pour coller au parametrage de sTICS

    def stepNuptakePlt(self, par, paramp=[{}], ls_lrac=None, ls_mWaterUptakePlt=None, ls_demandeN=None):
        """ calculation of actual N uptake by plant - if no plants (baresoil) -> let None in ls_rac,ls_mWaterUptakePlt, ls_demandeN """

        if ls_lrac is None or ls_mWaterUptakePlt is None or ls_demandeN is None: #si pas de plante (au moins fournir le paramp qui donne un nbre de plante)
            ActUpNtot = self.m_1*0.
            ls_Act_Nuptake_plt = [self.m_1*0.]*len(paramp)
            ls_DQ_N = [1.]*len(paramp)
            idmin = self.m_1*-1.
        else: #si plante
            PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt(self, par, paramp, ls_lrac, ls_mWaterUptakePlt)
            ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt(self, ls_Pot_Nuptake_plt, ls_demandeN)

        # retire les nitrates et ammomium rellement preleves du sol
        frac_NO3 =  self.m_NO3 / (self.m_NO3 + self.m_NH4 + 10**-15)
        self.m_NO3 = self.m_NO3 - frac_NO3*ActUpNtot
        self.m_NH4 = self.m_NH4 - (1. - frac_NO3)*ActUpNtot
        #bilan
        self.bilanN['cumUptakePlt'].append(ActUpNtot/self.soilSurface() *10000)
        self.bilanN['azomes'].append((sum(self.m_NO3)+sum(self.m_NH4))/self.soilSurface() *10000)

        return ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin


def VdistribResidues(SoilObj, Amount, Vprop, Wcontent=0.8,Ccontent=0.42):
    """initialise CRES (Amount of decompasble C in the residue) for a given amount/gratient """
    #Amount (T fresh matter.ha-1)
    #Wcontent (water content; proportion)
    #Ccontent (carbon content; propostion)
    #Vprop : list of proportion per horizon (distribution assumed homogeneous for a given horizon) -> same number of elements than dz
    Q = (Amount/10000.)*SoilObj.soilSurface()*1000.  #kg of fresh matter
    QC = Q*(1-Wcontent)*Ccontent #kg of C
    nbcase_layer = len(SoilObj.dxyz[0])*len(SoilObj.dxyz[1])

    m_CRES = []
    for z in range(len(SoilObj.dxyz[2])):
        v = []
        for x in range(len(SoilObj.dxyz[0])):
            vv = []
            for y in range(len(SoilObj.dxyz[1])):
                vv.append(Vprop[z]*QC/nbcase_layer)

            v.append(vv)
        m_CRES.append(v)

    m_CRES = array(m_CRES)
    return m_CRES






def VABSvox(paramp, CONCN):
    """ Specific absorbtion capacity by roots (micromole N.h-1.cm-1 root) eq. 8.36 p 161 """
    """ paramp = plant parameters """
    HATS = paramp['Vmax1']*CONCN / (paramp['Kmax1'] + CONCN)
    LATS = paramp['Vmax2']*CONCN / (paramp['Kmax2'] + CONCN)
    return HATS+LATS

def FLUXRACs(paramp, SN, ls_lrac):
    """ list of potential active uptake rate (kg N. voxel-1, day-1) per root system per voxel - fluxTot=total uptake potential by all roots ; 
    ls_frac_fluxrac = fraction de demande totale par voxel par systeme racinaire - eq. 8.37 p 161 """
    #ls_lrac = ls_ roots (en m) -> convert en cm (*100)

    MMA = 142.85 #mole d'N.kg-1 (g.mole-1)
    ls_flux_rac = []
    fluxTot = SN.m_1*0.
    for i in range(len(ls_lrac)):
        VABS = VABSvox(paramp[i], SN.ConcN())
        flux = VABS * ls_lrac[i] * 100. * 24. / (MMA * 10**6) #calcule car retombe pas sur coeff 33.6 indique p 161
        ls_flux_rac.append(flux)
        fluxTot = fluxTot + flux

    #calcul des fractions par racine
    ls_frac_fluxrac = []
    for i in range(len(ls_flux_rac)):
        ls_frac_fluxrac.append(ls_flux_rac[i] / (fluxTot + 10**-15) )#10**-15 pour eviter division par zero

    return fluxTot, ls_frac_fluxrac
    #renvoie demande flux total et fraction de demande par voxel par systeme racinaire?
    #passer les densite de racines (i.e. ls_roots) plutot que les longueurs?
    


def Convective_Nflux(SN, ls_mWaterUptakePlt):
    """ Eq. 8_34 p 160 - potential passive NO3 uptake with water flow"""
    ls_conv = []
    convtot = SN.m_1*0.
    for i in range(len(ls_mWaterUptakePlt)):
        flux_c = ls_mWaterUptakePlt[i] * SN.ConcNO3()
        ls_conv.append(flux_c)
        convtot = convtot+flux_c

    return convtot, ls_conv
    #renvoyer liste de fraction?

def Diffusive_Nflux(SN, parSN, ls_lrac):
    """ Eq. 8.35 p 161 """
    #ls_lrac = ls_ roots (en m) -> convert en cm (*100)

    ls_diff = []
    difftot = SN.m_1*0.

    FH = (SN.tsw_t - SN.m_QH20wp) / (SN.m_QH20fc - SN.m_QH20wp)
    filtre = FH>0
    filtre = filtre * 1. #remplace valeur negatives par zero
    DIFE = parSN['DIFNg'] * FH * filtre
    coeff = 4*pi**0.5
    for i in range(len(ls_lrac)):
        draci = (ls_lrac[i]*100.) / (SN.m_soil_vol*1000000.) #root density (cm.cm-3)
        flux_d = coeff * DIFE * (SN.m_NO3 + SN.m_NH4) *sqrt(draci)
        ls_diff.append(flux_d)
        difftot = difftot + flux_d

    return difftot, ls_diff
    #a priori ls_diff pas utilise



def Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt):
    """ Eq. 8.33 p 160 """
    convtot, ls_conv = Convective_Nflux(SN, ls_mWaterUptakePlt)
    difftot, ls_diff = Diffusive_Nflux(SN, parSN, ls_lrac)
    ls_suptot = []
    for i in range(len(ls_conv)):
        ls_suptot.append(ls_conv[i] + ls_diff[i])

    return convtot+difftot, ls_suptot
    #a priori ls_suptot pas utilise



def Potential_NuptakeTot(SN, parSN, paramp, ls_lrac, ls_mWaterUptakePlt):
    """ Eq. 8.38 """
    # philosophie: 
    # assimilation de l'azote est active, via equation de Devienne (HATS, LATS)-> FLUXRACs( qui est mal nomme!)
    # transport actif peut etre limitee par le flux d'azote a la racine et mobilite de l'N ds le sol (convectif + diffusif = soilNsupply) -> prend le min des 2
    # in fine, n'est rellement preleve que l'N necessaire a la croissance (actual_N uptake) -> borne

    Ntot = (SN.m_NO3 + SN.m_NH4)*SN.m_obstarac
    supNtot, ls_supNtot =  Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt)
    fluxTot, ls_frac_fluxrac = FLUXRACs(paramp, SN, ls_lrac)
    ActUptakeN = SN.m_1*0.
    idmin = SN.m_1*0 #code des facteur limitant uptake N (0=Ntot available, 1=passive soil supply, 2=Active root uptake)

    for z in range(len(SN.dxyz[2])):
        for x in range(len(SN.dxyz[0])):
            for y in range(len(SN.dxyz[1])):
                lslim = [Ntot[z][x][y], supNtot[z][x][y], fluxTot[z][x][y]] #mini des 3 sources d'N
                minN = min(lslim)
                id_min = lslim.index(minN)
                ActUptakeN[z][x][y] = minN
                idmin[z][x][y] = id_min

    return ActUptakeN, idmin, ls_frac_fluxrac
    # -> ajuste par voxel selon offre du sol et demande des plantes


def Distrib_Potential_Nuptake_Plt(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt):#ls_lrac, ls_mWaterUptakePlt):
    """ Eq. 8.38 - distibue uptake entre NO3/NH4 et entre les differentes plantes"""
    epsilon = 10e-12
    PotUpNtot, idmin, ls_frac_fluxrac = Potential_NuptakeTot(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt)

    #distribution de l'uptakeN par plante
    ##ls_rac_tot = SN.m_1*0.
    ##for rt in ls_lrac: 
    ##   ls_rac_tot += rt

    ls_Pot_Nuptake_plt = []
    ##for rt in ls_lrac: 
    ##    frac_rac_tot = rt / (ls_rac_tot +epsilon)
    ##    ls_Pot_Nuptake_plt.append( frac_rac_tot * PotUpNtot)    

    for rt in ls_frac_fluxrac:
        ls_Pot_Nuptake_plt.append( rt * PotUpNtot)    


    return PotUpNtot, ls_Pot_Nuptake_plt, idmin
    #rq: distribution entre plante se fait en fonction de densite relative de longueur uniquementif relatifs?
    #tenir compte des parametre de transport act: ls_frac_fluxrac de FLUXRACs?? -> OK =fait
    # !! Manque confrontation a demande totale des plantes (8.39): uptake n'excede pas demande des pantes integree sur tout le profil!! -> a reprendre
    # ls_demandeN a fournir!
    #doner l'option de renvoyer zero si paramp, ls_lrac, ls_mWaterUptakePlt sont a None (pour faire tourner en sol nu facilement)


def Actual_Nuptake_plt(SN, ls_Pot_Nuptake_plt, ls_demandeN):
    """ Eq. 8.39, 8.40 p162 """
    #caclule de demande sur offre par plante! -> ls_DQN
    ls_DQ_N = []
    for i in range(len(ls_demandeN)):
        DQ = min(ls_demandeN[i] / sum(ls_Pot_Nuptake_plt[i]), 1.) #plafonne a ratio a 1: preleve uniquement a hauteur de demande de la plante!
        ls_DQ_N.append(DQ)

    #calcul d'un actual uptake par plante et recalcul total
    ls_Act_Nuptake_plt = []
    ActUpNtot = SN.m_1*0.
    for i in range(len(ls_Pot_Nuptake_plt)):
        Act_Nuptake_i = ls_Pot_Nuptake_plt[i] * ls_DQ_N[i]
        ls_Act_Nuptake_plt.append(Act_Nuptake_i)
        ActUpNtot = ActUpNtot + Act_Nuptake_i

    ### retire les nitrates et ammomium rellement preleves du sol
    ##frac_NO3 =  SN.m_NO3 / (SN.m_NO3 + SN.m_NH4 + 10**-15)
    ##SN.m_NO3 = SN.m_NO3 - frac_NO3*ActUpNtot
    ##SN.m_NH4 = SN.m_NH4 - (1. - frac_NO3)*ActUpNtot
    ###bilan
    ##SN.bilanN['cumUptakePlt'].append(ActUpNtot/SN.soilSurface() *10000)

    return ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N


def critN (MS, a=4.8, b=-0.33):
    """ courbe critique de dilution de l'N - marche aussi pour array"""
    # MS = array od MS values (T.ha-1)
    vals = a*MS**b #en %
    if vals.size>1:
        for i in range(len(vals)): vals[i]=min(a, vals[i])
    else:
        vals = min(a, vals)
    return vals 


def demandeNdefaut(MSp,dMSp,Npc, surfsolref, a=4.8, b=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    NcritTot = critN (MStot, a, b)#N crit de MS new
    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN


#surfsolref = 0.05
#MSp = array([1.,1.2, 2.])
#dMSp = array([0.1,0.15,0.2])
#Npc = array([4., 3., 2.])
#demandeNdefaut(MSp,dMSp,Npc, surfsolref)

def demandeNdefaut2(MSp,dMSp,Npc, surfsolref, a=4.8, b1=-0.1 ,b2=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    if MStot>=1.:
        NcritTot = a*MStot**b2#critN (MStot, a, b2)#N crit de MS new dense
    else:
        NcritTot = a*MStot**b1#critN (MStot, a, b1)#N crit de MS new isole

    #filtre NcritTot
    NcritTot[NcritTot>9.]=9.#gN.plant-1

    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN



def demandeNroot(MSpiv,dMSpiv,Npcpiv, surfsolref, Noptpiv):
    """ demande N pour parties racinaire - suppose N critique constant - s'applique aux racines et aux pivots """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QNpiv = MSpiv * Npcpiv/100. #gN.plant-1
    PotNpiv = (MSpiv + dMSpiv) * Noptpiv/100. #gN.plant-1
    ls_demandeN = PotNpiv-QNpiv
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN

#def saveBilanN(bilanN, path):
#
#    dicout= {}
#    dicout['NRain'] = bilanN['cumRain']
#    dicout['NIrrig'] = bilanN['cumIrrig']
#    dicout['fertNO3'] = bilanN['cumfertNO3']
#    dicout['fertNH4'] = bilanN['cumfertNH4']
#    dicout['HumusNMin'] = bilanN['cumMinN']
#    dicout['Res1'] = bilanN['cumNRes1']
#    dicout['Res2'] = bilanN['cumNRes2']
#    dicout['Res3'] = bilanN['cumNRes3']
#    dicout['ResidueMinN'] = bilanN['cumNRes1']+bilanN['cumNRes2']+bilanN['cumNRes3']
#    dicout['Lix'] = bilanN['cumLix']
#    dicout['N2O'] = bilanN['cumN2O']
#    dicout['UptPlt'] = bilanN['cumUptakePlt']



