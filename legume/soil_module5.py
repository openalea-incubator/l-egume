import IOtable
from scipy import *
from copy import deepcopy
#from openalea.plantgl.all import *
#import sys
#path_ = r'H:\devel\grassland\grassland\luzerne'
#sys.path.insert(0, path_)
#from Obj3Dutils import *


##changements version 5
#ajout d'un patern en entree pour initialisation + stocke info surfsolref (facilite couplage racine dans impulse)
# ajout d'une entree obstarac (matrice 2D de profondeur d'obtacle en m ; None par defaut) et d'une matrice d'effet sur la disponibilite des ressource (m_obstrac)
# mise a jour d calcul de FTSW pour tenir compte de m_obstrac
# retire des ; dans le calcul de transpiration
# defini une fonction test_uni() pour faire un test uitaire
# retire dependence a opealea: fonctions de visu dans un module supplementaire:soil_modulevisu1.py

#A faire: retirer


## changements version 4
# distinguer les entree irrigation et rain pour les comptabiliser dans le bilan -> OK
# reprendre distrib_water_uptakeNC pour sortir les uptake par plante -> besoin pour les prelv d'N par plante -> OK



#a faire: map_Evap0 devrait etre une entree (pas evapo_tot) -> encore calcule!
#a faire: k_eq utilise pour evaporation: fixe = a sortir en parametre!
# pour eviter effet de marche de ftsw qd racine passe ds un nveau compartiment: considerer fraction du volume comme pour ZESX?
# gestion du drainage: 2 passe des precipitation -> simplifier (1 passe seulement suffit?) pour gere infiltration des nitrates?
# evapo_tot -> renvoie la somme (pas la moyenne)

# fonctions diverses
def sum3(mat):
    return mat.sum()#sum(sum(sum(mat)))
    ## pour gerer corectement les sommes totale un peu plus compliquee en muldim...


def sum_mat(ls_m):
    #pour gerer correctement les somme de matrices similaires dans une liste (avant sum le faisait bien, mais derniere version bug???)
    s = ls_m[0]
    if len(ls_m)>1:
        for i in range(1,len(ls_m)):
            s=s+ls_m[i]

    return s

def mask(mat, tresh=0.):
    """ retourne matrice masque de 0  et 1
    utile pour root/no roots; asw/no asw """
    #m = deepcopy(mat) #pas utile->mat cree nouvelle matrice
    return where(mat > tresh, 1, 0)
    #m = m*1. #pour convertir en float

#mask(R2)
#mask(asw_t)



##diverses fonction parametrage sol

def tetavol_pF_curve(par_sol, pF_):
    """ empirical relationship between volumic water content and pF (=ln(abs(h))where h is soil moisture succion )
    Driessen, 1986 cited by Penning de Cries et al., 1989
    WCST = volumic water content at saturation
    gamma = empirical parameter to account for bulk density and texture"""
    teta = float(par_sol['WCST'])*exp(-float(par_sol['gamma_theo'])*pF_*pF_)
    return teta
    #tetavol_pF_curve(par_sol, 4.2)

def default_tetaref(par_sol, fc=2., wp=4.2, ad=7.):
    par_sol['teta_sat'] = tetavol_pF_curve(par_sol, 0.)
    par_sol['teta_fc'] = tetavol_pF_curve(par_sol, fc)
    par_sol['teta_wp'] = tetavol_pF_curve(par_sol, wp)
    par_sol['teta_ad'] = tetavol_pF_curve(par_sol, ad)
    return par_sol
    #p = default_tetaref(par_sol)

def pF(h):
    return log10(abs(h))





class Soil(object):
    def __init__(self, par_sol, soil_number = [13]*10, dxyz = [[1.], [1.], [0.2]*10], vDA=[1.2]*10, ZESX=0.3, CFES=1., obstarac=None, pattern8=[[0,0],[100.,100.]]):
        """  
        dxyz : dimensions des voxels selon x, y, z (m)
        m_soil_vol : volume des voxels (m3)
        m_DA : densite apparente du sol par voxel (g.cm-3)
        m_soil_number: index du type de sol par voxel(lien avec dictionnaire des parametres hydraulique du sol) 
        m_vox_surf: surface de sol par voxel (m2)
        soilSurface(): surface de sol (m2)

        asw_t : quantite d'eau libre pour plante dans le voxel au temps t (mm)
        tsw_t : quantite d'eau totale dans le voxel au temps t (mm)
        ftsw_t : fraction d'eau libre disponible pour les plantes (Transpirable) par voxel (0-1, sans dimension)
        Kzi() : fraction d'eau disponible pour evaporation (Evaporable) par voxel (0-1, sans dimension)
        HRv() : humidite volumique du sol (%)
        HRp() : humidite ponderale du sol (%)

        m_QH20sat : quantite d'eau du sol a saturation par voxel (mm)
        m_QH20fc : quantite d'eau du sol a la capacite au champs par voxel (mm)
        m_QH20max : quantite d'eau libre dispo pour la plante  a la capacite au champs par voxel (mm)(field capacipi-wilting point)
        m_QH20wp : quantite d'eau du sol restante au wilting point (non accessible a la plante) par voxel (mm)
        m_QH20min : quantite d'eau du sol minimum a areet de l'evaporation (sol s'aseche pas en dessous) (mm)

        m_teta_sat : humidite volumique du sol a saturation par voxel (fraction 0-1)
        m_teta_fc : humidite volumique du sol a la capacite au champs par voxel  (fraction 0-1)
        m_teta_wp : humidite volumique du sol au wilting point par voxel  (fraction 0-1)
        m_teta_ad : humidite volumique minimum du sol sec  (fraction 0-1)

        m_frac_evapZ : distribution verticale relative de evaporation a la capacite au champ (qd Kzi=1) (0-1, sans dimension)
        m_corr_vol_evap : complement a la fraction de volume effectivement accessible a evaporation par voxel (0-1, sans dimension)

        ## SORTIES:
        evapo_tot : quantite d'eau evaporee totale au temps t (mm)
        m_evap : quantite d'eau evaporee au temps t par voxel (mm)
        m_transpi : quantite d'eau transpiree au temps t par voxel (mm) - cumul des n plantes par voxel
        ls_drainage : quantite d'eau drainee au temps t par voxel (mm) - liste de matrice par horizon
        ls_transp : quantite d'eau transpiree au temps t par plante (mm) - liste de cumul par plante
        ls_ftsw : liste des fraction d'eau disponible pour la transpiration par plante (fraction 0-1) - tient compte des distributions des racines

        """

        self.dxyz = dxyz
        size = list(map(len, dxyz))
        self.size = size[2:]+size[0:2] #z,x,y
        self.pattern = pattern8 #cm
        self.origin = array(pattern8[0]+[0.])/100. #m
        Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
        largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
        self.surfsolref = Lsol * largsol  # m2

        self.m_1 = []
        self.m_soil_vox = []#voxel positions
        self.m_soil_vol = []#voxel volume
        self.m_soil_number = []#strate de type de sol
        self.m_DA = []
        self.m_vox_surf = []
        for z in range(len(dxyz[2])):
            v, vvox, v2, v3, vd, vs = [], [], [], [], [], []
            for x in range(len(dxyz[0])):
                vv, vvvox, vv2, vv3, vvDA, vvs = [], [], [], [], [], []
                for y in range(len(dxyz[1])):
                    vv.append(1)
                    vol = dxyz[0][x] * dxyz[1][y] * dxyz[2][z]
                    dists = [sum(dxyz[0][0:x]), sum(dxyz[1][0:y]), sum(dxyz[2][0:z])]
                    vvvox.append(array(dists))
                    vv2.append(vol)
                    vv3.append(soil_number[z])
                    vvDA.append(vDA[z])
                    vvs.append(dxyz[0][x] * dxyz[1][y])

                v.append(vv)
                vvox.append(vvvox)
                v2.append(vv2)
                v3.append(vv3)
                vd.append(vvDA)
                vs.append(vvs)

            self.m_1.append(v)
            self.m_soil_vox.append(vvox)
            self.m_soil_vol.append(v2)
            self.m_soil_number.append(v3)
            self.m_DA.append(vd)
            self.m_vox_surf.append(vs)

        self.m_1, self.m_soil_vox, self.m_soil_vol, self.m_soil_number, self.m_DA, self.m_vox_surf = array(self.m_1), array(self.m_soil_vox), array(self.m_soil_vol), array(self.m_soil_number), array(self.m_DA), array(self.m_vox_surf)
        self.m_frac_evapZ = self.distrib_frac_evap(ZESX, CFES)
        self.m_corr_vol_evap = self.corr_vol_voxel_evap(ZESX)
        self.compute_teta_lim(par_sol)
        self.obstarac = obstarac
        self.m_obstarac  = self.set_mat_obstarac()
        self.init_asw(HRp_init=None)  # initialise par defaut a la capacite au champ...
        # self.update_ftsw()

    def build_teta_m(self, par_sol, key='teta_fc'):
        m_teta = []#m_teta = matrice de teneur en eau volumique
        for z in range(len(self.m_soil_number)):
            v = []
            for x in range(len(self.m_soil_number[z])):
                vv = []
                for y in range(len(self.m_soil_number[z][x])):
                    vv.append(par_sol[str(self.m_soil_number[z][x][y])][key])

                v.append(vv)
            m_teta.append(v)
        return array(m_teta)


    def compute_teta_lim(self, par_sol):
        self.m_teta_fc = self.build_teta_m(par_sol, key='teta_fc')
        self.m_teta_ad = self.build_teta_m(par_sol, key='teta_ad')
        self.m_teta_wp = self.build_teta_m(par_sol, key='teta_wp')
        self.m_teta_sat = self.build_teta_m(par_sol, key='teta_sat')
        self.m_QH20fc = multiply(self.m_soil_vol, self.m_teta_fc)*1000. #quantite totale max en mm d'eau dans le sol a la capacite au champ
        self.m_QH20max = multiply(self.m_soil_vol, self.m_teta_fc - self.m_teta_wp)*1000. #quantite max en mm dispo pour la plante (field capacipi-wilting point)
        self.m_QH20wp = multiply(self.m_soil_vol, self.m_teta_wp)*1000. #quantite eau restante au wp (non accessible a la plante)
        self.m_QH20min = multiply(self.m_soil_vol, self.m_teta_ad)*1000. #quantite min en mm dispo pour evaporation (sol s'aseche pas en dessous)
        self.m_QH20sat = multiply(self.m_soil_vol, self.m_teta_sat)*1000. #quantite totale max d'eau dans le sol a saturation 
        

    #def init_asw(self, frac=1.):
    #    """ par defaut a la capacite au champ """
    #    #avant, faut avoir lance compute_teta_lim(par_sol)
    #    F = self.m_1*frac
    #    self.tsw_t = multiply(F, self.m_QH20fc) #cree matrice tsw_t
    #    self.asw_t = multiply(F, self.m_QH20max) #cree matrice asw_t

    def init_asw(self, HRp_init=None):
        """ par defaut a la capacite au champ """
        'HRp_init : a 1D vector/list of initial HRp (%) in the vertical z axis'
        if HRp_init==None:
            #avant, faut avoir lance compute_teta_lim(par_sol)
            F = self.m_1*1.
            self.tsw_t = multiply(F, self.m_QH20fc) #cree matrice tsw_t
            self.asw_t = multiply(F, self.m_QH20max) #cree matrice asw_t
            self.update_ftsw()
            self.OpenWbalance() #initialise W balance
        else:#initialisation avec un vecteur d'humidite ponderales
            for z in range(len(self.dxyz[2])):
                self.tsw_t[z,:,:] = self.m_1[z,:,:]*HRp_init[z]*(10*self.m_DA[z,:,:])*self.m_soil_vol[z,:,:]

            self.update_asw()
            self.update_ftsw()
            self.OpenWbalance()#initialise W balance

    def update_asw(self):
        self.asw_t = self.tsw_t - self.m_QH20wp

    def set_mat_obstarac(self):
        """ definit une matrice de proportion de voxel concerne par obstrarac"""
        if self.obstarac is None or type(self.obstarac)!=type(array([0.])):#pas d'obstarac
            mat_obstarac = 1. * self.m_1
        else: #obstrac = matrice 2D de valeurs negatives (m)
            ls_prof = [0.]
            for i in self.dxyz[2]:
                ls_prof.append(ls_prof[-1] - i)

            mat_obstarac = 1. * self.m_1
            for x in range(len(self.dxyz[0])):
                for y in range(len(self.dxyz[1])):
                    obstXY = self.obstarac[x, y]
                    for z in range(len(self.dxyz[2])):
                        if ls_prof[z] > obstXY >= ls_prof[z + 1]:
                            ratio = (obstXY - ls_prof[z]) / (ls_prof[z + 1] - ls_prof[z])
                            mat_obstarac[z, x, y] = ratio
                        elif obstXY > ls_prof[z + 1]:
                            mat_obstarac[z, x, y] = 0.
        return mat_obstarac

    def update_ftsw(self):
        ftsw = self.asw_t / self.m_QH20max
        negs = ftsw >0.
        self.ftsw_t = negs * 1. *ftsw * self.m_obstarac #pour eviter valeurs negatives

    def HRv(self):
        """ compute relative soil humidity - Humidite volumique(%, 100 * vol eau par volume de sol)"""
        return divide(self.tsw_t, self.m_QH20fc)*100.
        # peut etre 100 pour cent?? -> OK par rapport aux parametetre de reponse de mineralisation a l'humidite du sol (proportion de capciate au champ)

    def HRp(self):
        """ compute relative soil humidity - Humidite ponderale (%, 100 * g H2O.g sol-1)"""
        return (self.tsw_t / self.m_soil_vol)/ (self.m_DA * 10.)

    def Kzi(self):
        """ coefficients K(Z,I) de l'equation 7.6 p 129 - proportion de l'eau dispo pour evaporation """
        return (self.tsw_t - self.m_QH20min) / (self.m_QH20fc - self.m_QH20min)

    def distrib_frac_evap(self, ZESX=0.3, CFES=1.):
        """ distribution verticale des fractions de l'evaporation en fonction de ZESX et CFES pour Kzi=1 p- integration sur dz=1cm - Eq. 7.6 p 129"""
        #calcul des ponderations par cm
        nbstratesZ = len(self.dxyz[2])
        v = [self.dxyz[2][0]]
        for i in range(1, nbstratesZ):
            v.append(v[-1]+self.dxyz[2][i])

        nbcm = int(v[-1]*100)
        calc = []
        for i in range(nbcm):
            if float(i)/100. <= ZESX:
                calc.append((ZESX-float(i)/100.)**CFES)
            else:
                calc.append(0.)

        calc = array(calc)/sum(calc)
        array(list(range(nbcm)))

        #repartition dans les n strates du sol selon dxyz
        dics = {}
        for i in range(nbstratesZ):
            dics[i] = array(list(range(nbcm)))/100./v[i] >= 1.
            dics[i] = dics[i]*1

        strate = array([0]*len(calc))
        for k in list(dics.keys()):
            strate += dics[k]

        vsum = [0.]*nbstratesZ
        for i in range(len(calc)):
            vsum[strate[i]] += calc[i]

        #creation de la matrice
        m = deepcopy(self.m_1)*1.
        for z in range(nbstratesZ):
            m[z,:,:] = m[z,:,:]*vsum[z]

        return m

    def corr_vol_voxel_evap(self, ZESX=0.3):
        """ correction pour tenir compe du fait que voxel de transition = toute l'eau n'est pas dispo (sous ZESX) -> a stocker """
        v = [self.dxyz[2][0]]
        for i in range(1, len(self.dxyz[2])):
            v.append(v[-1]+self.dxyz[2][i])

        epsilon = 0.000001 #pour gerer pb d'arrondi au limite des voxels
        test = array(v)/ZESX -1.  <= epsilon
        test = test *1.
        restes = (array(v)*test)%ZESX
        test = test.tolist()
        id_transition = test.index(0.)
        test = array(test)*1.
        test[id_transition] = (self.dxyz[2][id_transition] - restes[id_transition]) / self.dxyz[2][id_transition]
        m_cor_vol = deepcopy(self.m_1)*1.
        for z in range(m_cor_vol.shape[0]): m_cor_vol[z,:,:] = m_cor_vol[z,:,:]*test[z]

        comp = (1.-m_cor_vol)
        vide = comp <1.
        comp = comp * vide*1.

        return comp#m_cor_vol #renvoie le complement finalement, car utilise dans les calculs

    def distrib_evapSTICS(self, map_Evap0):
        """ map_Evap0 = liste des Evap par colone de voxel (matrice x, y) """
        # remplace sans bug check_soil_evap et distrib_evapNC
        coeffs = self.Kzi()*self.m_frac_evapZ
        coeffs1 = deepcopy(self.m_1)*1.
        for z in range(coeffs1.shape[0]): coeffs1[z,:,:] = sum(coeffs, axis=0) #etend en 3D la somme

        coeffsOK = coeffs/(coeffs1 + 0.0000000001) #coeffs equation 7.6, livre STICS

        demande = deepcopy(self.m_1)*1.
        for z in range(demande.shape[0]): demande[z,:,:] = map_Evap0 #etend en 3D la demande evapo

        demandeEvap = demande * coeffsOK
        dispo =  self.tsw_t - self.m_QH20min - self.m_corr_vol_evap*(self.m_QH20fc - self.m_QH20min)
        diff = dispo - demandeEvap
        negs = diff<=0.
        demandeEvapOK = demandeEvap + negs*diff
        map_EvapOK = sum(demandeEvapOK,axis=0)

        return map_EvapOK, demandeEvapOK
        #S.tsw_t - demandeEvapOK

    def soilSurface(self):
        """ compute soil surface (m2) """
        return sum(self.dxyz[0])*sum(self.dxyz[1])

    def stepWBmc(self, Et0, ls_roots, ls_epsi, Rain, Irrig, previous_state, ZESX=0.3, leafAlbedo=0.15, U=5., b=0.63, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1):
        """ calcul du bilan hydrique journalier """

        Precip = Rain + Irrig

        ## transpi 
        ls_roots_eff = effective_root_lengths(ls_roots, tresh = treshEffRoots)
        ls_mask = list(map(mask, ls_roots_eff))#liste de masques racine
        masked_aswt = multiply(mask(self.asw_t), self.asw_t) #pour sol, met a zero si teta sous wp
        ls_masked_asw = list(map(multiply, ls_mask, [masked_aswt]*len(ls_mask)))#pour chaque root system, asw(mm) si racine et teta> wp; zero sinon
        ls_masked_tswt = list(map(multiply, ls_mask, [self.m_QH20max]*len(ls_mask)))
        ls_ftsw = list(map(divide, list(map(float, list(map(sum3, ls_masked_asw)))), list(map(float, list(map(sum3, ls_masked_tswt))))))
        ls_transp = Transpi_NC(Et0, ls_epsi, ls_ftsw, leafAlbedo, FTSWThreshold)
        #m_transpi = distrib_water_uptakeNC(self.asw_t, ls_masked_asw, ls_roots_eff, ls_transp)
        m_transpi, ls_m_transpi = distrib_water_uptakeNC(self, ls_masked_asw, ls_roots_eff, ls_transp)
        ls_transp = list(map(sum, ls_m_transpi)) #correct ls_trans with actual transpiration

        ## evaporation
        #evapo_tot, state = soil_EV_1C(Et0, Precip, sum(ls_epsi), previous_state , leafAlbedo, U, b)
        evapo_tot, state = soil_EV_STICS(Et0, Precip, sum(ls_epsi), previous_state , leafAlbedo, U, b)
        map_Evap0 = self.m_1[0]*evapo_tot/float(len(self.dxyz[0])*len(self.dxyz[1])) #map_Evap0 -> devrait etre une entree? (interaction avec bilan radiatif?)

        ## bilan WB
        map_PI = self.m_1[0]*Precip/float(len(self.dxyz[0])*len(self.dxyz[1]))
        # 1) applique pluie et calcule D0 (drainage sans transpi pour saturer)
        ##m, D0 = distrib_PI(self.asw_t, self.m_QH20max, map_PI, self.dxyz,opt)
        m, ls_D0 = distrib_PI(self.tsw_t, self.m_QH20fc, map_PI, self.dxyz,opt)
        D0 = ls_D0[-1]
        # 2) retire evapo et transpi
        self.tsw_t = m 
        map_EvapOK, m_evap = self.distrib_evapSTICS(map_Evap0) # evap tient compte des pluies
        evapo_tot = sum(map_EvapOK)#mean(map_EvapOK) #corrige pour eviter valeur sous mini
        m = m - m_evap - m_transpi

        # 3) remet ce qui a drainer et recalcule drainage residuel apres transpi (=transpire en priorite ce qui est tombe pdt la journee)
        ##m, D1 = distrib_PI(m, self.m_QH20max, D0, self.dxyz)
        m, ls_drainage = distrib_PI(m, self.m_QH20fc, D0, self.dxyz)
        #D1 = ls_drainage[-1] #drainage profond

        self.tsw_t = m
        self.update_asw()
        self.update_ftsw()
        self.UpdateWbalance(Rain, Irrig, evapo_tot, ls_transp, ls_drainage) #Irrigj pas distingue en entree pour le moment!-> a faire

        return ls_transp, evapo_tot, ls_drainage, state,  ls_m_transpi, m_evap, ls_ftsw


    #def plot_soil_properties (self, vals, MaScene=Scene(), col_scale=5):#dxyz, m_soil_vox, asw_t pourraient etre remplacee par objet sol
    #    """ vals = matrice de valeur entre 0 et 1 de propriete de sol a visualiser / e.g S.ftsw_t """
    #    bx = Box(Vector3(1.,1.,1.))
    #    for z in range(len(self.dxyz[2])):
    #        for x in range(len(self.dxyz[0])):
    #            for y in range(len(self.dxyz[1])):
    #                dims = [self.dxyz[0][x]/2., self.dxyz[1][y]/2., self.dxyz[2][z]/2.]
    #                p_ini = self.m_soil_vox[z][x][y]
    #                col = couleur (col_scale, max(0.,vals[z][x][y]))
    #                b = transformation(bx, dims[0], dims[1], dims[2], 0,0,0, p_ini[0], p_ini[1], -p_ini[2]-self.dxyz[2][z]/2.)
    #                MaScene.add(Shape(b, Material(Color3(col[0],col[1],col[2]), transparency=0.65)))
    #
    #    return MaScene
    ## mis dans soil_modulevisu1.py


    def OpenWbalance(self):
        """ Dictionnary for soil Water balance (mm)
        Keys for Daily outputs: 'cumPP', 'cumIrrig', 'cumD', 'cumEV', 'cumTransp'
        Keys for Total Input: 'intialWC', 'Irrigtot','PPtot'
        Keys for Total outputs: 'FinalWC', 'EVtot', 'Tranptot','Drainagetot'
        Keys for totals: 'InputWtot', 'OutputWtot'
        """
        self.bilanW = {}
        self.bilanW['intialWC'] = sum3(self.tsw_t) / self.soilSurface()
        self.bilanW['cumEV'], self.bilanW['cumPP'], self.bilanW['cumIrrig'], self.bilanW['cumD'], self.bilanW['cumTransp'] = [],[],[], [], []
        self.bilanW['TSWt'] = []
        #cumET0 = []
        
    def UpdateWbalance(self, PPj, Irrigj, evapo_tot, ls_transp, Drainage):
        surfsolref = self.soilSurface()
        self.bilanW['cumPP'].append(PPj / surfsolref)
        self.bilanW['cumIrrig'].append(Irrigj / surfsolref)
        self.bilanW['cumTransp'].append(sum(ls_transp) / surfsolref)
        self.bilanW['cumEV'].append(evapo_tot/ surfsolref)
        self.bilanW['cumD'].append(sum(Drainage[-1]) / surfsolref)
        self.bilanW['TSWt'].append(sum3(self.tsw_t) / self.soilSurface())


    def CloseWbalance(self, print_=1):
        #input
        self.bilanW['PPtot'] = sum(self.bilanW['cumPP'])
        self.bilanW['Irrigtot'] = sum(self.bilanW['cumIrrig'])
        self.bilanW['InputWtot'] = self.bilanW['intialWC'] + self.bilanW['PPtot'] + self.bilanW['Irrigtot']
        #output
        self.bilanW['FinalWC'] = sum3(self.tsw_t) / self.soilSurface()
        self.bilanW['EVtot'] = sum(self.bilanW['cumEV'])
        self.bilanW['Tranptot'] = sum(self.bilanW['cumTransp'])
        self.bilanW['Drainagetot'] = sum(self.bilanW['cumD'])
        self.bilanW['OutputWtot'] = self.bilanW['FinalWC'] + self.bilanW['EVtot'] + self.bilanW['Tranptot'] + self.bilanW['Drainagetot']

        #print self.bilanW['OutputWtot'], self.bilanW['InputWtot']
        if print_==1:
            self.PrintWbalance()
        #pourrait le diriger vers un fichier de sortie texte?

    def PrintWbalance(self):
        bilanW = self.bilanW
        print ("")
        print ("Water Balance Input (mm)\t\t\t\t\t Water Balance Output (mm)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Initial Soil Water:\t {0:8.1f}\t\t\t\t Final Soil Water:\t {1:8.1f}".format(bilanW['intialWC'], bilanW['FinalWC'])))
        print(("Precipitation:\t\t {0:8.1f}\t\t\t\t Transpiration:\t\t {1:8.1f}".format(bilanW['PPtot'], bilanW['Tranptot'])))
        print(("Irrigation:\t\t\t {0:8.1f}\t\t\t\t Evaporation:\t\t {1:8.1f}".format(bilanW['Irrigtot'],bilanW['EVtot'])))
        print(("                            \t\t\t\t Deep infiltration:\t {0:8.1f}".format(bilanW['Drainagetot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t {0:8.1f}\t\t\t\t Total:\t\t\t\t {1:8.1f}".format(bilanW['InputWtot'], bilanW['OutputWtot'])))
        print ("")









#ls_transp, evapo_tot, D, state,  m_frac_transpi, m_frac_evap = S.stepWBmc(Et0, ls_roots, ls_epsi, Precip, previous_state, ZESX=0.3, leafAlbedo=0.15, U=5., b=0.63, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)


#def plot_soilWC (dxyz, m_soil_vox, ftsw_t, MaScene=Scene(), col_scale=5):#dxyz, m_soil_vox, asw_t pourraient etre remplacee par objet sol
#    """ """
#    for z in range(len(dxyz[2])):
#        for x in range(len(dxyz[0])):
#            for y in range(len(dxyz[1])):
#                dims = [dxyz[0][x]/2., dxyz[1][y]/2., dxyz[2][z]/2.]
#                p_ini = m_soil_vox[z][x][y] 
#                col = couleur (col_scale, max(0.,ftsw_t[z][x][y]))
#                b = transformation(bx, dims[0], dims[1], dims[2], 0,0,0, p_ini[0], p_ini[1], -p_ini[2]-dxyz[2][z]/2.)
#                #b.setName(luz_label(1, 0, 0, 0, 0))#ajout d'un label canestra a ajuster
#                MaScene.add(Shape(b, Material(Color3(col[0],col[1],col[2]), transparency=0.65)))
#
#    return MaScene
## a passer comme methode classe



## fonction bilan hydrique
# methode "init_asw"
#def ASWtot(m_teta, m_teta_fc, m_teta_wp, m_QH20max):
#    frac = divide(m_teta-m_teta_wp, m_teta_fc-m_teta_wp)
#    return multiply(frac, m_QH20max)


def water_uptakeVox(ls_masked_asw, ls_roots_eff, ls_transp, idx, idy, idz):
    """ STICS manual (p 140) adapte pour liste de root systems"""
    ls_uptake = []
    for p in range(len(ls_roots_eff)):
        if ls_masked_asw[p][idz][idx][idy]>0.:#y a de l'eau dispo
            frac_root = ls_roots_eff[p][idz][idx][idy]/sum3(ls_roots_eff[p])
            frac_asw = ls_masked_asw[p][idz][idx][idy]/sum3(ls_masked_asw[p])
            uptake_vox = (ls_transp[p] * (frac_root+frac_asw))/2.
        else:
            uptake_vox = 0.

        ls_uptake.append(uptake_vox)

    return ls_uptake
    #tester si prend plus... ajuster?

#water_uptakeVox(ls_masked_asw, ls_roots_eff, ls_Et0, 0, 0, 0)


def distrib_water_uptakeNC(S, ls_masked_asw, ls_roots_eff, ls_transp):
    asw_t = deepcopy(S.asw_t)
    upt = S.m_1*0.

    ls_upt = []
    for i in range(len(ls_roots_eff)):
        ls_upt.append(S.m_1*0.)

    for z in range(len(asw_t)):
        for x in range(len(asw_t[z])):
            for y in range(len(asw_t[z][x])):
                upvx = water_uptakeVox(ls_masked_asw, ls_roots_eff, ls_transp, x, y, z)
                upt[z][x][y] = sum(upvx) #gere negatif avec masque en amont
                for plt in range(len(ls_upt)):
                    ls_upt[plt][z][x][y] = upvx[plt]

    return upt, ls_upt
    #renvoyer info par plante?ici distribue slmt entre voxel -> oui!! besoin pour les flux d'N

#upt = distrib_water_uptakeNC(asw_t, ls_masked_asw, ls_roots_eff, ls_Et0)
#sum(distrib_water_uptakeNC(asw_t, ls_masked_asw, ls_roots_eff, ls_Et0)), sum(ls_Et0)


def Transpi_NC(Et0, ls_epsi, ls_FTSW, leafAlbedo=0.15, FTSWThreshold=0.4):
    ls_transp = []
    for i in range(len(ls_epsi)):

        ##Riou
        #coeffTV = ls_epsi[i]/ (1-leafAlbedo)   
        #potentialVineTranspiration = coeffTV * Et0

        ##adaptation eq 7.1 p127 STICS et 7.8 p131 book
        k_eq = 0.7 #LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
        LAIeq = -log(1-ls_epsi[i])/k_eq 
        potentialTranspiration = Et0 * (1 - exp(-(k_eq-0.2)*LAIeq))#complement a potentialSoilEvaporation
        #rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)
        # faut introduire crop coefficient KMAXp cf eq. 7.7 -> pour le moment limite a Et0 (reference gazon)

      
        if (ls_FTSW[i] > FTSWThreshold):#previousTSW/TTSW
            # la transpiration de la plante est faite a son potentielle
            Ks = 1.
        else:
            # regulation
            # la quantite d'eau presente dans le sol n'est pas suffisante pour
            # que la transpiration de plante se fasse a son maximum   
            Ks = ls_FTSW[i]/FTSWThreshold

        ls_transp.append(Ks*potentialTranspiration)

    return ls_transp

#Transpi_NC(2., [0.4], [0.8], leafAlbedo=0.15, FTSWThreshold=0.4)



##gerer sol evapo 
def soil_EV_1C(Et0, Precip, epsi, previous_state = [0.,0.,0.], leafAlbedo=0.15, U=5., b=0.63):
    """ U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
        b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
        meme formulation pour Lebon et al 2003 et STICS"""

    #recup previous state
    previousSES1 = previous_state[0]
    previousSES2 = previous_state[1]
    previousSEP = previous_state[2]

    coeffTV = epsi/ (1.-leafAlbedo)
    potentialSoilEvaporation = (1. - coeffTV) * Et0
    EP = (1. - coeffTV) * Et0
        
    SES1=max(0., previousSES1)
    SES2=max(0.,previousSES2)
    SEP=max(0.,previousSEP)
    
    P = Precip
    P1 = Precip 
    
    if (SES1<U):
        # phase 1
        if (P1>SES1):
            SES1 = 0.
            SES1 = SES1 + EP
        else:
            SES1 = SES1 - P
            SES1 = SES1 + EP
        
        if (SES1>U):
            ER = EP - 0.4*(SES1-U)
            SES2 = 0.6*(SES1-U)
        else:
            ER = EP
    else:
        # phase 2
        if (P1>SES2):
            SEP = 0
            P1 = P1 - SES2
            SES1 = U - P1
            
            if (P1>U):
                SES1 = EP
            else:
                SES1 = SES1 + EP
            
            if (SES1>U):
                ER = EP - 0.4*(SES1-U)
                SES2 = 0.6*(SES1-U)
            else:
                ER = EP
            
        else:
            SES21 = ((2*b*SEP + b**2)**0.5) - b
            SEP = SEP + EP
            SES20 = ((2*b*SEP + b**2)**0.5) - b
            ER = SES20 - SES21
            
            if (P>0):
                ESx = 0.8*P
                
                if (ESx<ER):
                    ESx = ER + P
                
                if (ESx>EP):
                    ESx = EP
                
                ER = ESx
                SES2 = SES2 + ER - P
            else:
                if (ER>EP):
                    ER = EP
                    SES2 = SES2 + ER - P
                else:
                    SES2 = SES2 + ER - P

    state = [SES1, SES2, SEP]
    soilEvaporation = ER 
    return soilEvaporation, state

# pourrait gerer par voxel de surface plutot que globalement (ne pas assumer homogeneite horizontale)


def soil_EV_STICS(Et0, Precip, epsi, previous_state = [0.,0.,0.], leafAlbedo=0.15, U=5., b=0.63):
    """ U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
        b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
        Et0: Evapotranspiration potenrielle de reference (gazon) en mm
        meme formulation pour Lebon et al 2003 et STICS
        SES1: EP cumule depuis la derniere pluie (previous_state[0])
        SES2: EP cumule depuis passage a phase 2  (previous_state[1])  
        ER SoilEvaporation Reel"""

    ##Riou -> donne des valeur negatives pour forts epsi
    #coeffTV = epsi/ (1.-leafAlbedo)
    #EP = (1. - coeffTV) * Et0#potentialSoilEvaporation

    ##adaptation eq 7.1 p127 STICS book
    k_eq = 0.7 #LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
    LAIeq = -log(1-epsi)/k_eq 
    EP = Et0 * exp(-(k_eq-0.2)*LAIeq)#potentialSoilEvaporation Eq. 7.1, p127
    #rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)


    if previous_state[0]+EP <U : #phase 1
        SES1 = previous_state[0]+EP
        SES2 = 0.
        ER=EP
        SEP = 0.
        #mise a jour state!
    elif previous_state[0] <U and previous_state[0]+EP >=U:#transition phase1->2
        fracEP2 = previous_state[0]+EP-U #precipitation evaporee en phase 2
        fracEP1 = EP-fracEP2
        SES1 = previous_state[0]+EP
        SES2 = fracEP2
        ER2 = ((2*b*SES2 + b**2)**0.5) - b
        ER = fracEP1 + (fracEP2/EP)*ER2 #verif
        SEP = ER2
    else: #phase 2
        SES1 = previous_state[0]+EP
        SES2 = previous_state[1]+EP
        ER = ((2*b*SES2 + b**2)**0.5) - b - previous_state[2] #fonction cumul, donc retire valeur de la date precedente
        SEP = ((2*b*SES2 + b**2)**0.5) - b

    if Precip>0.:#comme si pluie vient en fin de journee
        SES1 = 0.
        SES2 = 0.
        SEP = 0.

    return ER, [SES1, SES2, SEP]#0., [SES1, SES2, 0.]#
    #a tester

def bEV(ACLIMc, ARGIs, HXs):
    """ calcul de  b: coeff empirique sans dimension pour l'evaporation du sol a partir de (STICS manual p 128) : 
        ACLIMc : wheather parameter qui depend de la vitesse moyenne du vent (1m-1-> 20; 2m.s-1-> 14)
        ARGIs : clay content de la couche superficielle (%)
        HXs: volumetric moisture content at fc of the surface layer """
    HA = ARGIs/1500.
    return 0.5*ACLIMc*(HXs-HA)*power(0.63-HA, 5./3.)

#bEV(ACLIMc=20, ARGIs=20, HXs=0.4)

#def evap_frac_zone(dz = [0., 0.3], ZESX=0.3):
#    """ distribution en profondeur du prelevement d'eau evapore sous un voxel de surface
#        integrale de y = -ZESX*ZESX*Z/2 + ZESX -> relation lineaire de la figure 7.4b STICS manual (p130)
#        correspond a equation 7.6 pour CFES=1 (lineaire)"""
#    if dz[0]>=0. and dz[0]<=ZESX: #1ere borne OK
#        if dz[1]>0. and dz[1]<=ZESX: #2e borne OK
#            frac_interval = (-(1./ZESX**2)*(dz[1]**2)+(2./ZESX)*dz[1]) - (-(1./ZESX**2)*(dz[0]**2)+(2./ZESX)*dz[0])
#        elif dz[1]>0. and dz[1]>ZESX: #2e superieure-> remplace par profondeur max ZESX
#            frac_interval = (-(1./ZESX**2)*(ZESX**2)+(2./ZESX)*ZESX) - (-(1./ZESX**2)*(dz[0]**2)+(2./ZESX)*dz[0])
#        else:
#            frac_interval = 0.
#    else:
#        frac_interval = 0.
#
#    return frac_interval

#evap_frac_zone(dz = [0., 0.2], ZESX=0.3)


#def check_soil_evap(evapo_tot, m_frac_evap, tsw_t, m_QH20min, dxyz, ZESX = 0.3):
#    """ verifie que evaporation du sol n'entraine pas de tsw sous la tsw mini d'un sol sec (air) ; sinon ajuste distribution en profondeur(report sur la couche inferieure), puis evapo_tot"""
#    m_evap = deepcopy(m_frac_evap)
#
#    #recherche des id couches sur lesquelles evaporation joue
#    nz, cumz = 0, []
#    for i in range(len(dxyz[2])):
#        cumz.append(dxyz[2][i])
#        nz = nz+1
#        if sum(cumz)>=ZESX:
#            break 
#
#    for x in range(len(dxyz[0])):
#        for y in range(len(dxyz[1])):
#            for z in range(nz-1):#report par colone de sol
#                #res = tsw_t[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#                res = tsw_t[z][x][y] - m_QH20min[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#                if res<=0.:
#                    #m_evap[z+1][x][y] = m_evap[z+1][x][y] + m_evap[z][x][y] - tsw_t[z][x][y]
#                    #m_evap[z][x][y] = tsw_t[z][x][y]
#                    m_evap[z+1][x][y] = m_evap[z+1][x][y]-m_QH20min[z+1][x][y] + m_evap[z][x][y]-m_QH20min[z][x][y] - tsw_t[z][x][y]#??
#                    m_evap[z][x][y] = tsw_t[z][x][y]- m_QH20min[z][x][y]#??
#
#            #derniere couche affectee
#            z = nz-1
#            res = tsw_t[z][x][y] - m_QH20min[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#            if res<=0.:
#                m_evap[z][x][y] = tsw_t[z][x][y] - m_QH20min[z][x][y]
#
#    return sum3(m_evap), m_evap

#test = deepcopy(S.tsw_t)
#test[0][0][0]=0.
#test[1][0][0]=0.
#test[2][0][0]=0.
#test[3][0][0]=0.
#test[4][0][0]=0.#002
#test[5][0][0]=0.#002
#test[6][0][0]=0.#002
#test[7][0][0]=0.002
#evapo_tot,  m_evap = check_soil_evap(sum3(m_frac_evap), m_frac_evap, test, S.dxyz, ZESX = 0.3)


#def distrib_evapNC(m_soil_vox, m_soil_vol, map_Evap0, ZESX=0.3):
#    """ map_Evap0 = liste des Evap des voxel de la couche superieure (matrice x, y) """
#    m = deepcopy(m_soil_vol)
#    if len(m)>1:
#        m.fill(0.)
#        for z in range(len(m)-1):
#            dz = [abs(m_soil_vox[z][0][0][2]), abs(m_soil_vox[z+1][0][0][2])]
#            frac = evap_frac_zone(dz , ZESX)
#            for x in range(len(m[z])):
#                for y in range(len(m[z][x])):
#                    m[z][x][y] = frac*map_Evap0[x][y]
#    else:
#        m.fill(1.)
#
#    return m

#map_Evap0 = array([[1.5]]) ## faire une fonction qui le genere initialement
#m_frac_evap = distrib_evapNC(m_soil_vox, m_soil_vol, map_Evap0, ZESX=0.3)






##gerer infiltrations 
def ls_1storder_vox(dxyz, x,y,z, opt=1):
    """ definit liste de voxels concomittents de 1er ordre / considere sol thorique, opt=1 renvoie uniquement voxel pour ecoulement vertical direct"""
    lx,ly = len(dxyz[0]), len(dxyz[1])
    if lx<3 or ly<3 or opt==1: #pas au moins 9 voxel ou ecoulement direct (opt=1)
        return [[x,y,z]]
    else: 
        if x==0: #bord gauche
            if y>0 and y<ly-1:#pas coin
                return [[lx-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
            elif y==0:#coin bas
                return [[lx-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, ly-1, z], [x, ly-1, z], [x+1, ly-1, z]]
            elif y==ly-1:#coin haut
                return [[lx-1, 0, z], [x, 0, z], [x+1, 0, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
        elif x == lx-1:#bord droit
            if y>0 and y<ly-1:#pas coin
                return [[x-1, y+1, z], [x, y+1, z], [0, y+1, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, y-1, z], [x, y-1, z], [0, y-1, z]]
            elif y==0:#coin bas
                return [[x-1, y+1, z], [x, y+1, z], [0, y+1, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, ly-1, z], [x, ly-1, z], [0, ly-1, z]]
            elif y==ly-1:#coin haut
                return [[x-1, 0, z], [x, 0, z], [0, 0, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, y-1, z], [x, y-1, z], [0, y-1, z]]
        elif y==0:#cote bas sans les coins
            return [[x-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, ly-1, z], [x, ly-1, z], [x+1, ly-1, z]]
        elif y==ly-1:#cote haut sans les coins
            return [[x-1, 0, z], [x, 0, z], [x+1, 0, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
        else: #au centre
            return [[x-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]

#ls_1storder_vox(S.dxyz, 0,0,0, opt='1')

def infil_layer(tsw_t, m_QH20max, in_ , dxyz, idz, opt=1):#in_=map_PI, 
    new = tsw_t[idz]+in_
    out_ = deepcopy(in_)
    out_.fill(0.)
    for x in range(len(tsw_t[idz])):
        for y in range(len(tsw_t[idz][x])):
            if new[x][y] > m_QH20max[idz][x][y]:# si superieur a fc infiltration
                q_out = new[x][y]-m_QH20max[idz][x][y]
                ls_v = ls_1storder_vox(dxyz, x,y,idz, opt)#distribution entre les 1st order ; mettre opt=1 si veut forcer verticalement / 2 si
                if len(ls_v)>1:
                    ponder = [0.0416666, 0.0416666, 0.0416666, 0.0416666, 2/3., 0.0416666, 0.0416666, 0.0416666, 0.0416666]# 2/3 en dessous 1/3 au premier ordre
                else:
                    ponder = [1.]

                for i in range(len(ls_v)):#distribution dans les voxels du dessous
                    nx, ny = ls_v[i][0], ls_v[i][1]
                    out_[nx][ny] = out_[nx][ny]+ponder[i]*q_out #bug: out_[nx][ny] = ponder[i]*q_out Resolu! 
                 
                new[x][y] = m_QH20max[idz][x][y]# ce qui reste dans le voxel = fc

    return new, out_

#infil_layer(asw_t, m_QH20max, array([[5.]]), dxyz, 0)

def distrib_PI(tsw_t, m_QH20max, map_PI, dxyz, opt=1):
    """ distribution des precipitations+irrigation dans le sol 
    ecoulement vertical uniquement 
    tout ce qui est au dessus de la capacite au champ est transfere en 1 jour - drainage = sortie de la couche la plus profonde """
    in_ = map_PI
    ls_out = []
    for z in range(len(tsw_t)):
        new, out_ = infil_layer(tsw_t, m_QH20max, in_, dxyz, z, opt)
        ls_out.append(out_)
        in_ = out_
        tsw_t[z] = new

    return tsw_t, ls_out #tsw_t, out_
    # prevoir ditribution normale en fonction de distance au centre du voxel pour "homogeneiser" en profondeur et eviter chemin preferentiel vers drainage et limiter sensibilite a la taille des voxels
    # peut aussi ajouter stockage jusqua saturation puis relargage selon vitesse d'infiltration (STICS p224)

#map_PI = array([[5.]])
#distrib_PI(asw_t, m_QH20max, map_PI)

#Rq : distrib pas verticale buggee (opt=2) / coin extene (xn, yn, recoit plus d'infiltration) ;  ok??? Resolu dans infil_layer!


##  puis faire bilan
#def TSW_NC(previousTSW, map_PI, m_frac_evap, m_frac_transpi, m_QH20max, dxyz, opt=1):
#    """ """
#    # 1) applique pluie et calcule D0 (drainage sans transpi pour saturer)
#    m, D0 = distrib_PI(previousTSW, m_QH20max, map_PI, dxyz,opt)
#    # 2) retire evapo et transpi
#    m = m - m_frac_evap - m_frac_transpi
#    # 3) remet ce qui a drainer et recalcule drainage residuel apres transpi (=transpire en priorite ce qui est tombe pdt la journee)
#    m, D1 = distrib_PI(m, m_QH20max, D0, dxyz)
#    return m, D1



## gestion des racines
##rq: density above 0.5 cm.cm-3 are not taken into account for absorption of water and N in STICS (p90)
##rq2: valeur de 0.1 cm.cm-3 pour definir profondeur d'enracinement dans L2SU(p138)

def vert_roots(dxyz, lvopt):
    """ pour generer un profil vertical de densite racinaire a partir d'une liste - (homogene en x,y)"""
    m_roots = []
    for z in range(len(dxyz[2])):
        v = []
        for x in range(len(dxyz[0])):
            vv = []
            for y in range(len(dxyz[1])):
                vv.append(lvopt[z])

            v.append(vv)
        m_roots.append(v)

    return array(m_roots)

#R1 = vert_roots(dxyz, [0.5]*nstrate[2])
#R2 = vert_roots(dxyz, [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05])

#def effective_root_length(m_root, tresh = 0.5):#faire avec une liste de root systems ls_root_syst pour competition
#    """ for a single root system """
#    m = deepcopy(m_root)
#    for z in range (len(m_root)):
#        for x in range (len(m_root[z])):
#            for y in range (len(m_root[z][x])):
#                if m_root[z][x][y]>tresh:
#                    m[z][x][y] = tresh
#                else:
#                    m[z][x][y] = m_root[z][x][y]
#    return m

#effective_root_length(R2, 0.5)


def effective_root_lengths(ls_roots, tresh=0.5):
    """ for multiple root systems in competition
    H0 : perfect symetry locally (in a voxel) for ressource aquisition
    treshold of effective density similar to single species
    fraction of effective density to each species/plant = proportion of total root length density"""
    m = deepcopy(ls_roots)
    nb_r = len(ls_roots)
    tot_root_dens = sum_mat(ls_roots)  # sum bug
    for rt in range(nb_r):
        # frac_root = divide(ls_roots[rt], tot_root_dens)## rq: gere bien les voxels vides
        m[rt] = where(tot_root_dens > tresh, tresh * ls_roots[rt] / tot_root_dens, ls_roots[rt])

    return m


#RLprofil = {0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}#sortie l-system
def build_ls_roots(RLprofil, S):
    """ version 1 root system :a modifier pour prendre en compte une liste de RLprofil"""
    idz = list(RLprofil.keys())
    idz.sort()
    RLprofil_ls = []
    for i in idz: RLprofil_ls.append(RLprofil[i])

    R1 = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
    ls_roots = [R1]#[R1, R2]#[R3]#
    return ls_roots


#RLprofil = [{0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}, {0: 0.145, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}...]#sortie l-system

def build_ls_roots_mult(RLprofil, S):
    """ version 2: prends en compte une liste de RLprofil"""
    nbp = len(RLprofil)
    idz = list(RLprofil[0].keys())
    idz.sort()
    ls_roots = []
    for p in range(nbp):
        RLprofil_ls = []
        for i in idz: RLprofil_ls.append(RLprofil[p][i])
        R = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
        ls_roots.append(R)# = [R1]#[R1, R2]#[R3]#

    return ls_roots



def RLprof_t(t, ncouches_sol):
    """ pour generer un profil de densite qui varie au cours du temps / phase de test """
    RL_profil_max = [0.64]*31+[0.64, 0.6, 0.5, 0.32, 0.16, 0.08, 0.04, 0.02, 0.01]
    n0 = ncouches_sol-t
    if t<ncouches_sol:
        res  = RL_profil_max[-1-t:-1]+[0]*n0
    else:
        res = RL_profil_max
    return res

#for i in range(50):print RLprof_t(i, ncouches_sol)
    


################## TESTS ###############################

def test_uni1():
    ###############
    ## DEMO 1
    ## sol 1D (S)
    ## 1 seul systeme racinaire en dynamique
    ###############

    ## meteo
    Et0 = 4.
    precips = [0.]*29+[90.]# precipitations = 20mm tous les 30 jours
    precips = precips*20
    Precip = precips[0]#0.
    Irrig = 0.

    ## sol
    par_sol = {'13':{'KST': '5', 'teta_ad': 0.0078506705623363031, 'soil number': '13', 'teta_fc': 0.35816688969184857, 'WCST': '0.503', 'teta_wp': 0.11250451238171076, 'soil type': '13 - loam', 'teta_sat': 0.503, 'gamma_theo': '0.08489778'}}

    dz_sol = 5. #cm
    ncouches_sol = 40
    soil_type = 13 #"13 - loam"
    ZESX = 0.3#m
    CFES = 1.
    Uval = 3.#U quantite d'eau dans une couche superieure en mm (5 par default)
    b = 0.63

    ##vegetation
    espi_t = [0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006,0.995315225515,0.312710721209,0.430217175269,0.527633447259,0.608394373323,0.675347532642,0.745670599554,0.827513705872,0.884631623932,0.923898514961,0.950492197719,0.968236501678,0.979901799595,0.987458191053,0.992281421006]
    ls_epsi = [espi_t[0]]

    ## soil initialisation
    S = Soil(par_sol, soil_number = [soil_type]*ncouches_sol, dxyz = [[1.], [1.], [dz_sol/100.]*ncouches_sol], vDA=[1.25]*ncouches_sol, ZESX=ZESX, CFES=CFES)

    #RL_profil = table[0] #profil exporte de l-py via table
    RL_profil = RLprof_t(10, ncouches_sol)
    R1 = vert_roots(S.dxyz, RL_profil)#*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
    ls_roots = [R1]#[R1, R2]#[R3]#

    ##bouble journaliere
    ls_transp, evapo_tot, D, state,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig,previous_state=[0.,0.,0.], ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    for i in range(1,len(espi_t)):
        previous_asw = sum3(S.asw_t)
        Precip = 0.#precips[i]
        Irrig = 0.
        ls_epsi = [espi_t[i]]

        RL_profil = RLprof_t(10+i, ncouches_sol)#table[i]
        R1 = vert_roots(S.dxyz, RL_profil)#*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
        ls_roots = [R1]#[R1, R2]#[R3]#

        ls_transp, evapo_tot, D, state,  m_frac_transpi, m_frac_evap, ls_ftsw =  S.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, state, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=0.63, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)

        print(i, ls_ftsw, sum3(S.asw_t),  sum3(m_frac_evap), sum3(m_frac_transpi), state)

    #    ##visu
    #    #bx = Box(Vector3(1.,1.,1.))
    #    Monviewer=Viewer
    #    MaScene=Scene()
    #    MaScene = S.plot_soil_properties (vals=S.ftsw_t, MaScene=Scene(), col_scale=5)#plot_soilWC (S.dxyz, S.m_soil_vox, S.ftsw_t, MaScene)
    #    Monviewer.display(MaScene)

#test_uni1()


def test_uni2():
    #################
    ## DEMO 2
    ## avec sol 3D (S2)
    ## 2 systemes racinaires fixes en competition (homogenes horizontalement) R3/R4
    ## en statique
    #################


    ##meteo
    Et0 = 2.
    #precips = [0.]*29+[90.]# precipitations = 20mm tous les 30 jours
    #precips = precips*20
    precips = [0.]*500#pas de precipitations (dessechement)
    Precip = precips[0]#0.
    Irrig=0.

    ## sol
    par_sol = {'13':{'KST': '5', 'teta_ad': 0.0078506705623363031, 'soil number': '13', 'teta_fc': 0.35816688969184857, 'WCST': '0.503', 'teta_wp': 0.11250451238171076, 'soil type': '13 - loam', 'teta_sat': 0.503, 'gamma_theo': '0.08489778'}}

    dz_sol = 5. #cm
    ncouches_sol = 10
    soil_type = 13 #"13 - loam"
    ZESX = 0.3 #m
    CFES=1.
    b=0.63
    Uval = 3.#sum(S.asw_t[0])#U quantite d'eau dans une couche superieure en mm (5 par default)
    stateEV = [0.,0.,0.]

    ### soil initialisation
    S2 = Soil(par_sol, soil_number = [soil_type]*ncouches_sol, dxyz = [[0.2]*5, [0.2]*5, [0.2]*10], vDA=[1.25]*ncouches_sol, ZESX=ZESX, CFES=CFES)

    ## vegetation
    ls_epsi = [0.4, 0.4]##[0.4]# !! correspondance avec les nb de root systems!
    ## 3D root distributions (fixe)
    R3 = zeros(shape(S2.m_1))
    R3[0:5, 1:3, 1:3] = ones([5,2,2])*0.25
    R4 = zeros(shape(S2.m_1))
    R4[0:7, 0:4, 2:4] = ones([7,4,2])*0.15

    ls_roots = [R3, R4]

    ##bouble journaliere
    ls_transp, evapo_tot, D, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S2.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, stateEV, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S2.update_ftsw()

    for i in range(40):
        Precip = precips[i]
        ls_transp, evapo_tot, D, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S2.stepWBmc(Et0, ls_roots, ls_epsi, Precip, Irrig, stateEV, ZESX=ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
        S2.update_ftsw()

    #    ##visu
    #    Monviewer=Viewer
    #    MaScene=Scene()
    #    MaScene = S2.plot_soil_properties (vals=S2.ftsw_t, MaScene=Scene(), col_scale=5)
    #    Monviewer.display(MaScene)

    print((S2.ftsw_t))

#test_uni2()
