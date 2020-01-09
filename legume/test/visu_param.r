## pour visualiser / editer les fichiers parametres


library(readxl)
namef<-file.choose()
#"C:\\devel\\l-egume\\legume\\input\\Parametres_plante_exemple.xls"

param_name <- "Fix2"#"timbale"#'giga'#'leo'#'sevanskij'#'formica'#'canto'#

dat <- read_excel(namef,sheet=param_name, col_names = T, na="")




# profils de dimension des organes

ranks <- 1:40

##Longueur feuille
lmax <- as.numeric(dat[dat$name == 'Lfeuille', 4])
x <- dat[dat$name == 'profilLeafI_Rlen', ]
up <- lmax*(ranks*as.numeric(x[1,4]) + as.numeric(x[2,4]))
down <- lmax*(ranks*as.numeric(x[3,4]) + as.numeric(x[4,4]))
profil <- apply(data.frame(up, down), 1, min) #obttention des mini
plot(ranks, profil, main=paste(param_name, '- profils'), type='l', ylab='organ length (cm)', ylim=c(0,10))

##Largeur des feuilles

x <- as.numeric(dat[dat$name == 'profilLeafI_Rlarg', ])
up <- (ranks*as.numeric(x[1,4]) + as.numeric(x[2,4]))
down <- (ranks*as.numeric(x[3,4]) + as.numeric(x[4,4]))
data.frame(up, down)
profil_larg <- apply(data.frame(up, down), 1, max)*profil #obttention des mini
points(ranks, profil_larg, lty=2, type='l')
##?? verif? le *profil fai un truc en r2?


##Longueur petiole
lmax <- as.numeric(dat[dat$name == 'Lpet', 4])
x <- dat[dat$name == 'profilPetI', ]
up <- lmax*(ranks*as.numeric(x[1,4]) + as.numeric(x[2,4]))
down <- lmax*(ranks*as.numeric(x[3,4]) + as.numeric(x[4,4]))
profil <- apply(data.frame(up, down), 1, min) #obttention des mini
points(ranks, profil, type='l', col=2)

##Longueur en
lmax <- as.numeric(dat[dat$name == 'Len', 4])
x <- dat[dat$name == 'profilLeafI_Rlen', ]
up <- lmax*(ranks*as.numeric(x[1,4]) + as.numeric(x[2,4]))
down <- lmax*(ranks*as.numeric(x[3,4]) + as.numeric(x[4,4]))
profil <- apply(data.frame(up, down), 1, min) #obttention des mini
points(ranks, profil, type='l', col=3)




# developpement des tiges (I, II, delai, tallage)
TT <- 1:1000
vI <- as.numeric(1/dat[dat$name == 'phyllochron', 4])
vII <- as.numeric(1/dat[dat$name == 'phyllochronII', 4])
delai_deb <- as.numeric(dat[dat$name == 'delai_deb', 4])/(1/vI)
vtall <- as.numeric(dat[dat$name == 'RvitTallage', 4])*vI
nmax <- as.numeric(dat[dat$name == 'nshoots', 4])
deltal <- as.numeric(dat[dat$name == 'debTallage', 4])
up <- TT*vtall-deltal
down <- TT*0+nmax
nbsh <- apply(data.frame(up, down), 1, min) #obttention des mini

plot(TT, vI*TT, main=paste(param_name, '- development potentiel'), type='l', ylab='Nb', ylim=c(0,40))
points(TT,vI*TT-delai_deb, type='l', lty=2)
points(TT,vII*TT, col=2, type='l')
points(TT,nbsh,col=4,type='l')



# coordiantion de croissance des organes


# allocation (racine / pivot / racines fines)


# developpement des racines (mere/fille - vitesse/diametre)
Diams <- seq(0.,0.15,0.01)
Dmin <- as.numeric(dat[dat$name == 'Dmin', 4])
Dmax <- as.numeric(dat[dat$name == 'Dmax', 4])
DIDm <- as.numeric(dat[dat$name == 'DIDm', 4])
ELmax <- as.numeric(dat[dat$name == 'ELmax', 4])
ELD <- ELmax/(Dmax-Dmin) #dat[dat$name == 'ELD', 4]

Dfille <- Diams*DIDm + (-Dmin*(DIDm-1))
Vfille <- Diams*ELD+ (ELmax-ELD*Dmax)

layout(matrix(1:2,1,2))
plot(Diams, Diams, type='l', ylab='D Fille (cm)', xlab='D mere (cm)', main=paste(param_name, '- Root Mere:Fille'))
points(Dmin,Dmin,pch=16,col=2)
points(Diams, Dfille, type='l', lty=2)
lines(Diams*0+Dmax, Diams, col=2)

plot(Diams, Vfille, type='l', ylab='vitesse allongement (cm.°Cj-1)', xlab='Diametre (cm)', ylim=c(0,0.12), main=paste(param_name, '- vitesse~D root'))
points(Dmax,ELmax,col=2, pch=16)
points(Dmin,0,col=2, pch=16)



##stress hydrique
sigmo <- function(v,del,x)
{
  #sigmo croissante norme entre 0 et 1
  valeursig <- 1 - 1/(1+exp(v*(x-del)))
  return (valeursig)
}

FTSWrange <- seq(0,1,0.05)

#croissance
slope <- as.numeric(dat[dat$name == 'WaterTreshExpSurf', 4][1,])
ftsw50 <- as.numeric(dat[dat$name == 'WaterTreshExpSurf', 4][2,])
plot(FTSWrange, sigmo(slope, ftsw50, FTSWrange), type='l', xlab='FTSW', ylab='relative rate') #noir
#devI
slope <- as.numeric(dat[dat$name == 'WaterTreshDevI', 4][1,])
ftsw50 <- as.numeric(dat[dat$name == 'WaterTreshDevI', 4][2,])
points(FTSWrange, sigmo(slope, ftsw50, FTSWrange), type='l', col=3)#vert
#devII
slope <- as.numeric(dat[dat$name == 'WaterTreshDevII', 4][1,])
ftsw50 <- as.numeric(dat[dat$name == 'WaterTreshDevII', 4][2,])
points(FTSWrange, sigmo(slope, ftsw50, FTSWrange), type='l', col=2)#rouge
#RUE
slope <- as.numeric(dat[dat$name == 'WaterTreshRUE', 4][1,])
ftsw50 <- as.numeric(dat[dat$name == 'WaterTreshRUE', 4][2,])
points(FTSWrange, sigmo(slope, ftsw50, FTSWrange), type='l', col=4)#bleu
#fixation
slope <- as.numeric(dat[dat$name == 'WaterTreshFix', 4][1,])
ftsw50 <- as.numeric(dat[dat$name == 'WaterTreshFix', 4][2,])
points(FTSWrange, sigmo(slope, ftsw50, FTSWrange), type='l', col=5)#cyan


# azote
NNIrange <- seq(0,1,0.05)

#croissance
slope <- as.numeric(dat[dat$name == 'NTreshExpSurf', 4][1,])
nni50 <- as.numeric(dat[dat$name == 'NTreshExpSurf', 4][2,])
plot(NNIrange, sigmo(slope, nni50, NNIrange), type='l', xlab='NNI', ylab='relative rate') #noir
#dev
slope <- as.numeric(dat[dat$name == 'NTreshDev', 4][1,])
nni50 <- as.numeric(dat[dat$name == 'NTreshDev', 4][2,])
points(FTSWrange, sigmo(slope, nni50, NNIrange), type='l', col=3)#vert
#RUE
slope <- as.numeric(dat[dat$name == 'NTreshRUE', 4][1,])
nni50 <- as.numeric(dat[dat$name == 'NTreshRUE', 4][2,])
points(FTSWrange, sigmo(slope, nni50, NNIrange), type='l', col=4)#bleu




#reponses au rayonnement




#reponse a la temperature



date$name

