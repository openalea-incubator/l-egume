#visus sol
library("readxl")

dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\output"

setwd(dir)#



ls_files <- list.files(dir)#(dir0)#


#recupere la liste des outHR file names du dossier de travail
ls_outHR <- ls_files[grepl('outHR', ls_files)]


filHR <- ls_outHR[3]#ls_outHR[2]
#"outHR_52_l-egume_solnu-solnu_homogeneous0_scenario-1--1_LusignanAsso16_0_LusignanAsso16_.csv"

dat <- read.table(filHR, sep=';', header=T)
sdat <- split(dat, dat$var)



#plot spatio-temporel ftsw
DOYs <- sdat[['FTSW']]$DOY
vals <- sdat[['FTSW']][,c(-1,-2)]
lscol <- colorRampPalette(c("blue", "red"))( 101 ) #palette de couleur
nbcouches <- dim(vals)[2]
  

plot(-10,-10,ylim=c(-5*nbcouches,-0), xlim=c(min(DOYs),max(DOYs)), main='FTSW', xlab='DOY', ylab='soil depth')

#draw a sequence of recatngle
#jour 1
j <- 1
for (j in 1:dim(vals)[1])
{
  ftswj <- as.numeric(vals[j,])
  cols <- rev(lscol[round((1-ftswj)*100,0)+1])
  
  xleft <- rep(DOYs[j], nbcouches)
  ybottom <- seq(-5*nbcouches, -5, 5)
  xright <- rep(DOYs[j+1], nbcouches)
  ytop <- seq(-5*nbcouches, -5, 5)+5
  rect(xleft, ybottom, xright, ytop, col=cols, border=cols)
}



#plot spatio-temporel NO3
DOYs <- sdat[['m_NO3']]$DOY
valsNO3 <- sdat[['m_NO3']][,c(-1,-2)]
lscol <- colorRampPalette(c("brown", "green"))( 101 ) #palette de couleur
maxNO3 <- max(valsNO3)
vals <- valsNO3/maxNO3
nbcouches <- dim(vals)[2]


plot(-10,-10,ylim=c(-5*nbcouches,-0), xlim=c(min(DOYs),max(DOYs)), main='NO3', xlab='DOY', ylab='soil depth')

#draw a sequence of recatngle
#jour 1
j <- 1
for (j in 1:dim(vals)[1])
{
  ftswj <- as.numeric(vals[j,])
  cols <- rev(lscol[round((1-ftswj)*100,0)+1])
  
  xleft <- rep(DOYs[j], nbcouches)
  ybottom <- seq(-5*nbcouches, -5, 5)
  xright <- rep(DOYs[j+1], nbcouches)
  ytop <- seq(-5*nbcouches, -5, 5)+5
  rect(xleft, ybottom, xright, ytop, col=cols, border=cols)
}




#faire les moyennes humidite pour horizon de mesure

idH1 <- 1:6#0-30
idH2 <- 7:12#30-60
idH3 <- 13:18#60-90

DOYs <- sdat[['HRp']]$DOY
valsHR <- sdat[['HRp']][,c(-1,-2)]

HRpH1 <- as.numeric(rowMeans(valsHR[,idH1]))
HRpH2 <- as.numeric(rowMeans(valsHR[,idH2]))
HRpH3 <- as.numeric(rowMeans(valsHR[,idH3]))



#ouverture du fichier obs
library("readxl")
f_ <- file.choose()
#"C:\\devel\\l-egume\\legume\\test\\obs\\Asso-sol-nu_obs.xls"
onglet <- "Asso-sol-nu"
obs <- read_excel(f_, sheet = onglet, col_names = TRUE, na = "")

layout(matrix(1:3,3,1))
plot(DOYs, HRpH1, type='l', ylim=c(0,26), main='HRp H1', col="red")
points(obs$DOY, obs$`HR(1)`)
plot(DOYs, HRpH2, type='l', ylim=c(0,26), main='HRp H2', col="red")
points(obs$DOY, obs$`HR(2)`)
plot(DOYs, HRpH3, type='l', ylim=c(0,26), main='HRp H3', col="red")
points(obs$DOY, obs$`HR(3)`)

#NO3
valsNO3 <- sdat[['m_NO3']][,c(-1,-2)]
surfsolref <- 0.25*0.25

NO3H1 <- as.numeric(rowSums(valsNO3[,idH1]))/surfsolref*10000
NO3H2 <- as.numeric(rowSums(valsNO3[,idH2]))/surfsolref*10000
NO3H3 <- as.numeric(rowSums(valsNO3[,idH3]))/surfsolref*10000


layout(matrix(1:3,3,1))
plot(DOYs, NO3H1, type='l', ylim=c(0,20), main='NO3 H1', col="red")
points(obs$DOY, obs$`AZnit(1)`)
plot(DOYs, NO3H2, type='l', ylim=c(0,20), main='HRp H2', col="red")
points(obs$DOY, obs$`AZnit(2)`)
plot(DOYs, NO3H3, type='l', ylim=c(0,20), main='HRp H3', col="red")
points(obs$DOY, obs$`AZnit(3)`)




#faire les moyennes pour plusieurs fichier




