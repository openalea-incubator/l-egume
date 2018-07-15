dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\multisim\\sorties"
# "H:\\simul\\sorties" 
#"C:\\devel\\grassland\\sorties\\mixture Fix-NonFix scenario1"
setwd(dir)
ls_files <- list.files(dir)


#recupere la liste des toto file names
ls_toto <- ls_files[grepl('toto', ls_files)]

dtoto <- as.data.frame(t(as.data.frame(strsplit(ls_toto, '_'))))[,c(2,3,5,6,7,8,9)]
row.names(dtoto) <- 1: length(dtoto[,1])
names(dtoto) <- c('usm','lsystem','mix','damier','scenario','Mng', 'seed')
dtoto$name <- ls_toto
dtoto$seed <- substr(as.character(dtoto$seed), 1, 1)


#recuperer les fichiers toto dans une liste ltoto
ltoto <- vector('list', length(ls_toto))
names(ltoto) <- ls_toto

for (i in 1:length(ls_toto))
{
  name <- ls_toto[i]
  ltoto[[name]] <- read.table(name, header=T, sep=';')
}






#visu des rendement moyen m2 / a un DOY 
surfsolref <- NULL
nbplt <- NULL
nbplt1 <- NULL
nbplt2 <- NULL

#DOYScoupe <- c(165,199,231,271,334)#Avignon
DOYScoupe <- c(187,229,282,334)#Lusignan
DOYdeb <- 60
idDOYScoupe <- DOYScoupe - DOYdeb
Ytot <- NULL
Ycoupe <- NULL

YEsp1 <- NULL
YEsp2 <- NULL

QNfix <- NULL
QNupttot <- NULL
QNuptleg <- NULL

for (i in 1:length(ls_toto))
{
  name <- ls_toto[i]
  dat <- ltoto[[name]]
  s <- dat[dat$V1=='pattern',3]#m2
  surfsolref <- cbind(surfsolref, as.numeric(as.character(s)))
  nb <- length(dat)-2
  nbplt <- cbind(nbplt, nb)
  
  #Y Totaux
  MSaerien <- as.matrix(dat[dat$V1=='MSaerien' & dat$steps %in% DOYScoupe,3:(3+nb-1)], ncol=nb)
  ProdIaer <- rowSums(MSaerien) / s
  Ycoupe <- rbind(Ycoupe, ProdIaer)
  Ytot <- cbind(Ytot, sum(ProdIaer))#cumul des 5 coupes

  
  #N totaux et fixation
  Qfix <- as.matrix(dat[dat$V1=='Qfix',3:(3+nb-1)], ncol=nb)
  Qfix <- as.numeric(rowSums(Qfix) / s)
  Nuptake_sol_tot <- as.matrix(dat[dat$V1=='Nuptake_sol',3:(3+nb-1)], ncol=nb)
  Nuptake_sol_tot <- as.numeric(rowSums(Nuptake_sol_tot) / s)
  QNfix <-cbind(QNfix, sum(Qfix))
  QNupttot <- cbind(QNupttot, sum(Nuptake_sol_tot))

  #YEsp1
  esp <- 'Fix2'#'Fix3'#'Fix1'#'Fix' #pourquoi c'est ce nom au lieu de Fix???
  esp2 <- 'nonFixSimTest'#'nonFix1'#'nonFix0' #pourquoi c'est ce nom au lieu de Fix???
  
  nomcol <- names(ltoto[[name]])
  idcols <- grepl(esp, nomcol) & !grepl(esp2, nomcol)#contient esp1 et pas esp2
  dat1 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
  nb1 <- length(dat1)-2
  nbplt1 <- cbind(nbplt1, nb1)
  if (nb1>0)
  {
    MS1 <- as.matrix(dat1[dat1$V1=='MSaerien' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)
    ProdIaer1 <- rowSums(MS1) / s
    Nuptake_sol_leg <- as.matrix(dat1[dat1$V1=='Nuptake_sol',3:(3+nb1-1)], ncol=nb1)
    Nuptake_sol_leg <- as.numeric(rowSums(Nuptake_sol_leg) / s)
  }
  else
  {
    ProdIaer1 <- 0 #pas de plante de l'esp1
    Nuptake_sol_leg <- 0
  }
  YEsp1 <- cbind(YEsp1, sum(ProdIaer1))#cumul des 5 coupes
  QNuptleg <- cbind(QNuptleg, sum(Nuptake_sol_leg))
  
  #YEsp2
  idcols <- grepl(esp2, nomcol)
  dat2 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
  nb2 <- length(dat2)-2
  nbplt2 <- cbind(nbplt2, nb2)
  if (nb2>0)
  {
    MS2 <- as.matrix(dat2[dat2$V1=='MSaerien' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)
    ProdIaer2 <- rowSums(MS2) / s
  }
  else
  {
    ProdIaer2 <- 0 #pas de plante de l'esp2
  }
  YEsp2 <- cbind(YEsp2, sum(ProdIaer2))#cumul des 5 coupes
  
}

dtoto$surfsolref <- as.numeric(surfsolref)
dtoto$nbplt <- as.numeric(nbplt)
dtoto$nbplt1 <- as.numeric(nbplt1)
dtoto$nbplt2 <- as.numeric(nbplt2)
dtoto$Ytot <- as.numeric(Ytot)
dtoto$densite <- dtoto$nbplt/dtoto$surfsolref
dtoto$densite1 <- dtoto$nbplt1/dtoto$surfsolref
dtoto$YEsp1 <- as.numeric(YEsp1)
dtoto$densite2 <- dtoto$nbplt2/dtoto$surfsolref
dtoto$YEsp2 <- as.numeric(YEsp2)
dtoto$Semprop1 <- dtoto$densite1/dtoto$densite
dtoto$Yprop1 <- dtoto$YEsp1 / (dtoto$YEsp1 +dtoto$YEsp2)
dtoto$Yprop2 <- dtoto$YEsp2 / (dtoto$YEsp1 +dtoto$YEsp2)
dtoto$QNfix <- as.numeric(QNfix)
dtoto$QNupttot <- as.numeric(QNupttot)
dtoto$QNuptleg <- as.numeric(QNuptleg)
dtoto$QNtot <- dtoto$QNfix + dtoto$QNupttot

#?? tjrs meme denite de esp1??
#pb de nom de esp1 qui est dans celui de esp2!! -> bug du grepl

#hist(Yprop1)
#mean(Yprop1)




#calcul des moy et sd+graphs
# + ajustement carres

#dataframe des valeurs moyennes
x <- by(dtoto$YEsp1, as.factor(dtoto$densite1), mean)
tabmoy <- data.frame(densite1=as.numeric(names(x)), YEsp1=as.numeric(x))
x <- by(dtoto$YEsp2, as.factor(dtoto$densite1), mean)
tabmoy$YEsp2 <- as.numeric(x)
x <- by(dtoto$Ytot, as.factor(dtoto$densite1), mean)
tabmoy$Ytot <- as.numeric(x)
x <- by(dtoto$YEsp1, as.factor(dtoto$densite1), sd)
tabmoy$YEsp1sd <- as.numeric(x)
x <- by(dtoto$YEsp2, as.factor(dtoto$densite1), sd)
tabmoy$YEsp2sd <- as.numeric(x)
x <- by(dtoto$Ytot, as.factor(dtoto$densite1), sd)
tabmoy$Ytotsd <- as.numeric(x)
x <- by(dtoto$Semprop1, as.factor(dtoto$densite1), mean)
tabmoy$Semprop1 <- as.numeric(x)

x <- by(dtoto$QNtot, as.factor(dtoto$densite1), mean)
tabmoy$QNtot <- as.numeric(x)
x <- by(dtoto$QNupttot, as.factor(dtoto$densite1), mean)
tabmoy$QNupttot  <- as.numeric(x)
x <- by(dtoto$QNuptleg, as.factor(dtoto$densite1), mean)
tabmoy$QNuptleg  <- as.numeric(x)
x <- by(dtoto$QNfix, as.factor(dtoto$densite1), mean)
tabmoy$QNfix  <- as.numeric(x)



#plot des valeur moyennes Ytot
xx <- tabmoy$Semprop1
plot(xx, tabmoy$Ytot, ylim=c(0,2200), xlab='Sowing proportion (Esp. 1)', ylab='Shoot biomass (g.m-2)')
#segments(tabmoy$Semprop1, tabmoy$Ytot, tabmoy$Semprop1, tabmoy$Ytot+tabmoy$Ytotsd)
#segments(tabmoy$Semprop1, tabmoy$Ytot, tabmoy$Semprop1, tabmoy$Ytot-tabmoy$Ytotsd)
segments(xx[1], tabmoy$Ytot[1], xx[7], tabmoy$Ytot[7], lty=2)
modeltot <- smooth.spline(xx, tabmoy$Ytot)
lines(modeltot)
inttot = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) - (tabmoy$Ytot[1]+tabmoy$Ytot[7])/2
text(0.1,1250, round(inttot,2))

points(xx, tabmoy$YEsp1,col=2)
segments(xx[1], tabmoy$YEsp1[1], xx[7], tabmoy$YEsp1[7], lty=2, col=2)
modelesp1 <- smooth.spline(xx, tabmoy$YEsp1)
lines(modelesp1, col=2)
intesp1 = sum(predict(modelesp1, seq(0,1,0.001))$y*0.001) - (tabmoy$YEsp1[1]+tabmoy$YEsp1[7])/2
text(0.1,1200, round(intesp1,2),col=2)

points(xx, tabmoy$YEsp2,col=4)
segments(xx[1], tabmoy$YEsp2[1], xx[7], tabmoy$YEsp2[7], lty=2, col=4)
modelesp2 <- smooth.spline(xx, tabmoy$YEsp2)
lines(modelesp2, col=4)
intesp2 = sum(predict(modelesp2, seq(0,1,0.001))$y*0.001) - (tabmoy$YEsp2[1]+tabmoy$YEsp2[7])/2
text(0.1,1150, round(intesp2,2),col=4)

#ajouter une ponderation par ecart a 50% pour eviter de selectionner vers esp pure?
#ajouter surface au dessus de transgressive over-yilding? (au dessus de ligne horizontale?)
#ajouter overyieding max? et prop correcpondante?




#plot des valeurs moyennes Ntot
xx <- tabmoy$Semprop1
plot(xx, tabmoy$QNtot, ylim=c(0,80), xlab='Sowing proportion (Esp. 1)', ylab='Plant N (g N.m-2)')
segments(xx[1], tabmoy$QNtot[1], xx[7], tabmoy$QNtot[7], lty=2)
modeltot <- smooth.spline(xx, tabmoy$QNtot)
lines(modeltot)
intoverN = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) - (tabmoy$QNtot[1]+tabmoy$QNtot[7])/2
intQNtot = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) 
text(0.1,45, paste(round(intoverN,2), '(over)'))

points(xx, tabmoy$QNupttot,col=2)
#segments(xx[1], tabmoy$QNupttot[1], xx[7], tabmoy$QNupttot[7], lty=2, col=2)
modelesp1 <- smooth.spline(xx, tabmoy$QNupttot)
lines(modelesp1, col=2)
intNupt = sum(predict(modelesp1, seq(0,1,0.001))$y*0.001) 
intFix = intQNtot-intNupt
text(0.1,42, paste(round(intFix,2), '(Fix)'),col=1)
text(0.1,40, paste(round(intNupt,2), '(Nupt)'),col=2)

points(xx, tabmoy$QNuptleg,col=4)
segments(xx[1], tabmoy$QNuptleg[1], xx[7], tabmoy$QNuptleg[7], lty=2, col=4)
modeleg <- smooth.spline(xx, tabmoy$QNuptleg)
lines(modeleg, col=4)
intleg = sum(predict(modeleg, seq(0,1,0.001))$y*0.001) - (tabmoy$QNuptleg[1]+tabmoy$QNuptleg[7])/2
text(0.1,37, paste(round(intleg,2), '(leg)'),col=4)

#valeurs semblent basses???




########################################################################"



#ajouter les fit x2 ou x3
parab2 <- function(a,yy1,yy0,x)
{
  #parabole forcee par y0 et y1
  a*x*x+(yy1-yy0-a)*x+yy0
}


parab3 <- function(a,b,yy1,yy0,x)
{
  #polynome ordre 3 force par y0 et y1
  a*x*x*x+b*x*x+(yy1-yy0-a-b)*x+yy0
}


#poly2??

x1 <- tabmoy$Semprop1
y1 <- tabmoy$Ytot
y_0 <- tabmoy$Ytot[1]
y_1 <- tabmoy$Ytot[7]


startlist=list(a=1., yy0=y_0, yy1=y_1)
minis <- c(-10000, y_0, y_1)
maxis <- c(10000, y_0, y_1)

model1 <-nls(y1~parab2(a,yy0,yy1,x1),start=startlist,na.action = na.exclude,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
parameters1 <- summary(model1)[["parameters"]]
vals <- parab2(a=parameters1[1], yy0=parameters1[2], yy1=parameters1[3], seq(0,1.,0.05))


lines(seq(0,1.,0.05), vals, type='l')



x1 <- tabmoy$Semprop1
y1 <- tabmoy$YEsp1
y_0 <- tabmoy$YEsp1[1]
y_1 <- tabmoy$YEsp1[7]


startlist=list(a=1., yy0=y_0, yy1=y_1)
minis <- c(-10000, y_0, y_1)
maxis <- c(10000, y_0, y_1)

model1 <-nls(y1~parab2(a,yy0,yy1,x1),start=startlist,na.action = na.exclude,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
parameters1 <- summary(model1)[["parameters"]]
vals <- parab2(a=parameters1[1], yy0=parameters1[2], yy1=parameters1[3], seq(0,1.,0.05))


lines(seq(0,1.,0.05), vals, type='l', col=2)



x1 <- tabmoy$Semprop1
y1 <- tabmoy$YEsp2
y_0 <- tabmoy$YEsp2[1]
y_1 <- tabmoy$YEsp2[7]


startlist=list(a=1., yy0=y_0, yy1=y_1)
minis <- c(-10000, y_0, y_1)
maxis <- c(10000, y_0, y_1)

model1 <-nls(y1~parab2(a,yy0,yy1,x1),start=startlist,na.action = na.exclude,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
parameters1 <- summary(model1)[["parameters"]]
vals <- parab2(a=parameters1[1], yy0=parameters1[2], yy1=parameters1[3], seq(0,1.,0.05))


lines(seq(0,1.,0.05), vals, type='l', col=4)





#poly3??

x1 <- tabmoy$Semprop1
y1 <- tabmoy$Ytot
y_0 <- tabmoy$Ytot[1]
y_1 <- tabmoy$Ytot[7]


startlist=list(a=1000., b=-2000., yy0=y_0, yy1=y_1)
minis <- c(-100000, -100000, y_0, y_1)
maxis <- c(100000, 100000, y_0, y_1)

model1 <-nls(y1~parab3(a,b,yy0,yy1,x1),start=startlist,na.action = na.exclude,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
parameters1 <- summary(model1)[["parameters"]]
vals <- parab3(a=parameters1[1], b=parameters1[2],yy0=parameters1[3], yy1=parameters1[4], seq(0,1.,0.05))


lines(seq(0,1.,0.05), vals, type='l')

#marche pas??

















#utiliser smoothspline??

model2 <- smooth.spline(x1, y1)
lines(model2, col=4)
integrate(model2,0,1)

predict(model2,seq(0.,1.,0.01)) #cree une fonction empirique pour pouvoir l'intgrer
f <- function(x) predict(model2,newdata=x)

# perform integration
integrate(f,0,1)
#marche pas mais devrait pouvoir s'en sortir??







# plot des donnees brutes en fonction de densite de semis
plot(dtoto$densite1, dtoto$Ytot, ylim=c(0,1800))
points(dtoto$densite1, dtoto$YEsp1,col=2)
points(dtoto$densite1, dtoto$YEsp2,col=4)


# plot des donnees brutes en fonction de proportion de semis
plot(dtoto$Semprop1, dtoto$Ytot, ylim=c(0,1800))
points(dtoto$Semprop1, dtoto$YEsp1,col=2)
points(dtoto$Semprop1, dtoto$YEsp2,col=4)


# plot des donnees brutes en fonction de proportion de recolte
plot(dtoto$Yprop1, dtoto$Ytot, ylim=c(0,2200))
points(dtoto$Yprop1, dtoto$YEsp1,col=2)
points(dtoto$Yprop1, dtoto$YEsp2,col=4)

# plot des donnees brutes pour l'N en fonction de la densite de semis
plot(dtoto$Semprop1, dtoto$QNtot, ylim=c(0,65))
points(dtoto$Semprop1, dtoto$QNupttot,col=2)
points(dtoto$Semprop1, dtoto$QNuptleg,col=4)


