dir <- choose.dir()
#"C:\simul\PMA18\toto"
#"C:\\devel\\l-egume\\legume\\multisim\\sorties"
# "H:\\simul\\sorties" 
#"C:\\devel\\grassland\\sorties\\mixture Fix-NonFix scenario1"
setwd(dir)
ls_files <- list.files(dir)


#recupere la liste des toto file names du dossier de travail
ls_toto <- ls_files[grepl('toto', ls_files)]

#creation du dataFrame dtoto et recup des info fichier
dtoto <- as.data.frame(t(as.data.frame(strsplit(ls_toto, '_'))))[,c(2,3,6,7,8,9,10)]
row.names(dtoto) <- 1: length(dtoto[,1])
names(dtoto) <- c('usm','lsystem','mix','damier','scenario','Mng', 'seed')
dtoto$name <- ls_toto
dtoto$seed <- substr(as.character(dtoto$seed), 1, 1)
dtoto$scenario <- substr(as.character(dtoto$scenario), 9, nchar(as.character(dtoto$scenario)))
dtoto$keysc <- paste(dtoto$scenario, dtoto$mix, dtoto$Mng)# ajout d'une cle unique par scenario
#dtoto$damier <- as.numeric(substr(as.character(dtoto$damier), 7, 7))

#split de dtoto et stockage dans une liste de scenatios
sp_dtoto <- split(dtoto, dtoto$keysc)




read_ltoto <- function(ls_toto)
{
  #recuperer les fichiers toto du dossier e travail dans une liste ltoto
  ltoto <- vector('list', length(ls_toto))
  names(ltoto) <- ls_toto
  
  for (i in 1:length(ls_toto))
  {
    name <- ls_toto[i]
    ltoto[[name]] <- read.table(name, header=T, sep=';')
  }
  ltoto
}



#didcols <- as.data.frame(residcols[seq(1,21,3), ])
#didcols$damier <- unique(as.character(dtoto$damier))
#write.csv(didcols, "didcols.csv", row.names=F)

#fichier d'id colones esp1 pour damier 8 (pour les cas ou bug/oubli dans les noms de colonnes)
didcols <- read.csv("C:/devel/l-egume/legume/multisim/didcols.csv") 




for (key in names(sp_dtoto))#key <- names(sp_dtoto)[1]
{
  ls_toto_paquet <- sp_dtoto[[key]]$name
  
  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  ltoto <- read_ltoto(ls_toto_paquet)
  #version locale du paquet de doto
  dtoto <- sp_dtoto[[key]]
  
  #recup du nom des esp
  mix <- strsplit(ls_toto_paquet[1], '_')[[1]][6] #suppose paquet fait par traitement
  esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
  esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'
  
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
  
  for (i in 1:length(ls_toto_paquet))#(ls_toto))
  {
    name <- ls_toto_paquet[i]
    damier <- strsplit(name, '_')[[1]][7]
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
    #esp <- 'Fix2'#'Fix3'#'Fix1'#'Fix' #pourquoi c'est ce nom au lieu de Fix???
    #esp2 <- 'nonFixSimTest'#'nonFix1'#'nonFix0' #pourquoi c'est ce nom au lieu de Fix???
    
    nomcol <- names(ltoto[[name]])
    if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    {
      idcols <- as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
    } else
    {
      idcols <- grepl(esp, nomcol) & !grepl(esp2, nomcol)#contient esp1 et pas esp2
    }
    
    dat1 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
    nb1 <- length(dat1)-2
    nbplt1 <- cbind(nbplt1, nb1)
    if (nb1>0)
    {
      MS1 <- as.matrix(dat1[dat1$V1=='MSaerien' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)
      ProdIaer1 <- rowSums(MS1) / s
      Nuptake_sol_leg <- as.matrix(dat1[dat1$V1=='Nuptake_sol',3:(3+nb1-1)], ncol=nb1)
      Nuptake_sol_leg <- as.numeric(rowSums(Nuptake_sol_leg) / s)
    } else
    {
      ProdIaer1 <- 0 #pas de plante de l'esp1
      Nuptake_sol_leg <- 0
    }
    YEsp1 <- cbind(YEsp1, sum(ProdIaer1))#cumul des 5 coupes
    QNuptleg <- cbind(QNuptleg, sum(Nuptake_sol_leg))
    
    #YEsp2
    if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    {
      idcols <- !as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
      idcols[1:2] <- FALSE #remet a faux les deux premieres colonnes
    } else
    {
      idcols <- grepl(esp2, nomcol)#contient esp2
    }
    
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

  #remise du dtoto locl dans sp_dtoto
  sp_dtoto[[key]] <- dtoto
  
}

#reagrege dtoto
dtoto <- do.call("rbind", sp_dtoto)


#dtoto$keysc <- paste(dtoto$scenario, dtoto$mix, dtoto$Mng)# ajout d'une cle unique par scenario
#?? tjrs meme denite de esp1??
#pb de nom de esp1 qui est dans celui de esp2!! -> bug du grepl

#hist(Yprop1)
#mean(Yprop1)





Build_AverageScTable <- function(dtoto, keysc)
{
    sc <- strsplit(keysc," ")[[1]][1]
    mix <- strsplit(keysc," ")[[1]][2]
    mng <- strsplit(keysc," ")[[1]][3]
    #recup d'un scenario
    #sc <- '1-1'
    #mix <- 'Fix2-nonFixSimTest'
    #mng <- 'Lusignan30IrrN2'
    
    res <- dtoto[dtoto$scenario==sc & dtoto$mix==mix & dtoto$Mng==mng, ]
    
    #calcul des valeurs moyennes
    x <- by(res$YEsp1, as.factor(res$densite1), mean)
    tabmoy <- data.frame(densite1=as.numeric(names(x)), YEsp1=as.numeric(x))
    x <- by(res$YEsp2, as.factor(res$densite1), mean)
    tabmoy$YEsp2 <- as.numeric(x)
    x <- by(res$Ytot, as.factor(res$densite1), mean)
    tabmoy$Ytot <- as.numeric(x)
    x <- by(res$YEsp1, as.factor(res$densite1), sd)
    tabmoy$YEsp1sd <- as.numeric(x)
    x <- by(res$YEsp2, as.factor(res$densite1), sd)
    tabmoy$YEsp2sd <- as.numeric(x)
    x <- by(res$Ytot, as.factor(res$densite1), sd)
    tabmoy$Ytotsd <- as.numeric(x)
    x <- by(res$Semprop1, as.factor(res$densite1), mean)
    tabmoy$Semprop1 <- as.numeric(x)
    
    x <- by(res$QNtot, as.factor(res$densite1), mean)
    tabmoy$QNtot <- as.numeric(x)
    x <- by(res$QNupttot, as.factor(res$densite1), mean)
    tabmoy$QNupttot  <- as.numeric(x)
    x <- by(res$QNuptleg, as.factor(res$densite1), mean)
    tabmoy$QNuptleg  <- as.numeric(x)
    x <- by(res$QNfix, as.factor(res$densite1), mean)
    tabmoy$QNfix  <- as.numeric(x)
    
    x <- by(res$Yprop1, as.factor(res$densite1), mean)
    tabmoy$Yprop1  <- as.numeric(x)
    
    tabmoy$mix <- mix
    tabmoy$sc <- sc
    tabmoy$Mng <- mng
    tabmoy$keysc <- paste(mix, sc, mng)
    
    tabmoy
}



# calcul des valeurs d'interet complementaires
CalcOpt <- function(modeltot , xx, yy)
{
  ## calculla proportion et la valeur max de l'overyielding
  pred <- predict(modeltot, seq(0,1,0.001))
  #xx <- tabmoy$Yprop1
  lintot <- lsfit(c(xx[1], xx[7]), c(yy[1], yy[7]))#c(tabmoy$Ytot[1], tabmoy$Ytot[7]))
  ylin <- lintot$coefficients[["Intercept"]] + seq(0,1,0.001)*lintot$coefficients[["X"]]
  
  diff_predlin <- pred$y - ylin
  
  idopt <- which(abs(diff_predlin ) == max(abs(diff_predlin )))
  propOpt <- pred$x[idopt]
  OverMax <- diff_predlin[idopt]
  
  #calcul du max de rendement absolu e de la prop correspodante
  Ytotmax <- max(pred$y)
  idmax <- which(pred$y == Ytotmax)[1] #le premier si plusieurs
  propMax <- pred$x[idmax]
  
  c(propOpt, OverMax, idopt, Ytotmax, propMax)
}


CalcPropactu50 <- function (modelesp1, modelesp2, idopt)
{
  #calcul prop a laquelle biomasse fait 50/50 (2 modeles se croisent) et prop debiomase a l'otimum d'overyielding
  pred1 <- predict(modelesp1, seq(0,1,0.001))
  pred2 <- predict(modelesp2, seq(0,1,0.001))
  delta <- abs(pred1$y-pred2$y)
  idmin <- which(delta == min(delta))
  propsowing50 <- pred1$x[idmin]
  propLegOtp <- pred1$y[idopt]/(pred1$y[idopt]+pred2$y[idopt])
  c(propLegOtp, propsowing50)
}


YtotvsProp <- function(tabmoy, Ymax=2200, nom="", optProp="sowing",visuplot=T, visutext=T, ...)
{
  ## calcul des composante de l'overyielding biomasse et fait un plot (visutext=visualisation des valeurs; visuplot=visulaisation des )
  

  #actual or sowing proportions?
  if (optProp=="sowing")
  {
    xx <- tabmoy$Semprop1
    labx <- 'Sowing proportion (Esp. 1)'
  }
  if (optProp=="actual")
  {
    xx <- tabmoy$Yprop1
    labx <- 'Actual proportion (Esp. 1)'
  }

  #calcul des fits des valeurs moyennes
  modeltot <- smooth.spline(xx, tabmoy$Ytot)
  inttot = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) - (tabmoy$Ytot[1]+tabmoy$Ytot[7])/2
  
  modelesp1 <- smooth.spline(xx, tabmoy$YEsp1)
  intesp1 = sum(predict(modelesp1, seq(0,1,0.001))$y*0.001) - (tabmoy$YEsp1[1]+tabmoy$YEsp1[7])/2
  
  modelesp2 <- smooth.spline(xx, tabmoy$YEsp2)
  intesp2 = sum(predict(modelesp2, seq(0,1,0.001))$y*0.001) - (tabmoy$YEsp2[1]+tabmoy$YEsp2[7])/2
  
  #cacul des autres indices
  ids <- CalcOpt(modeltot , xx, tabmoy$Ytot)
  propOpt <- ids[1]
  OverMax <- ids[2]
  Ytotmax <- ids[4]
  propYtotmax <- ids[5]
  ids1 <- CalcPropactu50(modelesp1, modelesp2, ids[3])
  propsowing50 <- ids1[2]
  propLegOtp <- ids1[1]
  
  #plot des valeur moyennes Ytot si option activee
  if (visuplot==T)
  {
    plot(xx, tabmoy$Ytot, ylim=c(0,Ymax), xlab=labx, ylab='Shoot biomass (g.m-2)', main=nom, ...)
    #segments(tabmoy$Semprop1, tabmoy$Ytot, tabmoy$Semprop1, tabmoy$Ytot+tabmoy$Ytotsd)
    #segments(tabmoy$Semprop1, tabmoy$Ytot, tabmoy$Semprop1, tabmoy$Ytot-tabmoy$Ytotsd)
    segments(xx[1], tabmoy$Ytot[1], xx[7], tabmoy$Ytot[7], lty=2)
    lines(modeltot)
    
    points(xx, tabmoy$YEsp1,col=2)
    segments(xx[1], tabmoy$YEsp1[1], xx[7], tabmoy$YEsp1[7], lty=2, col=2)
    lines(modelesp1, col=2)
    
    points(xx, tabmoy$YEsp2,col=4)
    segments(xx[1], tabmoy$YEsp2[1], xx[7], tabmoy$YEsp2[7], lty=2, col=4)
    lines(modelesp2, col=4)
    
  }
  
  if (visutext==T & visuplot==T)
  {
    text(0.15, 0.97*Ymax, paste('overY: ' ,round(inttot,2)))
    text(0.15, 0.93*Ymax, paste('Esp1: ' , round(intesp1,2)),col=2)
    text(0.15,0.89*Ymax, paste('Esp2: ' ,round(intesp2,2)),col=4)
  }
  
  #renvoie valeurs calculees
  res <- as.list(c(inttot, intesp1, intesp2, propOpt, OverMax, propsowing50, propLegOtp, Ytotmax, propYtotmax))
  names(res) <- c("inttot", "intesp1", "intesp2", "propOpt", "OverMax", "propsowing50", "propLegOtp", "Ytotmax", "propYtotmax")
  res

}


QNtotvsProp <- function(tabmoy, Ymax=100, nom="", optProp="sowing", visuplot=T, visutext=T, ...)
{
  ## calcul des composante de l'overyielding Ntot et fait un plot (visutext=visualisation des valeurs; visuplot=visulaisation des plots)
  
  
  #actual or sowing proportions?
  if (optProp=="sowing")
  {
    xx <- tabmoy$Semprop1
    labx <- 'Sowing proportion (Esp. 1)'
  }
  if (optProp=="actual")
  {
    xx <- tabmoy$Yprop1
    labx <- 'Actual proportion (Esp. 1)'
  }
  
  #calcul des fits des valeurs moyennes
  modeltot <- smooth.spline(xx, tabmoy$QNtot)
  intoverN = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) - (tabmoy$QNtot[1]+tabmoy$QNtot[7])/2
  intQNtot = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) 
  
  modelesp1 <- smooth.spline(xx, tabmoy$QNupttot)
  intNupt = sum(predict(modelesp1, seq(0,1,0.001))$y*0.001) 
  intFix = intQNtot-intNupt
  
  modeleg <- smooth.spline(xx, tabmoy$QNuptleg)
  intleg = sum(predict(modeleg, seq(0,1,0.001))$y*0.001) - (tabmoy$QNuptleg[1]+tabmoy$QNuptleg[7])/2
  
  #cacul des autres indices
  ids <- CalcOpt(modeltot , xx, tabmoy$QNtot)
  propOptN <- ids[1]
  OverMaxN <- ids[2]
  QNmax <- ids[4]
  propQNmax <- ids[5]
  
  if (visuplot==T)
  {
    plot(xx, tabmoy$QNtot, ylim=c(0,Ymax), xlab=labx, ylab='Plant N (g N.m-2)', main=nom, ...)
    segments(xx[1], tabmoy$QNtot[1], xx[7], tabmoy$QNtot[7], lty=2)
    lines(modeltot)
    
    points(xx, tabmoy$QNupttot,col=2)
    #segments(xx[1], tabmoy$QNupttot[1], xx[7], tabmoy$QNupttot[7], lty=2, col=2)
    lines(modelesp1, col=2)

    points(xx, tabmoy$QNuptleg,col=4)
    segments(xx[1], tabmoy$QNuptleg[1], xx[7], tabmoy$QNuptleg[7], lty=2, col=4)
    lines(modeleg, col=4)

  }
  
  if (visutext==T)
  {
    text(0.15,Ymax, paste(round(intoverN,2), '(over)'))    
    text(0.15,0.97*Ymax, paste(round(intFix,2), '(Fix)'),col=1)
    text(0.15,0.94*Ymax, paste(round(intNupt,2), '(Nupt)'),col=2)
    text(0.15,0.91*Ymax, paste(round(intleg,2), '(leg)'),col=4)    
  }
  
  res <- as.list(c(intoverN, intQNtot, intNupt, intFix, intleg, propOptN, OverMaxN, QNmax, propQNmax))
  names(res) <- c("intoverN", "intQNtot", "intNupt", "intFix", "intleg", "propOptN", "OverMaxN", "QNmax", "propQNmax")
  res
  
}

OverYvsAll <- function(ls_tabmoys, key, Ymax=300, nom="", optProp="sowing", visuplot=T,...)
{
  #key <- ls_keysc[20]
  #figure de tous les overyielding
  ls_keysc = names(ls_tabmoys)
  
  if (optProp=="sowing")
  { labx <- 'Sowing proportion (Esp. 1)'}
  if (optProp=="actual")
  { labx <- 'Actual proportion (Esp. 1)'}
  
  if (visuplot==T)
  {
    plot(-100, -100, ylim=c(-Ymax,Ymax), xlim=c(0,1), main=nom, xlab=labx, ylab='Overyieding (g.m-2)', ...)
    segments(0, 0, 1, 0, col=1)
  }
  
  resx <- NULL
  resy <- NULL
  
  for (keysc in ls_keysc)
  {
    #keysc <- ls_keysc[3]
    tabmoy <- ls_tabmoys[[keysc]]
    
    #xx <- tabmoy$Semprop1#tabmoy$Yprop1#
    yy <- tabmoy$Ytot
    #actual or sowing proportions?
    if (optProp=="sowing")
    {
      xx <- tabmoy$Semprop1
      labx <- 'Sowing proportion (Esp. 1)'
    }
    if (optProp=="actual")
    {
      xx <- tabmoy$Yprop1
      labx <- 'Actual proportion (Esp. 1)'
    }
    
    lintot <- lsfit(c(xx[1], xx[7]), c(yy[1], yy[7]))
    ylin <- lintot$coefficients[["Intercept"]] + xx*lintot$coefficients[["X"]]
    overY <- yy - ylin
    
    if (keysc != key)
    {
      if (visuplot==T)
      { points(xx, overY, pch=16, col='light grey') }
      resx <- cbind(resx,xx)
      resy <- cbind(resy,overY)
    } else
    {
      savexx <- xx
      saveyy <- overY
    }
  }
  if (visuplot==T)
  { points(savexx, saveyy, pch=16, col='blue', type='b')}
  resx <- cbind(resx,savexx)
  resy <- cbind(resy,saveyy)
  data.frame(x=as.numeric(resx), y=as.numeric(resy))#renvoie les points du graphique (avec 7 derniers=ceux de l'id selectionne)
}





#faire une foonction, construire liste moyenne
#en amont recup les valeurs de parametres... ou apres pour etre generique

# recup de la liste des scenario et rempli liste avec des tabmoys
ls_keysc <- unique(dtoto$keysc)
ls_tabmoys <- vector("list", length(ls_keysc))
names(ls_tabmoys) <- ls_keysc

for (keysc in ls_keysc)
{
  ls_tabmoys[[keysc]] <- Build_AverageScTable(dtoto, keysc)
}





keysc <- ls_keysc[3]
tabmoy <- ls_tabmoys[[keysc]]
nom <- paste(tabmoy$mix[1], tabmoy$sc[1], tabmoy$Mng[1])
YtotvsProp(tabmoy, nom=nom, optProp="sowing")
YtotvsProp(tabmoy, nom=nom, optProp="actual")

QNtotvsProp(tabmoy, nom=nom)
QNtotvsProp(tabmoy, nom=nom, optProp="actual")

OverYvsAll(ls_tabmoys, keysc, nom="", optProp="sowing")


#sauve en csv tableau agrege
#ecriture fichier
write.csv(dtoto, "dtoto7ok.csv", row.names=F)
tabmoys <- do.call("rbind", ls_tabmoys) #merge a list of data.frames - do.call equalent de map
write.csv(tabmoys, "tabmoys7ok.csv", row.names=F)

ls_tabmoys


tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211965&authkey=ADR9j8ZkMzAXz8I") #serie 1
#marche aussi en lecture directe!
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)

ls_tabmoys <- ls_tabmoys2
ls_keysc <- names(ls_tabmoys)

#ecriture du merge
tabmoys_m <-do.call("rbind", c(ls_tabmoys, ls_tabmoys2))
write.csv(tabmoys_m, "tabmoys_merge5.csv", row.names=F)

tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211970&authkey=ABgO3gxUEz19IBE")#merge file (1-2)
#figure de tous les overyielding
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)

tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211974&authkey=APtBFCEffHjVZb0")#merge file (1-2-3)
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)

tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211975&authkey=AMhhE6ZN_dcGHzU")#merge file (1-2-3-4)
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)

tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211979&authkey=ALNOA6Fw6GMl6Lw")#merge file (1-2-3-4-5)
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)

tabmoys2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214039&authkey=APnHPpEOCjBL4Qg")#merge file (1-2-3-4-5-6)
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)


ls_tabmoys <- ls_tabmoys2
ls_keysc <- names(ls_tabmoys)




#construction et lecture des fichiers moyens de synthese
tabmoys2 <- read.csv("C:/simul/PMA18/tabmoys2.csv")#
tabmoys3 <- read.csv("C:/simul/PMA18/tabmoys3.csv")#
tabmoys4 <- read.csv("C:/simul/PMA18/tabmoys4.csv")#
tabmoys5 <- read.csv("C:/simul/PMA18/tabmoys5.csv")#
tabmoys6 <- read.csv("C:/simul/PMA18/tabmoys6.csv")#
tabmoys7 <- read.csv("C:/simul/PMA18/tabmoys7.csv")#
ls_tabmoys2 <- split(tabmoys2, tabmoys2$keysc)
ls_tabmoys3 <- split(tabmoys3, tabmoys3$keysc)
ls_tabmoys4 <- split(tabmoys4, tabmoys4$keysc)
ls_tabmoys5 <- split(tabmoys5, tabmoys5$keysc)
ls_tabmoys6 <- split(tabmoys6, tabmoys6$keysc)
ls_tabmoys7 <- split(tabmoys7, tabmoys7$keysc)

tabmoys_m <-do.call("rbind", c(ls_tabmoys2, ls_tabmoys3, ls_tabmoys4, ls_tabmoys5,ls_tabmoys6, ls_tabmoys7))
write.csv(tabmoys_m, "tabmoys_merge2-7.csv", row.names=F)

tabmoys_m <-do.call("rbind", c(ls_tabmoys2, ls_tabmoys3, ls_tabmoys4, ls_tabmoys5))
write.csv(tabmoys_m, "tabmoys_merge2-5.csv", row.names=F)

tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214078&authkey=ABMcPNVS8FNlM9E")#'tabmoys_merge2-7.csv'
tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214080&authkey=AHUcqW3ERDVAd5s")#'tabmoys_merge2-7.csv'
tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214082&authkey=AN-6wl3mf7PEeGM")#'tabmoys_merge2-7.csv'
ls_tabmoys2 <- split(tabmoys_m2, tabmoys_m2$keysc)

tabmoys_m2 


tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214078&authkey=ABMcPNVS8FNlM9E")#'tabmoys_merge2-7.csv'
#tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214080&authkey=AHUcqW3ERDVAd5s")#'tabmoys_merge2-5.csv'
tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214080&authkey=AHUcqW3ERDVAd5s")#'tabmoys_merge2-5.csv'


spmn <- split(tabmoys_m2, tabmoys_m2$Mng) #decoupe par niveau d'N
ls_tabmoys <- vector("list", length(spmn))
names(ls_tabmoys) <- names(spmn)

for (traitN in names(spmn))
{
  spmn[[traitN]]$keysc <- as.factor(as.character(paste(spmn[[traitN]]$mix , spmn[[traitN]]$sc)))#spmn[[traitN]]$keysc)) #! enelve Mng de la cle!!
  ls_tabmoys[[traitN]] <- split(spmn[[traitN]], spmn[[traitN]]$keysc)
}

length(ls_tabmoys[[1]])
names(ls_tabmoys) <- c('0N', '300N', '120N')

ls_keysc <- names(ls_tabmoys[[1]])#c(names(ls_tabmoys[[1]]), names(ls_tabmoys[[2]]), names(ls_tabmoys[[3]]))

#noms contiennent traitement N!!
# a refaire avant le split

names(ls_tabmoys[['120N']])

ls_tabmoys[['120N']][["nonFixSimTest-nonFixSimTest 55-55"]]#neutral situation = not neutral!!??



for (key in names(ls_tabmoys))
{
  OverYvsAll(ls_tabmoys, key)
  locator(1)
}
  



#ou le fait sans boucle!
  

#!!! attention: reussi sur onedrive en recuperant le lien onedrive pour 'embeded' et en remplacant dans l'url 'embeded' par download (cf: https://stackoverflow.com/questions/29579782/reading-onedrive-files-to-r)


#ecuperation des valeurs de parametre
library(readxl)

mix <- "Fix2-nonFixSimTest"

params <- vector("list", 2)
names(params) <- strsplit(mix, '-')[[1]]

path_param <- "C:/devel/l-egume/legume/input/liste_scenarios.xls"
params[[1]] <- read_excel(path_param, sheet = names(params)[1])
params[[2]] <- read_excel(path_param, sheet = names(params)[2])


#faire le lien entre les valeurs de difference de parametre et les scenaios

#liste unique des scenarios
dt <- as.data.frame(do.call("rbind", strsplit(ls_keysc, " ")))
lsc <- as.data.frame(do.call("rbind", strsplit(as.character((dt$V1)), "-")))
lsp <- as.data.frame(do.call("rbind", strsplit(as.character((dt$V2)), "-")))
dparams <- cbind(lsc, lsp)
names(dparams) <-c("esp2", "esp1","sc1", "sc2" )#c("sc1", "sc2", "esp2", "esp1")#verif si numero espece pas inverses?


#calcul des differences de valeurs de parametres
res <- NULL
for (i in 1:length(dparams$sc1))
{
  p1 <- params[[as.character(dparams$esp1[i])]][params[[as.character(dparams$esp1[i])]]$id_scenario == dparams$sc1[i] , c(2,3,4,5,6)]
  p2 <- params[[as.character(dparams$esp2[i])]][params[[as.character(dparams$esp2[i])]]$id_scenario == dparams$sc2[i] , c(2,3,4,5,6)]
  res <- rbind(res, p1-p2)
}

dparams <- cbind(dparams, res)
dparams$keysc <- ls_keysc

resnorm <- res
names(resnorm) <- c("normq", "normLen", "normVmax2", "normRUE", "normMaxFix")
resnorm$normq[resnorm$normq>0] <- 1
resnorm$normq[resnorm$normq<0] <- -1
resnorm$normLen[resnorm$normLen>0] <- 1
resnorm$normLen[resnorm$normLen<0] <- -1
resnorm$normVmax2[resnorm$normVmax2>0] <- 1
resnorm$normVmax2[resnorm$normVmax2<0] <- -1
resnorm$normRUE[resnorm$normRUE=='0.2'] <- -0.5
resnorm$normRUE[resnorm$normRUE=='0.6'] <- -1
resnorm$normMaxFix[resnorm$normMaxFix<0] <- 1

dparams <- cbind(dparams, resnorm)
#aller faire le lien avec les valeurs des boutons!


dparams[dparams$RUE>0.3 & ,]

normq <- 0
normLen <- 0
normVmax2 <- 0
normRUE <- 0

seletedkey <- dparams[dparams$normq == normq & dparams$normLen == normLen & dparams$normVmax2 == normVmax2 & dparams$normRUE == normRUE, "keysc"]





#difference semble pas coller?

#comment rentrer les interfaces (ecart de temperature opt?)


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




#tests d'ouverture des url de fichiers
load(url("https://sourcesup.renater.fr/frs/download.php/latestfile/2174/ls_tabmoys.Rdata"))
#marche pas???

https://sourcesup.renater.fr/frs/download.php/file/5630/ls_tabmoys.Rdata

download.file(url("https://sourcesup.renater.fr/frs/download.php/latestfile/2174/ls_tabmoys.Rdata"), "test.Rdata", method="wininet")
download.file("https://sourcesup.renater.fr/frs/download.php/latestfile/2174/ls_tabmoys.Rdata", "test.Rdata", method="libcurl")
load("test.Rdata")

#le fichier telecharge est une page html (probablement de login)-> secured urls


#install.packages("curl")
library(curl)
#install.packages("RCurl")
library(RCurl)

URL <- "https://sourcesup.renater.fr/frs/download.php/latestfile/2174/ls_tabmoys.Rdata"#"https://d396qusza40orc.cloudfront.net/getdata%2Fdata%2Fss06hid.csv"
x <- getURL(URL)
curl_download(url =URL ,destfile="TEST2.Rdata",quiet=FALSE, mode="wb")


curl_download("https://sourcesup.renater.fr/account/login.php?triggered=1&return_to=%2Ffrs%2Fdownload.php%2Flatestfile%2F2174%2Fls_tabmoys.Rdata", destfile="TEST3.Rdata",quiet=FALSE, mode="wb")

#f using RCurl you get an SSL error on the GetURL() function then set these options before GetURL(). This will set the CurlSSL settings globally. 
#continuer a creuser...



#adresse sur site unite
download.file("https://www6.nouvelle-aquitaine-poitiers.inra.fr/urp3f_admin/content/download/3346/33241/file/test.Rdata", "testinra.Rdata", method="auto")
load("testinra.Rdata")
#meme probleme depuis l'unite...

URL2 <- "https://www6.nouvelle-aquitaine-poitiers.inra.fr/urp3f_admin/content/download/3346/33241/file/test.Rdata"
x <- getURL(URL2)


library(httr)
URL2 <- "https://www6.nouvelle-aquitaine-poitiers.inra.fr/urp3f_admin/content/download/3346/33241/file/test.Rdata"
URL <- "https://sourcesup.renater.fr/account/login.php?triggered=1&return_to=%2Ffrs%2Fdownload.php%2Flatestfile%2F2174%2Fls_tabmoys.Rdata"
GET(URL, authenticate('glouarn', 'Tyt2omt_'), write_disk("TEST4.Rdata"), timeout(60))
load("TEST4.Rdata")


#meme probleme:: pas le fichier que je veux



#test depuis onedrive (public)

download.file("https://1drv.ms/u/s!AnATzWXkvRzDjyplVQJmGBZEzdsi", destfile= "test5.Rdata", method="auto")
load("test5.Rdata")

download.file("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211962&authkey=AP7KauBzwVGuU-w", destfile= "test6.Rdata", method="auto")
#fichier bien telecharge!
load("test6.Rdata")# arrive pas a me le charger!!


#rq: script en deux etapes marche pas pour un usage en ligne: peut pas faire le download local!
#essayer de sauver le fichier en tableau csv aggrege (que l'on split apres ouverture)... plus simple


download.file("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211963&authkey=AH1102XKYpbg654", destfile= "test7.csv", method="auto")
read.csv("test7.csv")


read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%211963&authkey=AH1102XKYpbg654")
#marche aussi en lecture directe!














### test du script 
tabmoys_m2 <- read.csv("https://onedrive.live.com/download?cid=C31CBDE465CD1370&resid=C31CBDE465CD1370%214082&authkey=AN-6wl3mf7PEeGM")#'tabmoys_merge2-7.csv'


spmn <- split(tabmoys_m2, tabmoys_m2$Mng) #decoupe par niveau d'N
ls_tabmoys <- vector("list", length(spmn))
names(ls_tabmoys) <- names(spmn)

for (traitN in names(spmn))
{
  spmn[[traitN]]$keysc <- as.factor(as.character(paste(spmn[[traitN]]$mix , spmn[[traitN]]$sc)))#spmn[[traitN]]$keysc)) #! enelve Mng de la cle!!
  ls_tabmoys[[traitN]] <- split(spmn[[traitN]], spmn[[traitN]]$keysc)

}

length(ls_tabmoys[[1]])
names(ls_tabmoys) <- c('0N', '300N', '120N')#c('0N',  '120N')#

ls_keysc <- names(ls_tabmoys[[1]])#clesscenario sans doublon management





library(readxl)

mix <- "Fix2-nonFixSimTest"
params <- vector("list", 2)
names(params) <- strsplit(mix, '-')[[1]]

path_param <- "C:/devel/l-egume/legume/input/liste_scenarios.xls"
params[[1]] <- read_excel(path_param, sheet = names(params)[1])
params[[2]] <- read_excel(path_param, sheet = names(params)[2])


#faire le lien entre les valeurs de difference de parametre et les scenaios

#liste unique des scenarios
dt <- as.data.frame(do.call("rbind", strsplit(ls_keysc, " ")))
lsc <- as.data.frame(do.call("rbind", strsplit(as.character((dt$V1)), "-")))
lsp <- as.data.frame(do.call("rbind", strsplit(as.character((dt$V2)), "-")))
dparams <- cbind(lsc, lsp)
names(dparams) <-c("esp2", "esp1","sc1", "sc2" )#c("sc1", "sc2", "esp2", "esp1")#verif si numero espece pas inverses?


#calcul des differences de valeurs de parametres
res <- NULL
for (i in 1:length(dparams$sc1))
{
  p1 <- params[[as.character(dparams$esp1[i])]][params[[as.character(dparams$esp1[i])]]$id_scenario == dparams$sc1[i] , c(2,3,4,5,6)]
  p2 <- params[[as.character(dparams$esp2[i])]][params[[as.character(dparams$esp2[i])]]$id_scenario == dparams$sc2[i] , c(2,3,4,5,6)]
  res <- rbind(res, p1-p2)
}

dparams <- cbind(dparams, res)
dparams$keysc <- ls_keysc

resnorm <- res
names(resnorm) <- c("normq", "normLen", "normVmax2", "normRUE", "normMaxFix")
resnorm$normq[resnorm$normq>0] <- 1
resnorm$normq[resnorm$normq<0] <- -1
resnorm$normLen[resnorm$normLen>0] <- 1
resnorm$normLen[resnorm$normLen<0] <- -1
resnorm$normVmax2[resnorm$normVmax2>0] <- 1
resnorm$normVmax2[resnorm$normVmax2<0] <- -1
resnorm$normRUE[resnorm$normRUE=='0.2'] <- -0.5
resnorm$normRUE[resnorm$normRUE=='0.6'] <- -1
resnorm$normMaxFix[resnorm$normMaxFix<0] <- 1#inverse

dparams <- cbind(dparams, resnorm)


#inputs
nfert_level <- '300N'
bw_q <- 0
bw_Len <- 0
bw_Vmax2 <- 0
bw_RUE <- 0
bw_Fix <- 0

seletedkey <- dparams[dparams$normq == bw_q & dparams$normLen == bw_Len & dparams$normVmax2 == bw_Vmax2 & dparams$normRUE == bw_RUE & dparams$normMaxFix == bw_Fix, "keysc"]
# 


seletedkey <- dparams[dparams$normq == 0 & dparams$normLen == 0 & dparams$normVmax2 == 0 & dparams$normRUE == 0 & dparams$normMaxFix == 0, "keysc"]


tabmoy <- ls_tabmoys[[nfert_level]][[seletedkey]]
YtotvsProp(tabmoy, nom="", optProp="actual", visutext=F, cex.lab=1.2)
OverYvsAll(ls_tabmoys[[nfert_level]], seletedkey, nom="", optProp="sowing", cex.lab=1.2, Ymax=400)

#sowing prop
bw_Fix <- 0#1#
layout(matrix(1:3,1,3))
for (nfert_level in c('0N','120N','300N'))
{
  seletedkey <- dparams[dparams$normq == bw_q & dparams$normLen == bw_Len & dparams$normVmax2 == bw_Vmax2 & dparams$normRUE == bw_RUE & dparams$normMaxFix == bw_Fix, "keysc"]
  isfix <- grepl('Fix2-', names(ls_tabmoys[[nfert_level]]))
  if (bw_Fix == 0) 
  {ls_tab <- ls_tabmoys[[nfert_level]][!isfix]
  } else 
  {ls_tab <- ls_tabmoys[[nfert_level]][isfix]}
  OverYvsAll(ls_tab, seletedkey, nom=paste('sowing',nfert_level, bw_Fix), optProp="sowing", cex.lab=1.2, Ymax=500)
}


#actual prop
bw_Fix <- 1#0#
layout(matrix(1:3,1,3))
for (nfert_level in c('0N','120N','300N'))
{
  seletedkey <- dparams[dparams$normq == bw_q & dparams$normLen == bw_Len & dparams$normVmax2 == bw_Vmax2 & dparams$normRUE == bw_RUE & dparams$normMaxFix == bw_Fix, "keysc"]
  isfix <- grepl('Fix2-', names(ls_tabmoys[[nfert_level]]))
  if (bw_Fix == 0) 
  {ls_tab <- ls_tabmoys[[nfert_level]][!isfix]
  } else 
  {ls_tab <- ls_tabmoys[[nfert_level]][isfix]}
  res <- OverYvsAll(ls_tab, seletedkey, nom=paste('actual',nfert_level, bw_Fix), optProp="actual", cex.lab=1.2, Ymax=500)
}






#sowing histo
bw_Fix <- 0#1#
layout(matrix(1:3,1,3))
for (nfert_level in c('0N','120N','300N'))
{
  seletedkey <- dparams[dparams$normq == bw_q & dparams$normLen == bw_Len & dparams$normVmax2 == bw_Vmax2 & dparams$normRUE == bw_RUE & dparams$normMaxFix == bw_Fix, "keysc"]
  isfix <- grepl('Fix2-', names(ls_tabmoys[[nfert_level]]))
  if (bw_Fix == 0) 
  {ls_tab <- ls_tabmoys[[nfert_level]][!isfix]
  } else 
  {ls_tab <- ls_tabmoys[[nfert_level]][isfix]}
  res <- OverYvsAll(ls_tab, seletedkey, nom=paste('sowing',nfert_level, bw_Fix), optProp="sowing", visuplot = F)
  x <- hist(res$y, main=paste('sowing',nfert_level, bw_Fix), col=c(rep(0,20),rep(2,20)), xlim=c(-350,350), breaks = seq(-1000,1000,50), ylim=c(0,400), xlab="OverY")
  #segments(0,0,0,400, col=2)
  propsup0 <- sum(x$counts[21:40])/sum(x$counts)*100
  propsup50 <- sum(x$counts[22:40])/sum(x$counts)*100
  text(200, 350, round(propsup0, 2))
  text(200, 320, round(propsup50, 2))
}





#actual histo
bw_Fix <- 0#1#
layout(matrix(1:3,1,3))
for (nfert_level in c('0N','120N','300N'))
{
  seletedkey <- dparams[dparams$normq == bw_q & dparams$normLen == bw_Len & dparams$normVmax2 == bw_Vmax2 & dparams$normRUE == bw_RUE & dparams$normMaxFix == bw_Fix, "keysc"]
  isfix <- grepl('Fix2-', names(ls_tabmoys[[nfert_level]]))
  if (bw_Fix == 0) 
  {ls_tab <- ls_tabmoys[[nfert_level]][!isfix]
  } else 
  {ls_tab <- ls_tabmoys[[nfert_level]][isfix]}
  res <- OverYvsAll(ls_tab, seletedkey, nom=paste('actual',nfert_level, bw_Fix), optProp="actual", visuplot = F)
  x <- hist(res$y, main=paste('actual',nfert_level, bw_Fix), col=c(rep(0,20),rep(2,20)), xlim=c(-350,350), breaks = seq(-1000,1000,50), ylim=c(0,400), xlab="OverY")
  #segments(0,0,0,400, col=2)
  propsup0 <- sum(x$counts[21:40])/sum(x$counts)*100
  propsup50 <- sum(x$counts[22:40])/sum(x$counts)*100
  text(200, 350, round(propsup0, 2))
  text(200, 320, round(propsup50, 2))
}



#pour faire la comparaison avec modele neutre 
#ligne : caluler fit sur derniere serie / faire les calcul de difference en fontion des proportion et du modele spleen fite
#recalule OverY-neutral




#passer tous les tab et faire de data.frame de synthse des indicateurs

# A faire: recuperer un plan: (e.g.) 
#retirer les bons id avec un %in% dans dparams

#de la meme facon faire graph triangle...

#ou graph sctter 3D avec pojetions comme dans leaf spectrum?


resindice <- vector("list",2)
names(resindice) <- c('actual','sowing')
aa <-  vector("list",3)
names(aa) <- c('0N','120N','300N')
resindice[['actual']] <- aa
resindice[['sowing']] <- aa

for (outype in c('actual','sowing'))
{
  for (nfert_level in c('0N','120N','300N'))
  {
  
    resall <- NULL
    for (seletedkey  in names(ls_tabmoys[[nfert_level]]))
    {
      tabmoy <- ls_tabmoys[[nfert_level]][[seletedkey]]
      res <- as.data.frame(YtotvsProp(tabmoy, nom="", optProp=outype, visutext=F, visuplot=F, cex.lab=1.2))
      res$nfert <- nfert_level
      res$keysc <- seletedkey
      resall <- rbind(resall, res)
    }
    
    #y ajouter le decodage des valeurs de parametes de dparams
    #length(dparams$keysc)
    #length(resall$keysc)
    #?? ps la meme longueur?? manque 2 simul!
    
    
    resall2 <- NULL
    for (key  in resall$keysc)
    {
      res <- dparams[dparams$keysc == key, -10]
      resall2 <- rbind(resall2, res)
    }
    
    resindice[[outype]][[nfert_level]] <- cbind(resall, resall2)
  }
}

nfert_level <- '120N'#"300N"#'0N'#
outype <- "actual"
isfix <- grepl('Fix2-', resindice[[outype]][[nfert_level]]$keysc)

resindice[[outype]][[nfert_level]][isfix,]#pour les asso
resindice[[outype]][[nfert_level]][!isfix,]#pour les purs




#install.packages("scatterplot3d") # Install
library("scatterplot3d")
scatterplot3d(resindice[[outype]][[nfert_level]][isfix,c("normq","normLen","inttot")], pch = 16, color="steelblue")



#plot triangle
library(ggtern)
require(akima)#pour interpollation


# The biggest difficulty in the making of a ternary plot is to transform triangular coordinates into cartesian coordinates, here is a small function to do so:
tern2cart <- function(coord)
{
  coord[1]->x
  coord[2]->y
  coord[3]->z
  x+y+z -> tot
  x/tot -> x  # First normalize the values of x, y and z
  y/tot -> y
  z/tot -> z
  (2*y + z)/(2*(x+y+z)) -> x1 # Then transform into cartesian coordinates
  sqrt(3)*z/(2*(x+y+z)) -> y1
  return(c(x1,y1))
}



nfert_level <-'0N'#"300N"#'120N'# 
outype <- "sowing"#"actual"#
normRUE <- 0.0
asso <- 'pur'#'asso'#
var <- "propOpt"#"propsowing50"#"OverMax"#"inttot"#
breaks_ <- c(-300,-200,-100,0,50,100,150,200)#"inttot"/actual
breaks_ <- c(-400,-250,-100,0,100,200,300,400)#"inttot"/actual
breaks_ <- c(0,0.1,0.3,0.4,0.5,0.6,0.7,0.9)#"propOpt / sowing

layout(matrix(1:3,1,3))
for (nfert_level in c('0N','120N','300N'))
{
  
  #mise en forme donnee
  isfix <- grepl('Fix2-', resindice[[outype]][[nfert_level]]$keysc)
  islevRUE <- resindice[[outype]][[nfert_level]]$normRUE == 0.0
  
  if (asso=='asso')
  {
  plan <- resindice[[outype]][[nfert_level]][isfix & islevRUE, c("normq","normLen","normVmax2",var)]  
  } else
  {plan <- resindice[[outype]][[nfert_level]][!isfix & islevRUE, c("normq","normLen","normVmax2",var)]}
  
  df <- (plan[,c(1,2,3)]+1)/2*100 #passge en %
  #df[26,] <- df[26,]+0.0001 #pour le point zero qui bug ->non! normal que null d'y soit pas (serait au meme endroit que 1,1,1)-> repérer nul sur l'échelle avec une fleche
  #!! en fait 0,0,0 ps nul c'est 0.5,0.5,0.5) qui est null; mais aussi au meme endroit que 1,1,1!! probleme de cette visu
  #df <- df[-26,]#enleve le zero (neutre)
  retire <- c(-6,-10,-12,-21,-22,-23,-25,-26)#doublons a retirer (en dur avec modele nul)
  df <- df[retire,]#enleve le zero (neutre)
  
  q <- df[,1]
  Len <- df[,2]
  Vmax <- df[,3] 
  OverY <- plan[retire,4]#c(500,2324.90,2551.44,1244.50, 551.22,-644.20,-377.17,-100, 2493.04) 
  
  df <- data.frame(q, Len, Vmax, OverY)
  
  coords <- as.data.frame(t(apply(df,1,tern2cart)))
  names(coords) <- c('X', 'Y')
  #plot(coords)
  #pas utilise dans le tri finalement
  cols <- heat.colors(100, alpha = 1)
  #idcols <- cols[round((OverY+250)/(150+250)*100, 0)]#cols[round((OverY+350)/(350+350)*100, 0)]#
  #plot(coords, col=idcols, pch=16)
  
  coords$OverY <- plan[retire,4]
  #library(lattice)
  #levelplot(OverY~X*Y, coords)
  #coords
  
  
  #plot de la figure triangle
  # First create the limit of the ternary plot:
  titre <- paste(outype, nfert_level, asso, normRUE)
  plot(NA,NA,xlim=c(-0.1,1.1),ylim=c(-0.1,sqrt(3)/2)+0.1,asp=1,bty="n",axes=F,xlab="",ylab="",main=titre)
  
  
  resolution <- 0.001
  
  #interp(tern[,1],tern[,2],z=d, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid
  interp(coords[,1],coords[,2],z=coords[,3], xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid
  
  # And then plot:
  image(tern.grid,breaks=breaks_,col=rev(heat.colors(7)),add=T)
  contour(tern.grid,levels=breaks_,add=T)
  points(coords[,c(1,2)],pch=19)
  
  #graph personalise super! a finir (gradient de couleur negatif...)
  #titre, nom des parametres
  #fleches
  segments(0,0,0.5,sqrt(3)/2)
  segments(0.5,sqrt(3)/2,1,0)
  segments(1,0,0,0)
  text(0.5,(sqrt(3)/2),names(df)[3], pos=3)#"c"
  text(-0.08,0.05,names(df)[1], pos=1)#"a"
  text(1.08,0.05,names(df)[2], pos=1)#"b"
  
  
  legend(-0.1, 0.9, rev(c(-200,-100,0,50,100,150,200)), fill = heat.colors(7), col = rev(heat.colors(7)), border = NULL, bty="n", cex=0.8)
  text(0., 0.95, var)
  points(0.5, 0.288, col=4, pch=16, cex=1.3)#modele null

}






#plot de correlations
pairs(resindice[['sowing']][[nfert_level]][,c(1,4,5)])

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


dat <- rbind(resindice[['sowing']][['0N']][,c(1,4,5)], resindice[['sowing']][['120N']][,c(1,4,5)], resindice[['sowing']][['300N']][,c(1,4,5)])
pairs(dat[dat$OverMax<410,], lower.panel = panel.cor, upper.panel = panel.smooth)



pairs(resindice[['actual']][[nfert_level]][,c(1,4,5,6,7,9)])
pairs(resindice[['actual']][[nfert_level]][,c(1,4,5,6,7,9)])

plot(resindice[['sowing']][[nfert_level]]$inttot, resindice[['actual']][[nfert_level]]$inttot)
#??

#pb dans la boucme??


ggtern(df,aes(df[,1],df[,2],df[,3])) + 
  geom_point(aes(fill=plan[retire,4]))



#retirer de df les points doublons

df <- cbind(df,coords)

#faire une cle unique pour tier et enlever les doublons

#df$charx <- as.character(coords['X'])






library(ggtern)


#Orignal Data as per Question
q <- df[,1]#c(0.1, 0.5,0.5, 0.6, 0.2, 0          , 0         , 0.004166667, 0.45) 
Len <- df[,2]#c(0.75,0.5,0  , 0.1, 0.2, 0.951612903,0.918103448, 0.7875     , 0.45)
Vmax <- df[,3]#c(0.15,0  ,0.5, 0.3, 0.6, 0.048387097,0.081896552, 0.208333333, 0.10) 
OverY <- plan[retire,4]#c(500,2324.90,2551.44,1244.50, 551.22,-644.20,-377.17,-100, 2493.04) 
df <- data.frame(q, Len, Vmax, OverY)

#For labelling each point.
df$id <- 1:nrow(df)
ggtern(data=df,aes(x=Vmax,y=q,z=Len),aes(x,y,z)) + 
  stat_density_tern(geom="polygon",color='black',
                    n=400,h=0.75,expand = 1.1,
                    base='identity',
                    aes(fill   = ..level..,weight = d),
                    na.rm = TRUE) + 
  geom_point(color="black",size=5,shape=21) +
  geom_text(aes(label=id),size=3) + 
  scale_fill_gradient(low="yellow",high="red") + 
  scale_color_gradient(low="yellow",high="red") + 
  theme_rgbw() + 
  theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(  title= "Ternary Plot and Filled Contour",fill = "Value, OverY")




#seulement 19 points et as 27???
#7 doublons dont je ne comprends pas trop l'origine?
cbind(plan, coords, df)
#comment lit on le graph?
#is a barycentric plot on three variables which sum to a constant!!! -> ma somme n'est pas constante!
#comment fait tilman?
#2 variables ou cube....

#piste : diagramme chromatique (aussi en triangle) -> transformation differente
#comme pour les couleur la somme des range est constante


ggtern(plan,aes(plan[,1],plan[,2],plan[,3])) + 
  geom_interpolate_tern(aes(value=d,fill=..level..),
                        binwidth=500,
                        colour="white") +
  geom_point(aes(fill=plan[,4])),color="black",shape=21,size=3) + 
  scale_fill_gradient(low="yellow",high="red") +
  theme(legend.position=c(0,1),legend.justification=c(0,1)) + 
  labs(fill="Value, d")


ggtern(plan,aes(plan[,1],plan[,2],plan[,3])) + 
  geom_point(aes(fill=plan[,4]))








#exemple1
a<- c (0.1, 0.5, 0.5, 0.6, 0.2, 0, 0, 0.004166667, 0.45) 
b<- c (0.75,0.5,0,0.1,0.2,0.951612903,0.918103448,0.7875,0.45)
c<- c (0.15,0,0.5,0.3,0.6,0.048387097,0.081896552,0.208333333,0.1) 
d<- c (500,2324.90,2551.44,1244.50, 551.22,-644.20,-377.17,-100, 2493.04) 
df<- data.frame (a, b, c)


# First create the limit of the ternary plot:
plot(NA,NA,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,bty="n",axes=F,xlab="",ylab="")
segments(0,0,0.5,sqrt(3)/2)
segments(0.5,sqrt(3)/2,1,0)
segments(1,0,0,0)
text(0.5,(sqrt(3)/2),"c", pos=3)
text(0,0,"a", pos=1)
text(1,0,"b", pos=1)

# The biggest difficulty in the making of a ternary plot is to transform triangular coordinates into cartesian coordinates, here is a small function to do so:
tern2cart <- function(coord)
{
  coord[1]->x
  coord[2]->y
  coord[3]->z
  x+y+z -> tot
  x/tot -> x  # First normalize the values of x, y and z
  y/tot -> y
  z/tot -> z
  (2*y + z)/(2*(x+y+z)) -> x1 # Then transform into cartesian coordinates
  sqrt(3)*z/(2*(x+y+z)) -> y1
  return(c(x1,y1))
}

# Apply this equation to each set of coordinates
t(apply(df,1,tern2cart)) -> tern

# Intrapolate the value to create the contour plot
resolution <- 0.001
require(akima)
interp(tern[,1],tern[,2],z=d, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid

# And then plot:
image(tern.grid,breaks=c(-1000,0,500,1000,1500,2000,3000),col=rev(heat.colors(6)),add=T)
contour(tern.grid,levels=c(-1000,0,500,1000,1500,2000,3000),add=T)
points(tern,pch=19)




#exemple 2
#Orignal Data as per Question
a <- c(0.1, 0.5,0.5, 0.6, 0.2, 0          , 0         , 0.004166667, 0.45) 
b <- c(0.75,0.5,0  , 0.1, 0.2, 0.951612903,0.918103448, 0.7875     , 0.45)
c <- c(0.15,0  ,0.5, 0.3, 0.6, 0.048387097,0.081896552, 0.208333333, 0.10) 
d <- c(500,2324.90,2551.44,1244.50, 551.22,-644.20,-377.17,-100, 2493.04) 
df <- data.frame(a, b, c, d)

#For labelling each point.
df$id <- 1:nrow(df)

#Build Plot
ggtern(data=df,aes(x=c,y=a,z=b),aes(x,y,z)) + 
  stat_density2d(geom="polygon",
                 n=400,
                 aes(fill=..level..,
                     weight=d,
                     alpha=abs(..level..)),
                 binwidth=100) + 
  geom_density2d(aes(weight=d,color=..level..),
                 n=400,
                 binwidth=100) +
  geom_point(aes(fill=d),color="black",size=5,shape=21) + 
  geom_text(aes(label=id),size=3) + 
  scale_fill_gradient(low="yellow",high="red") + 
  scale_color_gradient(low="yellow",high="red") + 
  theme_tern_rgbw() + 
  theme(legend.justification=c(0,1), legend.position=c(0,1)) + 
  guides(fill = guide_colorbar(order=1),
         alpha= guide_legend(order=2),
         color="none") + 
  labs(  title= "Ternary Plot and Filled Contour",
         fill = "Value, V",alpha="|V - 0|")


#exemple3
df <- data.frame(a, b, c, d)
ggtern(df,aes(a,c,b)) + 
  geom_interpolate_tern(aes(value=d,fill=..level..),
                        binwidth=500,
                        colour="white") +
  geom_point(aes(fill=d),color="black",shape=21,size=3) + 
  scale_fill_gradient(low="yellow",high="red") +
  theme(legend.position=c(0,1),legend.justification=c(0,1)) + 
  labs(fill="Value, d")



#exemple 4 (updae ggtern https://stackoverflow.com/questions/36071015/how-to-get-ternary-contour-plots-with-ggtern-2-1-0)ok!!
library(ggtern)

#Orignal Data as per Question
a <- c(0.1, 0.5,0.5, 0.6, 0.2, 0          , 0         , 0.004166667, 0.45) 
b <- c(0.75,0.5,0  , 0.1, 0.2, 0.951612903,0.918103448, 0.7875     , 0.45)
c <- c(0.15,0  ,0.5, 0.3, 0.6, 0.048387097,0.081896552, 0.208333333, 0.10) 
d <- c(500,2324.90,2551.44,1244.50, 551.22,-644.20,-377.17,-100, 2493.04) 
df <- data.frame(a, b, c, d)

#For labelling each point.
df$id <- 1:nrow(df)
ggtern(data=df,aes(x=c,y=a,z=b),aes(x,y,z)) + 
  stat_density_tern(geom="polygon",color='black',
                    n=400,h=0.75,expand = 1.1,
                    base='identity',
                    aes(fill   = ..level..,weight = d),
                    na.rm = TRUE) + 
  geom_point(color="black",size=5,shape=21) +
  geom_text(aes(label=id),size=3) + 
  scale_fill_gradient(low="yellow",high="red") + 
  scale_color_gradient(low="yellow",high="red") + 
  theme_rgbw() + 
  theme(legend.justification=c(0,1), legend.position=c(0,1)) +
  theme_gridsontop() + 
  guides(fill = guide_colorbar(order=1),color="none") + 
  labs(  title= "Ternary Plot and Filled Contour",fill = "Value, V")



OverYvsAll <- function(ls_tabmoys, key, Ymax=300, nom="", optProp="sowing", visuplot=T,...)
{
  #key <- ls_keysc[20]
  #figure de tous les overyielding
  ls_keysc = names(ls_tabmoys)
  
  if (optProp=="sowing")
  { labx <- 'Sowing proportion (Esp. 1)'}
  if (optProp=="actual")
  { labx <- 'Actual proportion (Esp. 1)'}
  
  if (visuplot==T)
  {
    plot(-100, -100, ylim=c(-Ymax,Ymax), xlim=c(0,1), main=nom, xlab=labx, ylab='Overyieding (g.m-2)', ...)
    segments(0, 0, 1, 0, col=1)
  }
  
  resx <- NULL
  resy <- NULL
  
  for (keysc in ls_keysc)
  {
    #keysc <- ls_keysc[3]
    tabmoy <- ls_tabmoys[[keysc]]
    
    #xx <- tabmoy$Semprop1#tabmoy$Yprop1#
    yy <- tabmoy$Ytot
    #actual or sowing proportions?
    if (optProp=="sowing")
    {
      xx <- tabmoy$Semprop1
      labx <- 'Sowing proportion (Esp. 1)'
    }
    if (optProp=="actual")
    {
      xx <- tabmoy$Yprop1
      labx <- 'Actual proportion (Esp. 1)'
    }
    
    lintot <- lsfit(c(xx[1], xx[7]), c(yy[1], yy[7]))
    ylin <- lintot$coefficients[["Intercept"]] + xx*lintot$coefficients[["X"]]
    overY <- yy - ylin
    
    if (keysc != key)
    {
      if (visuplot==T)
      { points(xx, overY, pch=16, col='light grey') }
      resx <- cbind(resx,xx)
      resy <- cbind(resy,overY)
    } else
    {
      savexx <- xx
      saveyy <- overY
    }
  }
  if (visuplot==T)
  { points(savexx, saveyy, pch=16, col='blue', type='b')}
  resx <- cbind(resx,savexx)
  resy <- cbind(resy,saveyy)
  data.frame(x=as.numeric(resx), y=as.numeric(resy))
}






#modif liste fichier
nom<-file.choose()#"C:\\simul\\PMA18\\dtoto3ok.csv"

dtoto <- read.table(nom, header = TRUE,row.names=F, sep=";", dec=".")

dtoto <- read.table("clipboard", header = TRUE, sep="\t", dec=".")




names(dtoto)
dtoto[dtoto$scenario =="55-55" & dtoto$nbplt1==48, ]

#usm 3415 a retirer = ligne211
dtoto <- dtoto[-211,]

tabmoys <- do.call("rbind", ls_tabmoys) #merge a list of data.frames - do.call equalent de map
write.csv(tabmoys, "tabmoys3ok_.csv", row.names=F)

