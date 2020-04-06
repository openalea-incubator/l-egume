

##############
## fonctions criteres stat eval model
##############


rmse <- function(O,P){
  sqrt(mean((O-P)**2,na.rm=T))
}

rmsesCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egs<-0
  resreg<-lm(P~O)
  Preg<-fitted.values(resreg)
  for(k in 1:length(P)){
    error.egs<-(Preg[k]-O[k])^2 # difference entre predit par la regression lineaire et observe
    sum.error.egs<-sum.error.egs+error.egs
  }
  sum.error.egs=unname(sum.error.egs) #permet de rendre le vecteur sans-nom (sans quoi il est avec-nom "1" pour une raison qui m'echappe...)
  return(((1/Nb.O)*sum.error.egs)^0.5)
}


rmseuCoucheney=function(O,P){
  Nb.O=length(O)
  sum.error.egu<-0
  resreg<-lm(P~O)
  Preg<-fitted.values(resreg)
  for(k in 1:length(P)){
    error.egu<-(P[k]-Preg[k])^2 # difference entre predit par la regression lineaire et simule
    sum.error.egu<-sum.error.egu+error.egu
  }
  sum.error.egu=unname(sum.error.egu) #permet de rendre le vecteur sans-nom (sans quoi il est avec-nom "1" pour une raison qui m'echappe...)
  return(((1/Nb.O)*sum.error.egu)^0.5)
}


rrmseCoucheney=function(O,P){
  O.mean=mean(O)
  return(100*rmse(O,P)/O.mean)
}

pRMSEs=function(rmse_, rmses_)
{
  rmses_**2 / rmse_**2
}

pRMSEu=function(rmse_, rmseu_)
{
  rmseu_**2 / rmse_**2
}


efficiencyCoucheney=function(O,P){
  sum.error<-0
  sum.mean.dev<-0
  O.mean=mean(O)
  for(k in 1:length(P)){
    dif<-(O[k]-P[k]) # difference entre valeur simulee et observee
    error<-dif^2
    sum.error<-sum.error+error
    mean.dev<-(O[k]-O.mean)^2 # deviation des observes / moyenne des observes
    sum.mean.dev<-sum.mean.dev+mean.dev
  }
  EF=1-sum.error/sum.mean.dev # efficience du modele ; max 1; the closer to 1 the better
  EF[!is.finite(EF)]<-NA
  return(EF)
}




##############
## fonctions plot melanges binaires
##############


Build_AverageScTable <- function(dtoto, keysc)
{
  sc <- strsplit(keysc," ")[[1]][1]
  mix <- strsplit(keysc," ")[[1]][2]
  mng <- strsplit(keysc," ")[[1]][3]
  sd_ <- strsplit(keysc," ")[[1]][4]
  #recup d'un scenario
  #sc <- '1-1'
  #mix <- 'Fix2-nonFixSimTest'
  #mng <- 'Lusignan30IrrN2'
  
  res <- dtoto[dtoto$scenario==sc & dtoto$mix==mix & dtoto$Mng==mng & dtoto$sd==sd_, ]
  
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


YtotvsProp <- function(tabmoy, Ymax=2200, nom="", optProp="sowing",visuplot=T, visutext=T, labx=NA,col2=2,...)
{
  ## calcul des composante de l'overyielding biomasse et fait un plot (visutext=visualisation des valeurs; visuplot=visulaisation des )
  
  
  #actual or sowing proportions?
  if (optProp=="sowing")
  {
    xx <- tabmoy$Semprop1
    if (is.na(labx))
    {labx <- 'Sowing proportion (Sp. 1)'}
  }
  if (optProp=="actual")
  {
    xx <- tabmoy$Yprop1
    if (is.na(labx))
    {labx <- 'Actual proportion (Sp. 1)'}
  }
  
  #esp pures
  moyesp1_pur <- mean(tabmoy[tabmoy$Semprop1==1., c("Ytot")])
  moyesp2_pur <- mean(tabmoy[tabmoy$Semprop1==0., c("Ytot")])
  
  #calcul des fits des valeurs moyennes
  #modeltot <- smooth.spline(xx, tabmoy$Ytot)
  modeltot <- tryCatch(smooth.spline(xx, tabmoy$Ytot), error=function(e) smooth.spline(xx, tabmoy$Ytot, nknots =5))
  inttot = sum(predict(modeltot, seq(0,1,0.001))$y*0.001) - (moyesp1_pur + moyesp2_pur)/2
  
  #modelesp1 <- smooth.spline(xx, tabmoy$YEsp1)
  modelesp1 <- tryCatch(smooth.spline(xx, tabmoy$YEsp1), error=function(e) smooth.spline(xx, tabmoy$YEsp1, nknots =5))
  intesp1 = sum(predict(modelesp1, seq(0,1,0.001))$y*0.001) - (moyesp1_pur + 0)/2
  
  #modelesp2 <- smooth.spline(xx, tabmoy$YEsp2)
  modelesp2 <- tryCatch(smooth.spline(xx, tabmoy$YEsp2), error=function(e) smooth.spline(xx, tabmoy$YEsp2, nknots =5))
  intesp2 = sum(predict(modelesp2, seq(0,1,0.001))$y*0.001) - (0 + moyesp2_pur)/2
  
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
    #segments(xx[1], tabmoy$Ytot[1], xx[7], tabmoy$Ytot[7], lty=2)
    segments(xx[1], moyesp2_pur, xx[7], moyesp1_pur, lty=2)
    lines(modeltot)
    
    points(xx, tabmoy$YEsp1,col=col2)
    #segments(xx[1], tabmoy$YEsp1[1], xx[7], tabmoy$YEsp1[7], lty=2, col=col2)
    segments(xx[1], 0, xx[7], moyesp1_pur, lty=2, col=col2)
    lines(modelesp1, col=col2)
    
    points(xx, tabmoy$YEsp2,col=4)
    #segments(xx[1], tabmoy$YEsp2[1], xx[7], tabmoy$YEsp2[7], lty=2, col=4)
    segments(xx[1], moyesp2_pur, xx[7], 0, lty=2, col=4)
    lines(modelesp2, col=4)
    
  }
  
  if (visutext==T & visuplot==T)
  {
    text(0.15, 0.97*Ymax, paste('overY: ' ,round(inttot,2)))
    text(0.15, 0.93*Ymax, paste('Sp1: ' , round(intesp1,2)),col=2)
    text(0.15,0.89*Ymax, paste('Sp2: ' ,round(intesp2,2)),col=4)
  }
  
  #renvoie valeurs calculees
  res <- as.list(c(inttot, intesp1, intesp2, propOpt, OverMax, propsowing50, propLegOtp, Ytotmax, propYtotmax))
  names(res) <- c("inttot", "intesp1", "intesp2", "propOpt", "OverMax", "propsowing50", "propLegOtp", "Ytotmax", "propYtotmax")
  res
  
}


QNtotvsProp <- function(tabmoy, Ymax=100, nom="", optProp="sowing", visuplot=T, visutext=T, labx=NA,...)
{
  ## calcul des composante de l'overyielding Ntot et fait un plot (visutext=visualisation des valeurs; visuplot=visulaisation des plots)
  
  
  #actual or sowing proportions?
  if (optProp=="sowing")
  {
    xx <- tabmoy$Semprop1
    if (is.na(labx))
    {labx <- 'Sowing proportion (Sp. 1)'}
  }
  if (optProp=="actual")
  {
    xx <- tabmoy$Yprop1
    if (is.na(labx))
    {labx <- 'Actual proportion (Sp. 1)'}
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



OverYvsAll <- function(ls_tabmoys, key, Ymax=300, nom="", optProp="sowing", visuplot=T,labx=NA,laby=NA,...)
{
  #key <- ls_keysc[20]
  #figure de tous les overyielding
  ls_keysc = names(ls_tabmoys)
  
  if (optProp=="sowing" & is.na(labx))
  { labx <- 'Sowing proportion (Sp. 1)'}
  if (optProp=="actual" & is.na(labx))
  { labx <- 'Actual proportion (Sp. 1)'}
  if (optProp=="sowing" & is.na(laby))
  { laby <- 'Apparent Overyielding (g.m-2)'}
  if (optProp=="actual" & is.na(laby))
  { laby <- 'Overyieding (g.m-2)'}
  
  if (visuplot==T)
  {
    plot(-100, -100, ylim=c(-Ymax,Ymax), xlim=c(0,1), main=nom, xlab=labx, ylab=laby, ...)
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




Which_decile <- function(valparams)
{
  #to find in which decile is a value in a distribution
  qt <- quantile(valparams, probs=seq(0, 1, 0.1))
  qt1 <- as.numeric(valparams<=qt[[2]])*1
  qt2 <- as.numeric(valparams>qt[[2]] & valparams<=qt[[3]])*2
  qt3 <- as.numeric(valparams>qt[[3]] & valparams<=qt[[4]])*3
  qt4 <- as.numeric(valparams>qt[[4]] & valparams<=qt[[5]])*4
  qt5 <- as.numeric(valparams>qt[[5]] & valparams<=qt[[6]])*5
  qt6 <- as.numeric(valparams>qt[[6]] & valparams<=qt[[7]])*6
  qt7 <- as.numeric(valparams>qt[[7]] & valparams<=qt[[8]])*7
  qt8 <- as.numeric(valparams>qt[[8]] & valparams<=qt[[9]])*8
  qt9 <- as.numeric(valparams>qt[[9]] & valparams<=qt[[10]])*9
  qt10 <- as.numeric(valparams>qt[[10]])*10
  qtn <- qt1+qt2+qt3+qt4+qt5+qt6+qt7+qt8+qt9+qt10
  qtn
}


My_AreaPlot <- function(don, titre="", xlab="x", ylab="y", lscol="")
{
  #area plot
  #prends un dataframe don avec x en colonne 1 et les n colonnes de y a mettre en ordre decroissnt (+couleur)
  
  xmin <- min(don[,1])
  xmax <- max(don[,1])
  cumtot <- as.numeric(rowSums(as.matrix(don[,2:dim(don)[2]])))
  ymax <- max(cumtot)
  
  plot(-100,-100, ylim=c(0,ymax), xlim=c(xmin,xmax), xlab=xlab, ylab=ylab, main=titre)
  
  for (i in 2:dim(don)[2])
  {
    #i <- 2
    cumi <- as.numeric(rowSums(as.matrix(don[,i:dim(don)[2]])))
    x <- c(xmin, don[,1], xmax)
    y <- c(0, cumi,0)
    col <- if(lscol != "") lscol[i] else i #genere warnings
    polygon(x,y,col=col)
  }
  
}


PlotDynMStot <- function(MStot, sp_tabSD, sp, lscol="", titre="", ymax=28)
{
  # plot dynamique de MStot au cours du temps par plante avec couleur selon decile
  plot(-10, -10, xlim=c(1,dim(MStot)[1]), ylim=c(0,ymax), main=titre, xlab="t", ylab="MStot")
  for (i in 1:length(sp_tabSD[[sp]]$nump))
  {
    nump <- sp_tabSD[[sp]]$nump[i]
    col <- if(lscol != "") lscol[sp_tabSD[[sp]]$decile[i]] else i #genere warnings
    points(1:dim(MStot)[1], MStot[,nump+1], col=col, type='l')
  }
}


Build_EvolProportions <- function(MStot, sp_tabSD, sp)
{
  #consrtuction d'un tableau res des proportion par decile d'une espece
  
  # 1 MStot esp au cour du temps
  dynMtotsp <- as.numeric(rowSums(as.matrix(MStot[,sp_tabSD[[sp]]$nump+1])))
  res <- data.frame(dynMtotsp, t=1:dim(MStot)[1])
  # 2 ajout des proportion pour chaque decile
  for (dec in 10:1)
  {
    #dec <-9 #numero de decile
    lsp <- sp_tabSD[[sp]][sp_tabSD[[sp]]$decile==dec, c("nump")]
    #lsp+1
    
    frac <- as.numeric(rowSums(as.matrix(MStot[,lsp+1])))*100 / dynMtotsp
    res <- cbind(res, frac)
  }
  names(res) <- c("MStot_esp","t", "dec10", "dec9", "dec8", "dec7", "dec6", "dec5", "dec4", "dec3", "dec2", "dec1")
  
  res
}


