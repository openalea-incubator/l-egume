

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




############## fonction analyse/plot diversite intra


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


PlotDynMStot <- function(MStot, sp_tabSD, sp, lscol="", titre="", ymax=28, append=F)
{
  # plot dynamique de MStot au cours du temps par plante avec couleur selon decile
  if (append==F)
  {
    plot(-10, -10, xlim=c(1,dim(MStot)[1]), ylim=c(0,ymax), main=titre, xlab="t", ylab="MStot")
  }
  for (i in 1:length(sp_tabSD[[sp]]$nump))
  {
    nump <- sp_tabSD[[sp]]$nump[i]
    col <- if(lscol != "") lscol[sp_tabSD[[sp]]$decile[i]] else i #genere warnings
    points(1:dim(MStot)[1], MStot[,nump+1], col=col, type='l')
  }
}


Build_EvolProportions <- function(MStot, sp_tabSD, sp, var="decile")
{
  #consrtuction d'un tableau res des proportion par decile d'une espece
  
  # 1 MStot esp au cour du temps
  dynMtotsp <- as.numeric(rowSums(as.matrix(MStot[,sp_tabSD[[sp]]$nump+1])))
  res <- data.frame(dynMtotsp, t=1:dim(MStot)[1])
  # 2 ajout des proportion pour chaque decile
  for (dec in 10:1)
  {
    #dec <-9 #numero de decile
    lsp <- sp_tabSD[[sp]][sp_tabSD[[sp]][,c(var)]==dec, c("nump")]
    #lsp+1
    
    frac <- as.numeric(rowSums(as.matrix(MStot[,lsp+1])))*100 / dynMtotsp
    res <- cbind(res, frac)
  }
  names(res) <- c("MStot_esp","t", "dec10", "dec9", "dec8", "dec7", "dec6", "dec5", "dec4", "dec3", "dec2", "dec1")
  
  res
}



BuildResDecil <- function(MStot, sp_tabSD)
{
  #fonction pour calculer les decile /sp et a partir des tabSP (avec decile pour 1 parametre donne)
  
  #sp_tabSD <- split(resread[["ls_tabSD"]][[1]], resread[["ls_tabSD"]][[1]]$name)
  #MStot <- resread[["ls_MStot"]][[1]]
  
  sp <- names(sp_tabSD)[1]#"Fix0"
  x1 <- Build_EvolProportions(MStot, sp_tabSD, sp)
  IDlastday <- dim(x1)[1]
  #dec1 <- sum(x1[IDlastday,c(3)])
  dec3_1 <- sum(x1[IDlastday,c(3:5)])
  dec5_1 <- sum(x1[IDlastday,c(3:7)])
  
  sp <- names(sp_tabSD)[2]#"Fix1"
  x2 <- Build_EvolProportions(MStot, sp_tabSD, sp)
  IDlastday <- dim(x2)[1]
  dec3_2 <- sum(x2[IDlastday,c(3:5)])
  dec5_2 <- sum(x2[IDlastday,c(3:7)])
  
  res <- data.frame(dec3_1, dec5_1, dec3_2, dec5_2)
  names(res) <- paste(c("dec3", "dec5", "dec3", "dec5"), param_name, c(names(sp_tabSD)[1], names(sp_tabSD)[1], names(sp_tabSD)[2], names(sp_tabSD)[2]),  sep="_")
  
  res
  #
}


plotMean_Div <-function (matind, xval, title="", xlab="",ylab="",ylim=c(0,100), lscol=NULL)
{
  # plot de la moyenne versus les valeurs par individu
  nbplt <- dim(matind)[2]
  moy_ <-rowMeans(matind, na.rm=T)
  plot(xval, moy_, col="dark grey", lwd=3, main=title, xlab=xlab, ylab=ylab,ylim=ylim)
  for (idp in 1:nbplt)
  {
    if (is.null(lscol))
    {col<-idp}
    else
    {col <- lscol[idp]}
    
    points(xval, matind[,idp], col=idp, type="l")
  }
  
}




########### functions competition coefficient indices
# Sackeville Hamilton 2001

Yresp_densite1 <- function(a ,b, densite)
{
  # eq 2.3 - sackeville hamilton exprimee en reponse au Ytot
  Ytot = densite/(a+b*densite)
  Ytot
}


Calc_Beta_coeff_Jul <- function(x)
{
  #x = tableau dtoto des culture pure avec Ytot, densite, nbplt et surfsolref
  # eq 2.3 - sackeville hamilton
  #par fit lineaire sur l'inverse plant
  MY_plant1 <- x$Ytot/x$nbplt*x$surfsolref
  mod1 <- lm(1/MY_plant1 ~ x$densite)
  a1 <- as.data.frame(summary(mod1)[["coefficients"]])$Estimate[1] #intercept
  b1 <- as.data.frame(summary(mod1)[["coefficients"]])$Estimate[2] #slope
  beta1 <- b1/a1
  res <- list(mod1, a1,b1,beta1)
  names(res) <- c("model", "a","b","beta")
  res
} 
#inconvenient: sensible aux point pHD, parfois valeur negative de beta
#avantage: plan simple en 2 points de densite utilisable (isole/dense pur)
#Calc_Beta_coeff_Jul(pur1)
#Calc_Beta_coeff_Jul(pur2)


Calc_Beta_coeff <- function(x)
{
  #x = tableau dtoto des culture pure avec Ytot, densite, nbplt et surfsolref
  # eq 2.3 - sackeville hamilton sous forme non lineaire
  #par fit non lineaire
  
  startlist <- list(a=0.01, b=0.0004)
  model1 <- nls(Ytot~Yresp_densite1(a ,b, densite), data=x, start=startlist )
  parameters1 <- summary(model1)[["parameters"]]
  
  a1 <- parameters1[1] #intercept
  b1 <- parameters1[2] #slope
  beta1 <- b1/a1
  res <- list(model1, a1,b1,beta1)
  names(res) <- c("model", "a","b","beta")
  res
}


#Calc_Beta_coeff(pur1)
#Calc_Beta_coeff(pur2)


Yresp_densite2 <- function(a ,beta, gamma, densite1, densite2)
{
  # eq 2.4 - sackeville hamilton exprimee en reponse au Yespi a densite de deux espece
  #par modele lineaire sur inverse rendement
  inv_Yi = a + a*beta*densite1 + a*gamma*densite2
  inv_Yi
}

#Yresp_densite2_diag <- function(a ,beta, gamma, dtot, densite1)
#{
#  # eq 2.4 - sackeville hamilton exprimee en reponse au Yespi a densite de deux espece
#  # une seule variable dans dipositif de DeWit car dtot=cst (a forcer dans l'optimisation)
#  inv_Yi = a + a*beta*densite1 + a*gamma*(dtot-densite1)
#  inv_Yi
#}
#marche pas plus... seule sur diag

Yresp_densite2bis <- function(a ,beta, gamma, densite1, densite2)
{
  # eq 2.4 - sackeville hamilton exprimee en reponse au Yespi a densite de deux espece
  #par modele non lineaire
  Yesp1 = densite1 / (a + a*beta*densite1 + a*gamma*densite2)
  Yesp1
}
#Yresp_densite2bis(a ,beta, gamma, densite1, densite2)


Calc_Gamma_coeffesp12 <- function(x, iso1, iso2, res1, res2, free=F)
{
  # x = tableau dtoto des culture pure et associee avec Ytot, densite1, densite2, Yesp1, Yesp2, nbplt et surfsolref
  # iso1 et iso2: ligne equivalente avec valeur des plantes isolee
  #res1 et res2: resultats des fits des reponses en pur des especes 
  # !!: actuellement: force a et beta values lors du fit de Yi~Yresp_densite2
  # fit eq 2.4 
  #free = T -> laisse fitter les 3 params
  #avec reponse non lineaire
  
  #Esp1
  df <- x[,c("YEsp1", "densite1", "densite2")]
  df$inv_Yi <- x$densite1/x$YEsp1
  df <- df[df$densite1!=0,] #retire les purs 
  
  #ajout point isole 50/50
  iso50 <- data.frame(YEsp1 = iso1$YEsp1/(iso1$YEsp1+iso2$YEsp2)*iso1$YEsp1, densite1=iso1$densite1/2, densite2=iso1$densite1/2)
  iso50$inv_Yi <- iso50$densite1/iso50$YEsp1
  df <- rbind(df,iso50)
  
  #avec nls en forcant beta et a!
  startlist <- list(a=res1[["a"]], beta=res1[["beta"]], gamma=0.002)
  if (free==F)
  {
    #force a et beta
    minis <- c(res1[["a"]],res1[["beta"]],-1.)
    maxis <- c(res1[["a"]],res1[["beta"]],1.)
  } else
  {
    #fit libre borne
    #minis <- c(0.,res1[["beta"]],0.)
    #maxis <- c(1.,res1[["beta"]],1.)
    minis <- c(0.,0.,-1.)
    maxis <- c(10.,1.,1.)
  }
  
  #model1 <- nls(inv_Yi~Yresp_densite2(a ,beta, gamma, densite1, densite2), data=df, start=startlist ,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
  model1 <- nls(YEsp1~Yresp_densite2bis(a ,beta, gamma, densite1, densite2), data=df, start=startlist ,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
  parameters1 <- summary(model1)[["parameters"]]
  #yes! rq: forcage de a change pas grand chose
  
  #Esp2
  df <- x[,c("YEsp2", "densite1", "densite2")]
  df$inv_Yi <- x$densite2/x$YEsp2
  df <- df[df$densite2!=0,] #retire les purs 
  
  #ajout point isole 50/50
  iso50 <- data.frame(YEsp2 = iso2$YEsp2/(iso1$YEsp1+iso2$YEsp2)*iso2$YEsp2, densite1=iso1$densite1/2, densite2=iso1$densite1/2)
  iso50$inv_Yi <- iso50$densite1/iso50$YEsp2
  df <- rbind(df,iso50)
  
  #inverse les noms pour le fit
  names(df) <- c("YEsp1", "densite2", "densite1", "inv_Yi")
  
  startlist <- list(a=res2[["a"]], beta=res2[["beta"]], gamma=0.002)
  if (free==F)
  {
    minis <- c(res2[["a"]],res2[["beta"]],-1.)
    maxis <- c(res2[["a"]],res2[["beta"]],1.)
  } else
  {
    #minis <- c(0.,res2[["beta"]],0.)
    #maxis <- c(1.,res2[["beta"]],1.)
    minis <- c(0.,0.,-1.)
    maxis <- c(10.,1.,1.)
  }
  
  #model2 <- nls(inv_Yi~Yresp_densite2(a ,beta, gamma, densite1, densite2), data=df, start=startlist ,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
  model2 <- nls(YEsp1~Yresp_densite2bis(a ,beta, gamma, densite1, densite2), data=df, start=startlist ,trace=TRUE,algorithm="port",lower=minis,upper=maxis)
  parameters2 <- summary(model2)[["parameters"]]
  
  res <- list(model1, model2, parameters1, parameters2)
  names(res) <- c("model1", "model2", "parameters1", "parameters2")
  res
}
#A faire: rendre facultatif les iso1 et 2!
#rendre facultatif res11 et res2
#x <- diag 
#resg <- Calc_Gamma_coeffesp12(x, iso1, iso2, res1, res2)
#Calc_Gamma_coeffesp12(all, iso1, iso2, res1, res2)
#Calc_Gamma_coeffesp12(tribas, iso1, iso2, res1, res2)



Calc_Sij_coefficients <- function(res1, res2, parameters1, parameters2)
{
  # substitution rates Sij Eq 2.7
  #esp1
  S12 <- parameters1[3] / res1[["beta"]] 
  #esp2
  S21 <- parameters2[3] / res2[["beta"]] 
  res <- list(S12, S21)
  names(res) <- c("S12", "S21")
  res
}
#Calc_Sij_coefficients(res1, res2, parameters1, parameters2)


Calc_Eij_coefficients <- function(res1, res2, parameters1, parameters2)
{
  # relative competitive effect  Eij Eq 2.8
  #esp1
  E12 <- 1/(parameters1[3] / res1[["beta"]])
  #esp2
  E21 <- 1/(parameters2[3] / res2[["beta"]])
  res <- list(E12, E21)
  names(res) <- c("E12", "E21")
  res
}
#Calc_Eij_coefficients(res1, res2, parameters1, parameters2)


Calc_Rij_coefficients <- function(res1, res2, parameters1, parameters2)
{
  # relative competitive effect  Eij Eq 2.8
  #esp1
  S12 <- parameters1[3] / res1[["beta"]] 
  #esp2
  S21 <- parameters2[3] / res2[["beta"]]
  res <- 1 / (S12*S21) #le meme pour les deux
  names(res) <- c("R")
  res
}
#Calc_Rij_coefficients(res1, res2, parameters1, parameters2)
# R bien superieur a zero!


Calc_CompetitionIntensity<- function(YEsp1, densite1, iso_ind)
{
  #calculate competition intensity relative to isolated plant of sp (with vectors Yesp and density)
  #eq 2.5 
  Yi <- YEsp1/densite1
  #iso_ind <- iso$YEsp1/iso$densite
  res <- log10(iso_ind/Yi)
  res
}

#!! prendre en charge le calcul des moyennes d'isole
#Calc_CompetitionIntensity(x$YEsp1, x$densite1, iso1$YEsp1/iso1$densite1)
#Calc_CompetitionIntensity(x$YEsp2, x$densite2, iso2$YEsp2/iso2$densite2)

Calc_RCI_coeff <- function(YEsp1, densite1, iso_ind)
{
  # Relative compeition index (Maamouri et al2017)
  Yi <- YEsp1/densite1
  #iso_ind <- iso$YEsp1/iso$densite
  res <- (iso_ind-Yi)/iso_ind
  res
}
#Calc_RCI_coeff(x$YEsp1, x$densite1, iso1$YEsp1/iso1$densite1)
#Calc_RCI_coeff(x$YEsp2, x$densite2, iso2$YEsp2/iso2$densite2)


Calc_CEi <- function(Nb,M,deltaRY)
{
  #complementarity effect (Hector et Loreau, 2001)
  #Nb: nb sp du melange
  #M: vecteur des rendement en pur par sp du melange
  # deltaRY: vecteur delta de RY par sp du melange
  CEi <- Nb*mean(M)*mean(deltaRY)
  CEi
}

Calc_SEi <- function(Nb,M,deltaRY)
{
  #selection effect (Hector et Loreau, 2001)
  #Nb: nb sp du melange
  #M: vecteur des rendement en pur par sp du melange
  # deltaRY: vecteur delta de RY par sp du melange
  SEi <- Nb*cov(deltaRY, M)
  SEi
}


Calc_CESE_diag <- function(diag)
{
  #calculate CE and SE coeffcicient (Hector et Loreau 2001) 
  # for a tabmoy or diag of a scenario (1 seed), with minimum pure controls ans 1 mixture
  
  x <- diag
  
  x$M1 <- x[x$Semprop1==1.,c("YEsp1")] #yield esp pur1 (meme densite) 
  x$M2 <- x[x$Semprop1==0.,c("YEsp2")] #yield esp pur2 (meme densite) 
  
  #x$deltaRY1 <- x$Yprop1 - x$Semprop1
  #x$deltaRY2 <- x$Yprop2 - (1-x$Semprop1)
  x$deltaRY1 <- x$YEsp1/x$M1-x$Semprop1
  x$deltaRY2 <- x$YEsp2/x$M2 - (1-x$Semprop1)
  
  x$Ytheo1 <- x$M1*x$Semprop1
  x$Ytheo2 <- x$M2*(1-x$Semprop1)
  x$Yteo <- x$Ytheo1  + x$Ytheo2
  x$Ytot - x$Yteo
  x$OYa <- x$Ytot - x$Yteo#x$YEsp1 - x$Semprop1*x$M1 + x$YEsp2 - (1-x$Semprop1)*x$M2
  
  #somme des RYi*Mi obs - somme des RYi/Mi theo
  #(x$M1 * x$YEsp1/x$M1 + x$M2 * x$YEsp2/x$M2) - (x$M1 * x$Semprop1 + x$M2 *(1-x$Semprop1))
  # bon!
  #somme avec les deltaRY
  #x$M1 * (x$YEsp1/x$M1-x$Semprop1) + x$M2*(x$YEsp2/x$M2 - (1-x$Semprop1))
  #x$M1 * x$deltaRY1 + x$M2*x$deltaRY2
  # bon!
  
  #selection des liste de vecteur des asso (retire purs)
  ls_vdeltaRY <- x[, c("deltaRY1", "deltaRY2")]
  ls_vM <- x[, c("M1", "M2")]
  
  #calcul par couvert
  CE <- NULL
  SE <- NULL
  for (i in 1:length(ls_vdeltaRY[,1]))
  {
    #i <-2 #numero ligne
    Nb <- length(ls_vdeltaRY[i,])
    M <- as.numeric(ls_vM[i,])
    deltaRY <- as.numeric(ls_vdeltaRY[i,])
    
    CEi <- Calc_CEi(Nb,M,deltaRY) # Nb*mean(M)*mean(deltaRY)
    SEi <- Calc_SEi(Nb,M,deltaRY)# Nb*cov(deltaRY, M)
    CE <- rbind(CE, CEi)#rbind(CE, CEi/Nb)
    SE <- rbind(SE, SEi)#rbind(SE, SEi/Nb)
    #facteur 2 (Nb qui traine) -> #c'est reference qu'est somme pas moy des cultures pures ! 
    #?OY = CE+0.5*SE
    
  }
  
  x$CE <- CE
  x$SE <- SE
  x
}




#### fonctions de plot

Plt_Yresp_densite1 <- function(x, res, titre="")
{
  #plot de la reponse a la densite d'une espece pure avec le coeff beta
  plot(x$densite, x$Ytot, ylim=c(0, h= 1/res[["b"]]), main=titre)
  parameters1 <- summary(res[["model"]])[["parameters"]]
  vals <- Yresp_densite1(a = parameters1[1], b=parameters1[2], densite=seq(0,400,1.))
  lines(seq(0,400,1.), vals, type='l',col=2)
  abline(h= 1/res[["b"]],col=3, lty=2)
  abline(h= 0.5/res[["b"]],col=3, lty=2)
  abline(v= 1/res[["beta"]],col=4, lty=3)
  text(1/res[["beta"]]+30, 50, "1/beta", col=4)
  text(1/res[["beta"]]+30, 0.5/res[["b"]]-50, "b/2", col=3)
  text(200,50,paste ( "beta = b/a = ", round(res[["beta"]], 4)))
  text(200, 200, paste ( "1/a = ", round(1/res[["a"]], 4)))
}



Plot_diag_respFitsd1d2 <- function(x, parameters1, parameters2, dmax=400., titre="")
{
  # Plot des ajustement de gamma sur une diagonale de dispositif de deWit (substitution)
  # x=data.frame de donnee, parameters1 et 2: fits des 2 especes
  
  plot(x$densite1, x$Ytot, ylim=c(0,2300), main=titre)
  points(x$densite1, x$YEsp1, col=2)
  points(x$densite1, x$YEsp2, col=4)
  
  vals <- Yresp_densite2(a = parameters1[1], beta=parameters1[2], gamma=parameters1[3], densite1=seq(0,dmax,1.), densite2=rev(seq(0,dmax,1.)))
  Yesp1 <- (1/vals)*seq(0,dmax,1.)
  lines(seq(0,dmax,1.), Yesp1, type='l',col=2)
  #OK!!
  
  vals2 <- Yresp_densite2(a = parameters2[1], beta=parameters2[2], gamma=parameters2[3], densite1=rev(seq(0,dmax,1.)), densite2=seq(0,dmax,1.))
  Yesp2 <- (1/vals2)*rev(seq(0,dmax,1.))
  lines(seq(0,dmax,1.), Yesp2, type='l',col=4)
  
  lines(seq(0,dmax,1.), Yesp1+Yesp2, type='l')
  text(220,1000, paste("gamma1: ",round(parameters1[3], 5)), col=2)
  text(220,900, paste("gamma2: ",round(parameters2[3], 5)), col=4)
  
}
#Plot_diag_respFitsd1d2(x, parameters1, parameters2, dmax=400.)


Plot_OneOne_Resp_dtot <- function(x, parameters1, parameters2, dmax=400.,titre="")
{
  # Plot des ajustement de gamma sur la diagonale 1:1 (50/50 semis) (additif) - pour Ytot
  # x=data.frame de donnee, parameters1 et 2: fits des 2 especes
  
  plot(x$densite1+x$densite2, x$Ytot, ylim=c(0,2300), main=titre, xlim=c(0, max(x$densite1+x$densite2)))
  
  vals <- Yresp_densite2(a = parameters1[1], beta=parameters1[2], gamma=parameters1[3], densite1=seq(0,dmax,1.), densite2=seq(0,dmax,1.))
  Yesp1 <- (1/vals)*seq(0,dmax,1.)
  #lines(seq(0,400,1.)*2, Yesp1, type='l',col=2)
  
  vals2 <- Yresp_densite2(a = parameters2[1], beta=parameters2[2], gamma=parameters2[3], densite1=seq(0,dmax,1.), densite2=seq(0,dmax,1.))
  Yesp2 <- (1/vals2)*seq(0,dmax,1.)
  
  Ytot <- Yesp1+Yesp2
  lines(seq(0,dmax,1.)*2, Ytot, type='l')
}
#Plot_OneOne_Resp_dtot(x, parameters1, parameters2, dmax=400.)

#lines(seq(0,400,1.)*2, Ytot-Yesp1, type='l',col=4)
#graph pas bon?? -> si 


Plot_OneOne_prop <- function(x, parameters1, parameters2, dmax=400., titre="")
{
  # Plot des ajustement de gamma sur la diagonale 1:1 (50/50 semis) (additif) - pour proportions d'especes
  # x=data.frame de donnee, parameters1 et 2: fits des 2 especes
  
  plot(x$densite1+x$densite2, x$YEsp1/x$Ytot, ylim=c(0,1), main=titre, xlim=c(0, max(x$densite1+x$densite2)), xlab="",ylab="sp prop", col=2, cex.axis=0.8)
  points(x$densite1+x$densite2, x$YEsp2/x$Ytot, col=4)
  
  vals <- Yresp_densite2(a = parameters1[1], beta=parameters1[2], gamma=parameters1[3], densite1=seq(0,dmax,1.), densite2=seq(0,dmax,1.))
  Yesp1 <- (1/vals)*seq(0,dmax,1.)
  #lines(seq(0,400,1.)*2, Yesp1, type='l',col=2)
  
  vals2 <- Yresp_densite2(a = parameters2[1], beta=parameters2[2], gamma=parameters2[3], densite1=seq(0,dmax,1.), densite2=seq(0,dmax,1.))
  Yesp2 <- (1/vals2)*seq(0,dmax,1.)
  Ytot <- Yesp1+Yesp2
  
  points(seq(0,dmax,1.)*2, Yesp1/Ytot, col=2, type='l')
  points(seq(0,dmax,1.)*2, Yesp2/Ytot, col=4, type='l')
  text(500,0.6, paste("gamma1: ",round(parameters1[3], 5)), col=2, cex=0.75)
  text(500,0.4, paste("gamma2: ",round(parameters2[3], 5)), col=4, cex=0.75)
  
}
#Plot_OneOne_prop(x, parameters1, parameters2)




Plot_resCE_SE <- function(resCE_SE, titre="")
{
  #faire le plot - CE-SE (indices Loreau) pour une diag
  #resCE_SE <- data.frame(Semprop1=x$Semprop1, OYa=x$OYa, CE, SE)
  
  plot(resCE_SE$Semprop1, resCE_SE$OYa, ylim=c(-400,400), type='b', ylab="OYa", xlab="Prop1", main=titre)
  #ajouter les surfaces! (un peu transparentes?)
  
  x <- c(0., resCE_SE$Semprop1, 1.)
  y <- c(0, resCE_SE$SE*0.5,0)#+resCE_SE$SE
  col <- rgb(0,0,1,1/4)#4
  polygon(x,y,col=col)
  x <- c(0., resCE_SE$Semprop1, 1.)
  y <- c(0, resCE_SE$CE,0)
  col <- rgb(1,0,0,1/4)#2
  polygon(x,y,col=col)
  points(resCE_SE$Semprop1, resCE_SE$OYa)
  segments(0,0, 1,0, lty=2)
  
  text(0.5,-200,paste("mean CE: ", round(mean(resCE_SE$CE),1)), col=2)
  text(0.5,-250,paste("mean SE: ", round(0.5*mean(resCE_SE$SE),1)), col=4)
  text(0.5,-300,paste("mean aOY: ", round(mean(resCE_SE$OYa),1)))
  #pourrait calculer les max!
  
}
#Plot_resCE_SE(resCE_SE, titre="")




#fonction des exemple de pairs
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
#df <- data.frame(retard,Val_param,ParaMvois,PARivois,MScumvois, MStot_ini, MStot_fin, MStot_coupe1, MStot_coupe2, MStot_coupe3, MStot_coupe4, MStot_coupe5)
#pairs(df, lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=key)



plot_ranking_lines <- function(tab, id_refdec=1, decHaut=9, decBas=2, id_refline=NA, titre='', ymax=1, xlab='',ylab='')
{
  # plot de trajectoires d'interaction/ranking entre modalite classees en colonnes dans un tableau (tab)
  # visualise selection de decile ou d'individus particuliers dans une population
  
  #id_refdec <- 3#1 #colone de ref pour les deciles
  #id_refline <- NA#10#
  #decHaut <- 9
  #decBas <- 2
  #titre <- paste(trait, sp)
  #ymax <- 15000
  
  nb_dates <- dim(tab)[2]
  plot(-10,-10, xlim=c(0,nb_dates+1), ylim=c(0, ymax),main=titre,xlab=xlab,ylab=ylab)
  for (i in 1:dim(dd)[1])
  {
    points(1:nb_dates, as.numeric(tab[i,]), type="b", col="grey")
  }
  
  #ajout decile de la date 'id_refdec'
  Dec_ <- quantile(tab[,id_refdec], na.rm=T, probs = seq(0, 1, 0.1))
  D9 <- as.numeric(Dec_[decHaut+1]) #decile 1 haut
  D2 <- as.numeric(Dec_[decBas+1]) #decile 2 bas
  for (i in 1:dim(tab)[1])
  {
    if (! is.na(tab[i,id_refdec]))
    {
      
      if (tab[i,id_refdec] > D9)
      {
        points(1:nb_dates, as.numeric(tab[i,]), type="b", col="red")
      }
      
      if (tab[i,id_refdec] < D2)
      {
        points(1:nb_dates, as.numeric(tab[i,]), type="b", col="blue")
      }
    }
  }
  
  #ajout de la ligne a surligner
  if (! is.na(id_refline))
  {
    points(1:nb_dates, as.numeric(tab[id_refline,]), type="b", col=1, lwd=2)
  }
}

#decile sur base 1ere date
#plot_ranking_lines(dd, id_refdec=1, titre=paste(trait, sp), ymax=15000)
#decile sur base derniere date
#plot_ranking_lines(dd, id_refdec=3, titre=paste(trait, sp), ymax=15000)
#visu plante 50 et pas les deciles
#nbp <- 50
#plot_ranking_lines(dd, id_refdec=1, decHaut=10, decBas=0, id_refline=nbp, titre=paste(trait, sp, nbp), ymax=15000)
  




#fonction pour determiner l'id des voisins d'ordre 1
ls_idvois_ordre1 <- function(n, cote, nblignes)
{
  # pour une plante n, dans un dispocitif regulier arrange en colonnes croissantes de cote indiv
  nbindiv <- cote * nblignes
  ls_defaut <- c(n - (cote + 1), n - cote, n - (cote - 1), n - 1, n + 1, n + (cote - 1), n + cote, n + (cote + 1))
  
  if (n %% cote == 0)  # bord haut
  {
    ls_defaut[1] <- ls_defaut[1] + cote
    ls_defaut[4] <- ls_defaut[4] + cote
    ls_defaut[6] <- ls_defaut[6] + cote
  }
  
  if ((n + 1) %% cote == 0)  # bord bas
  {
    ls_defaut[3] <- ls_defaut[3] - cote
    ls_defaut[5] <- ls_defaut[5] - cote
    ls_defaut[8] <- ls_defaut[8] - cote
  }
   
  
  for (i in 1:length(ls_defaut))
  {
    if (ls_defaut[i] < 0)  # bord gauche
    {ls_defaut[i] <- ls_defaut[i] + nbindiv}
  
    if (ls_defaut[i] >= nbindiv)  # bord droit
    {ls_defaut[i] <- ls_defaut[i] - nbindiv}
  }
  
  ls_defaut
}
#ls_idvois_ordre1(9,6,4)#marche pas pour zero
#pour ordre 2: voisin d'ordre 1 de tous tes voisins!
#!! prevu pour id python commencant a zero




calc_norm_par <- function(tabpar,lspar, plot_=F, main_="")
{
  # fonction pour Calcul des valeur normalisee (par la moyenne) des parametresSD et la moyenne des valeurs normalisee
  # avec plot_ a True et les coord x,y,  fait un graph de visu
  
  nbpar <- length(lspar)
  ls_resNorm <- vector("list", length=(nbpar+1))
  names(ls_resNorm) <- c(lspar, "mean_norm_par")
  
  ls_col_ <- 1:nbpar #a passer en argument eventuellement
  
  
  Val_par <- tabpar[,c(lspar[1])]
  #normalise par la moyenne (de ce qui est donne en entree: communaute ou population)
  norm_par <- Val_par/mean(Val_par) 
  ls_resNorm[[lspar[1]]] <- norm_par
  mean_norm_par <- norm_par
  
  if (plot_ == T)
  {plot(tabpar$x, tabpar$y, cex=1.5*norm_par,col="blue", main=main_,xlab="",ylab="")}
  
  for (i in 2:nbpar)
  {
    Val_par <- tabpar[,c(lspar[i])]
    norm_par <- Val_par/mean(Val_par)
    ls_resNorm[[lspar[i]]] <- norm_par
    mean_norm_par <- mean_norm_par+norm_par
    
    if (plot_ == T)
    {points(tabpar$x, tabpar$y, cex=1.5*norm_par,col=ls_col_[i])}
    
  }
  ls_resNorm[["mean_norm_par"]] <- mean_norm_par/nbpar
  names(ls_resNorm)[1:nbpar] <- paste(lspar,"Norm", sep="")
  ls_resNorm <- as.data.frame(ls_resNorm)
  ls_resNorm
}

#ls_resNorm <- calc_norm_par(x ,lspar, plot_=F)
#ParamNorm <- calc_norm_par(temptab[,lspar] ,lspar, plot_=F)$mean_norm_par



def_indice_vois5050 <- function(cote, nblignes)
{
  #id des voisins pour un damier a 50/50 (1 sur 2 est un Kin/nonKin)
  
  #cote <- 16
  #nblignes <- 16
  
  ls_idvois <- vector("list", length=(cote*nblignes))
  names(ls_idvois) <- 1:(cote*nblignes)
  ls_idKin <- vector("list", length=(cote*nblignes))
  names(ls_idKin) <- 1:(cote*nblignes)
  ls_idnonKin <- vector("list", length=(cote*nblignes))
  names(ls_idnonKin) <- 1:(cote*nblignes)
  
  for (i in 1:(cote*nblignes))
  {
    idvois <- ls_idvois_ordre1(i-1, cote, nblignes) +1 # appel avec i-1 (pour comme nump python) # ajout 1 a sortie pour rang R
    ls_idvois[[i]] <- idvois 
    ls_idKin[[i]] <- idvois[c(1,3,6,8)]
    ls_idnonKin[[i]] <- idvois[c(2,4,5,7)]
  }
  
  list(ls_idvois, ls_idKin, ls_idnonKin)
}
#ls_idv <- def_indice_vois5050(cote=16, nblignes=16)
#ls_idvois <- ls_idv[[1]]
#ls_idKin <- ls_idv[[2]]
#ls_idnonKin <- ls_idv[[3]]

#a generaliser pour autres configuations (ou a lire qqs part!)




calc_neighb_param <- function(tabpar,lspar, ls_idvois, ls_idKin, ls_idnonKin)
{
  #calculate average parameter value of order 1 neighbours, with kin and non kin in a binary mixture
  #can be used for any vector and any list of lspar (not only parameters)
  
  nbpar <- length(lspar)
  ls_res <- vector("list", length=nbpar)
  names(ls_res) <- c(lspar)
  
  for (param_name in lspar)
  {
    
    #param_name <- lspar[1]
    Val_param <- tabpar[,c(param_name)]
    
    
    ParaMvois <- NULL
    ParaMKin <- NULL
    ParaMnonKin <- NULL
    for (i in 1:(cote*nblignes))
    {
      #ALL ordre 1
      ParaMvois <- cbind(ParaMvois, mean(Val_param[ls_idvois[[i]]]))
      #Kin/NonKin
      ParaMKin <- cbind(ParaMKin, mean(Val_param[ls_idKin[[i]]]) )
      ParaMnonKin <- cbind(ParaMnonKin, mean(Val_param[ls_idnonKin[[i]]]) )
    }
    ParaMvois <- as.numeric(ParaMvois)
    ParaMKin <- as.numeric(ParaMKin)
    ParaMnonKin <- as.numeric(ParaMnonKin)
    
    res <- data.frame(ParaMvois, ParaMKin, ParaMnonKin)
    names(res) <- paste(param_name, c("Mvois", "MKin", "MnonKin"), sep="")
    
    
    ls_res[[param_name]] <- res
    
  }
  
  res <- as.data.frame(ls_res)
  names(res) <- as.character(as.data.frame(t(as.data.frame(strsplit(names(res),"\\."))))$V2)
  
  res
}
#pourrait prevoir de mettre ls_idvois, ls_idKin, ls_idnonKin dans la table d'entree au prelablable
#calc_neighb_param(x,lspar, ls_idvois, ls_idKin, ls_idnonKin)



Calc_MSindiv_Corr <- function(ltoto, ls_toto_paquet, ls_paramSD, lspar=c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh"))
{
  ## fonction pour mettre en forme valeurs de parametres normalise, valeur d'effet de voisinnage, et calculer les correlations entre indices
  
  #key <- names(sp_dtoto)[260]#[330]#[16]#[3]#[31]#[19]#
  #ls_toto_paquet <- sp_dtoto[[key]]$name
  #ltoto <- read_ltoto(ls_toto_paquet)
  #names(ltoto[[ls_toto_paquet]])
  #ltoto[[ls_toto_paquet]]$V1
  dat <- ltoto[[ls_toto_paquet]]
  nb <- dim(dat)[2]-2
  
  #lspar <-  c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")
  param_name <- "phyllochron"#"Len"#ls_par[1] #pour exemple, pas utilise
  resread <- read_lsSD_MStot(ltoto, ls_paramSD, param_name)
  sp_tabSD <- split(resread[["ls_tabSD"]][[1]], resread[["ls_tabSD"]][[1]]$name)
  MStot <- resread[["ls_MStot"]][[1]]
  #res <- BuildResDecil(MStot, sp_tabSD)
  
  #c(187,229,282,334) #dates de coupes fixes
  MStot_ini <- as.numeric(MStot[60,])#30
  MStot_coupe1 <- as.numeric(MStot[127,])#65
  MStot_coupe2 <- as.numeric(MStot[169,])#100
  MStot_coupe3 <- as.numeric(MStot[222,])#150
  MStot_fin <- as.numeric(MStot[dim(MStot)[1],])
  #MStot_coupe4 <- as.numeric(MStot[200,])
  #MStot_coupe5 <- as.numeric(MStot[250,])
  #hist(MStot_fin, main=key)
  #hist(MStot_coupe1, main=key)
  
  #temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","x","y","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  
  
  #ordonne dans l'ordre des nump!!
  temptab <- temptab[order(temptab$nump),]
  
  
  #Val_param <- temptab[,c(param_name)]#temptab$phyllochron
  #Val_param <- temptab$Len
  #hist(Val_param, main=key)
  
  #calcul de la valeur normalisee des parametres (multi-trait)
  temptab$phyllochron[temptab$phyllochron<8] <- 8 #pour les valeur <0 mise a 10-10!
  temptab$phyllochron <- 1/(temptab$phyllochron)
  temptab$PPtreshh <- 24-temptab$PPtreshh
  ParamAllNorm <- calc_norm_par(temptab[,lspar] ,lspar, plot_=F)$mean_norm_par
  temptab$ParamAllNorm <- ParamAllNorm 
  
  #agrege par Light / N (specifique papier beatrice)
  lightPar <- c("Len","Lfeuille","phyllochron")
  ParamLightNorm <- calc_norm_par(temptab[,lightPar] ,lightPar, plot_=F)$mean_norm_par
  temptab$ParamLightNorm <- ParamLightNorm
  NPar <- c("Vmax2", "ELmax", "PPtreshh")
  ParamNNorm <- calc_norm_par(temptab[,NPar] ,NPar, plot_=F)$mean_norm_par
  temptab$ParamNNorm <- ParamNNorm
  
  
  
  #calcul des moyenne des voisins
  x <- temptab[,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm")]
  #transforme param phyllochrone et PPtreshh pour avoir effet positif pour valeur croissante
  
  
  resN <- calc_neighb_param(x,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm"), ls_idvois, ls_idKin, ls_idnonKin)
  temptab <- cbind(temptab, resN)
  #caluler les difference pour sp1 et sp2
  temptab$diffMvoisNorm <- temptab$ParamAllNormMvois - temptab$ParamAllNorm
  temptab$diffMKinNorm <- temptab$ParamAllNormMKin - temptab$ParamAllNorm
  temptab$diffMnonKinNorm <- temptab$ParamAllNormMnonKin - temptab$ParamAllNorm
  temptab$diffMvoisLightNorm <- temptab$ParamLightNormMvois - temptab$ParamLightNorm
  temptab$diffMvoisNNorm <- temptab$ParamNNormMvois - temptab$ParamNNorm
  
  
  #recup PARiPlante et N uptake plante et faire cumul
  PARi <- dat[dat$V1=='PARiPlante',3:(3+nb-1)] #
  for (i in 1:nb) {PARi[,i] <- cumsum(PARi[,i])}
  Nuptake <- dat[dat$V1=='Nuptake_sol',3:(3+nb-1)] #sans fixation!!!
  for (i in 1:nb) {Nuptake[,i] <- cumsum(Nuptake[,i])}
  PARi_fin <- as.numeric(PARi[dim(PARi)[1],])
  Nuptake_fin <- as.numeric(Nuptake[dim(Nuptake)[1],])
  
  
  #calcul du cumul de biomasse, note moyenne, uptake des voisins
  MScumvois <- NULL
  MScumKin <- NULL
  MScumnonKin <- NULL
  PARivois <- NULL
  PARiKin <- NULL
  PARinonKin <- NULL
  Nuptakevois <- NULL
  NuptakeKin <- NULL
  NuptakenonKin <- NULL
  for (i in 1:(cote*nblignes))
  {
    #ALL ordre 1
    MSvois <- sum(MStot_fin[ls_idvois[[i]]])
    MScumvois <- cbind(MScumvois, MSvois)
    PARivois <- cbind(PARivois, sum(PARi_fin[ls_idvois[[i]]]) )
    Nuptakevois <- cbind(Nuptakevois, sum(Nuptake_fin[ls_idvois[[i]]]) )
    
    #Kin/NonKin
    MScumKin <- cbind(MScumKin, sum(MStot_fin[ls_idKin[[i]]]))
    MScumnonKin <- cbind(MScumnonKin, sum(MStot_fin[ls_idnonKin[[i]]]))
    PARiKin <- cbind(PARiKin, sum(PARi_fin[ls_idKin[[i]]]) )
    NuptakeKin <- cbind(NuptakeKin, sum(Nuptake_fin[ls_idKin[[i]]]) )
    PARinonKin <- cbind(PARinonKin, sum(PARi_fin[ls_idnonKin[[i]]]) )
    NuptakenonKin <- cbind(NuptakenonKin, sum(Nuptake_fin[ls_idnonKin[[i]]]) )
    
  }
  MScumvois <- as.numeric(MScumvois)
  PARivois <- as.numeric(PARivois)
  Nuptakevois <- as.numeric(Nuptakevois)
  MScumKin <- as.numeric(MScumKin)
  MScumnonKin <- as.numeric(MScumnonKin)
  PARiKin <- as.numeric(PARiKin)
  NuptakeKin <- as.numeric(PARiKin)
  PARinonKin <- as.numeric(PARinonKin)
  NuptakenonKin <- as.numeric(PARinonKin)
  
  
  dfMS <- data.frame(nump=temptab$nump, MStot_fin, MStot_ini, MStot_coupe1,MStot_coupe2,MStot_coupe3,PARi=PARi_fin, Nuptake=Nuptake_fin, MScumvois, MScumKin, MScumnonKin, PARivois, PARiKin, PARinonKin, Nuptakevois, NuptakeKin, NuptakenonKin)
  #ratio de capture des ressources avec voisins
  dfMS$ratioLight <- PARi_fin/PARivois
  dfMS$ratioNupt <- Nuptake_fin/Nuptakevois
  #EcardPotentiel <-  MStot_fin/mean(MStot_fin) - ParamAllNorm #pas tres logique en multitrait / simple trait
  
  
  temptab <- merge(temptab, dfMS, by="nump")
  
  
  #correlation MSindiv avec valeur des parametres / valeur des voisins / ecart des voisins / ressources / ressources des voisins
  #subx <- temptab[,5:dim(temptab)[2]]#new: avec x,y
  subx <- temptab[,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorAll <- rescor$MStot_fin
  #barplot(valcorAll, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)
  
  #faire un data.frame de ca
  res <- data.frame(t(valcorAll))
  names(res) <- paste("Cor_", row.names(rescor),sep="")
  
  #Corr par espece
  s_temp <- split(temptab, temptab$name)
  sp <- names(s_temp)[1]#"Fix0"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#new: avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp1 <- rescor$MStot_fin
  res1 <- data.frame(t(valcorSp1))
  names(res1) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  
  sp <- names(s_temp)[2]#"Fix1"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#new: avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp2 <- rescor$MStot_fin
  res2 <- data.frame(t(valcorSp2))
  names(res2) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  ##barplot(valcorSp2, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)
  
  res <- cbind(res,res1,res2)
  res$key <- key
  
  res
  
  #renvoie aussi du tableau des donnees : temptab
  ls_resOK <- list(res, temptab)
  names(ls_resOK) <- c("tabCorMSindiv", "datIndices")
  ls_resOK
  
}
#distinguer 2 fonctions? voir 3?: mef et calcul des correlations?
#ls_res_cor_i <- Calc_MSindiv_Corr(ltoto, ls_toto_paquet, ls_paramSD, lspar=c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh"))
#ls_res_cor_i[["tabCorMSindiv"]]
#ls_res_cor_i[["datIndices"]]


