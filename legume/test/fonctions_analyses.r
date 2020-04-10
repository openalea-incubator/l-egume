

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
  y <- c(0, resCE_SE$CE+resCE_SE$SE,0)
  col <- rgb(0,0,1,1/4)#4
  polygon(x,y,col=col)
  x <- c(0., resCE_SE$Semprop1, 1.)
  y <- c(0, resCE_SE$CE,0)
  col <- rgb(1,0,0,1/4)#2
  polygon(x,y,col=col)
  points(resCE_SE$Semprop1, resCE_SE$OYa)
  segments(0,0, 1,0, lty=2)
  
  text(0.5,-200,paste("mean CE: ", round(mean(resCE_SE$CE),1)), col=2)
  text(0.5,-250,paste("mean SE: ", round(mean(resCE_SE$SE),1)), col=4)
  text(0.5,-300,paste("mean aOY: ", round(mean(resCE_SE$OYa),1)))
  #pourrait calculer les max!
  
}
#Plot_resCE_SE(resCE_SE, titre="")

