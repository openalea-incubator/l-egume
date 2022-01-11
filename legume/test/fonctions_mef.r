#########
## fonctions de lecture et de mise en forme R pour analyse des sorties de simul l-egume
#########
library(ggplot2)
library(ineq)


read_ltoto <- function(ls_toto)
{
  #recuperer par paquet les fichiers toto du dossier de travail dans une liste ltoto
  ltoto <- vector('list', length(ls_toto))
  names(ltoto) <- ls_toto
  
  for (i in 1:length(ls_toto))
  {
    name <- ls_toto[i]
    ltoto[[name]] <- read.table(name, header=T, sep=';')
  }
  ltoto
}



read_lsSD_MStot <- function(ltoto, ls_paramSD, param_name = "Len")
{
  
  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  #ltoto <- read_ltoto(ls_toto_paquet)
  
  #lit la liste des fichier Sd et les MStot pour une liste de ltoto
  
  ls_MStot <- vector("list",length(ltoto))
  names(ls_MStot) <- names(ltoto)
  
  
  ls_tabSD <- vector("list",length(ltoto))
  names(ls_tabSD) <- names(ltoto)
  
  for (nomfichier in names(ltoto))
  {
    dat <- ltoto[[nomfichier]]
    
    num_usm <- strsplit(nomfichier, '_')[[1]][2]
    scenar <- strsplit(nomfichier, '_')[[1]][6]
    graine <- strsplit(nomfichier, '_')[[1]][8]
    secenarSD <- strsplit(nomfichier, '_')[[1]][10]
    esps <- strsplit(nomfichier, '_')[[1]][4]
    damier <- strsplit(nomfichier, '_')[[1]][5]
    titre <- paste(num_usm, scenar, secenarSD,  damier, graine)#esps, 
    
    #lecture fichier paramSD de l'USM dans tabSD
    nomSD <- ls_paramSD[grepl(paste("paramSD_",num_usm,"_",sep=""), ls_paramSD)] 
    #param_name <- "Len"
    tabSD <- read.table(nomSD, header=T, sep=';')
    
    nb <- dim(dat)[2]-2
    MStot <- dat[dat$V1=='MStot',3:(3+nb-1)] #ajout de MStot
    tabSD$MStotfin <- as.numeric(MStot[dim(MStot)[1],])#derniere ligne
    tabSD$id <- titre
    tabSD$graine <- graine
    
    #split de tabSD par espece et ajout des decile
    sp_tabSD <- split(tabSD, tabSD$name)
    
    sp <- unique(as.character(tabSD$name))[1]#"Fix2"#"nonFixSimTest"#
    valparams <- sp_tabSD[[sp]][,c(param_name)]
    sp_tabSD[[sp]]$decile <- Which_decile(valparams)
    sp <- unique(as.character(tabSD$name))[2]#"nonFixSimTest"#
    valparams <- sp_tabSD[[sp]][,c(param_name)]
    sp_tabSD[[sp]]$decile <- Which_decile(valparams)
    
    tabSD <- do.call("rbind", sp_tabSD)
    
    #stocke dans ls_tabSD et ls_MStot
    ls_tabSD[[nomfichier]] <- tabSD
    ls_MStot[[nomfichier]] <- MStot
  }
  res <- list(ls_tabSD, ls_MStot)
  names(res) <- c("ls_tabSD","ls_MStot")
  res
}



## fonction de mise en formse des simule

moysimval <- function(ltoto, lsusm, var,esp=NA, optSD=F)
{
  # Fait moyenne de la somme pour toute les plantes d'une variable var pour une liste d'usm simulee
  #utilise pour construire le tableau simmoy
  #version GL adapt lucas (v4)
  #optSD=T renvoie standard deviation de la somme des individus
  #esp = NA pour tous le couvert
  #esp pour definir pour une espece du couvert
  
  res <- vector("list",length(lsusm))
  names(res) <- lsusm
  for (usm in lsusm)
  {
    
    if (is.na(esp))
    {dat <- ltoto[[usm]]
    } else
    {
      #garde uniquement col esp
      nomcol <- names(ltoto[[usm]])
      idcols <- grepl(esp, nomcol)
      dat <- cbind(ltoto[[usm]][,c(1:2)], ltoto[[usm]][,idcols])
    }
    
    nbplt <- length(dat)-2
    xplt <- as.matrix(dat[dat$V1==var,3:(3+nbplt-1)], ncol=nbplt)
    xsum <- rowSums(xplt)
    res[[usm]] <- xsum
  }
  if (optSD==F)
  {
    #fait moyenne des sim
    xav <- rowSums(as.data.frame(res))/length(lsusm)
  }else
  {
    #calcule standard deviation des sim
    xav <- apply(as.data.frame(res),MARGIN=1,sd)
  }
  
  xav
}
#LAI <- moysimval(ltoto, lsusm=names(ltoto), var='SurfPlante')/ surfsolref
#LAIsd <- moysimval(ltoto, lsusm=names(ltoto), var='SurfPlante',optSD=T)/ surfsolref






build_simmoy <- function(ltoto, lsusm, esp=NA, optSD=F)
{
  #moy des simul des differentes graines d'un meme usm avec moysimval (pour variables dynamiques)
  
  #recup info generale sur la premier usm
  #dat <- ltoto[[lsusm[1]]]
  if (is.na(esp))
  {dat <- ltoto[[lsusm[1]]]
  } else
  {
    #garde uniquement col esp
    nomcol <- names(ltoto[[lsusm[1]]])
    idcols <- grepl(esp, nomcol)
    dat <- cbind(ltoto[[lsusm[1]]][,c(1:2)], ltoto[[lsusm[1]]][,idcols])
  }
  
  TT <- dat[dat$V1=='TT',3] #peut changer selon les plantes!
  STEPS <- dat[dat$V1=='TT',2]
  nbplt <- length(dat)-2
  surfsolref <- dat[dat$V1=='pattern',3] #m2

  LAI <- moysimval(ltoto, lsusm, var='SurfPlante', esp, optSD)/ surfsolref
  MSA <- moysimval(ltoto,lsusm, var='MSaerien', esp, optSD)/ surfsolref
  MSArec <- moysimval(ltoto,lsusm, var='MSaerienRec', esp, optSD)/ surfsolref
  MSAnonrec <- moysimval(ltoto,lsusm, var='MSaerienNonRec', esp, optSD)/ surfsolref
  MSpiv <- moysimval(ltoto,lsusm, var='MS_pivot', esp, optSD)/ surfsolref
  MSracfine <- moysimval(ltoto,lsusm, var='MS_rac_fine', esp, optSD)/ surfsolref
  MSrac <- MSpiv + MSracfine
  NBI <- moysimval(ltoto,lsusm, var='NBI', esp, optSD)/ nbplt
  NBI <- pmax(0, NBI - 0.75) #correction des simuls pour les comptages decimaux
  #NBIquart <- quantsimval(ltoto,lsusm, var_='NBI',esp=esp)
  NBphyto <- moysimval(ltoto, lsusm, var='NBphyto', esp, optSD)/ surfsolref
  Nbapex <- moysimval(ltoto, lsusm, var='NBapexAct', esp, optSD)/ surfsolref
  NBphyto <- pmax(0,NBphyto - 0.5*Nbapex) #correction simuls pour les comptages decimaux
  NBsh <- moysimval(ltoto, lsusm, var='NBsh', esp, optSD)/ surfsolref
  
  RDepth <- moysimval(ltoto,lsusm, var='RDepth', esp, optSD)/ nbplt
  Hmax <- moysimval(ltoto,lsusm, var='Hplante', esp, optSD)/ nbplt
  FTSW <- moysimval(ltoto,lsusm, var='FTSW', esp, optSD)/ nbplt
  NNI <- moysimval(ltoto,lsusm, var='NNI', esp, optSD)/ nbplt
  R_DemandC_Root <- moysimval(ltoto,lsusm, var='R_DemandC_Root', esp, optSD)/ nbplt
  cutNB <- moysimval(ltoto,lsusm, var='cutNB', esp, optSD)/ nbplt
  Npc_aer <- moysimval(ltoto,lsusm, var='Npc_aer', esp, optSD)/ nbplt
  Ndfa <- moysimval(ltoto,lsusm, var='Ndfa', esp, optSD)/ nbplt
  Epsi <- moysimval(ltoto,lsusm, var='epsi', esp, optSD)
  
  simmoy <- data.frame(STEPS, TT, NBI, NBphyto, LAI, MSA, MSArec, MSAnonrec, MSpiv, MSracfine, MSrac, RDepth, Hmax, FTSW, NNI, R_DemandC_Root, cutNB, Npc_aer,Ndfa,Epsi,NBsh)
  simmoy
}#version revue par Lucas tient cmpte du nom de l'espece dans les assos

#simmoy <- build_simmoy(ltoto, lsusm=names(ltoto))
#simmoy <- build_simmoy(ltoto, lsusm=names(ltoto), esp="timbale")


dynamic_graphs <- function(simmoy, name, obs=NULL, surfsolref=NULL)
{
  #serie de figure pour rapport dynamique d'une simulation Avec ou sans ajouts de points observe
  #pour obs ajoute si dataframe fourni au format obs (teste morpholeg)
  
  op <- par(mfrow = c(3,1), #lignes, colonnes
            oma = c(5.,2.,3,0) + 0.1, #outer margins c(bottom, left, top, right)
            mar = c(0,4,0,2) + 0.1) #marges internes a chaque compartiment c(bottom, left, top, right)
  
  
  #1) Leaf area components
  plot(simmoy$STEPS, simmoy$NBI, type='l', xlab='Time',ylab='Nb phytomere I', labels=F, ylim=c(0,1.5*max(simmoy$NBI)))
  axis(2,labels=T) #remet tick labels y
  title(main=name, outer=T)
  if (!is.null(obs) & 'NBI' %in% names(obs)) 
  {  points(obs$DOY, obs$NBI, pch=16) }
  
  
  plot(simmoy$STEPS, simmoy$NBphyto, type='l', xlab='Time',ylab='Nb phytomere tot', labels=F, ylim=c(0,1.5*max(simmoy$NBphyto)))
  axis(2,labels=T) #remet tick labels y
  if (!is.null(obs) & 'nb_phyto_tot' %in% names(obs))
  { points(obs$DOY, obs$nb_phyto_tot/surfsolref, pch=16) }
  
  
  plot(simmoy$STEPS, simmoy$LAI, type='l', xlab='Time',ylab='LAI', labels=F, ylim=c(0,1.5*max(simmoy$LAI)))
  axis(2,labels=T) #remet tick labels y
  axis(1,labels=T) #remet tick labels x
  title(xlab='DOY', outer=T)
  if (!is.null(obs) & 'LAI' %in% names(obs))
  { points(obs$DOY, obs$LAI, pch=16) } 
  if (!is.null(obs) & 'surf_tot' %in% names(obs))
  { points(obs$DOY, obs$surf_tot/ (10000*surfsolref), pch=16) } #a reprendre fait 2 courbes actuellement pour eviter bug
  
  #2)MS et taille
  plot(simmoy$STEPS, -simmoy$RDepth, type='l', xlab='Time',ylab='RDepth', labels=F, ylim=c(-1.5*max(simmoy$RDepth),0))
  axis(2,labels=T) #remet tick labels y
  title(main=name, outer=T)
  if (!is.null(obs) & 'long_pivot' %in% names(obs))
  { points(obs$DOY, -obs$long_pivot, pch=16)}
  
  
  plot(simmoy$STEPS, simmoy$Hmax, type='l', xlab='Time',ylab='Hmax', labels=F, ylim=c(0,1.5*max(simmoy$Hmax)))
  axis(2,labels=T) #remet tick labels y
  if (!is.null(obs) & 'Hmax' %in% names(obs))
  { points(obs$DOY, obs$Hmax, pch=16) }
  
  plot(simmoy$STEPS, simmoy$MSA, type='l', xlab='Time',ylab='MS', labels=F, ylim=c(0,1.5*max(simmoy$MSA)))
  axis(2,labels=T) #remet tick labels y
  axis(1,labels=T) #remet tick labels x
  title(xlab='DOY', outer=T)
  points(simmoy$STEPS, simmoy$MSrac, type='l', lty=2)
  if (!is.null(obs) & 'MSaerien' %in% names(obs) & 'MSroot_tot' %in% names(obs))
  { 
    points(obs$DOY, obs$MSaerien/surfsolref, pch=16)
    points(obs$DOY, obs$MSroot_tot/surfsolref)
  }
  
  
  #3) fonctions de stress
  plot(simmoy$STEPS, simmoy$FTSW, type='l', xlab='Time',ylab='FTSW', labels=F, ylim=c(0,1.1))
  axis(2,labels=T) #remet tick labels y
  title(main=name, outer=T)
  
  plot(simmoy$STEPS, simmoy$NNI, type='l', xlab='Time',ylab='NNI', labels=F, ylim=c(0, 1.2*max(simmoy$NNI)))
  axis(2,labels=T) #remet tick labels y
  
  plot(simmoy$STEPS, simmoy$R_DemandC_Root, type='l', xlab='Time',ylab='R_DemandC_Root', labels=F, ylim=c(0,1.1))
  axis(2,labels=T) #remet tick labels y
  axis(1,labels=T) #remet tick labels x
  title(xlab='DOY', outer=T)
}
#dynamic_graphs(simmoy, onglet, obs, surfsolref)

#sans points d'obsevation
#dynamic_graphs(simmoy, name=names(ltoto)[1]) 

#avec points d'observation...
#namexl <- "morpholeg14_obs.xls"
#onglet <- "morpholeg14_ISO_timbale"
#obs <- read_excel(paste(pathobs,namexl,sep="\\"), sheet = onglet, col_names = TRUE, na = "")
#dynamic_graphs(simmoy, name=names(ltoto)[1], obs=obs, surfsolref=1) 



#fonction graph ggplot2
gg_plotsim <- function(varsim, simmoy, simsd, name="", col="blue")
{
  #fait le line plot avec ecart type a partir des tableau moyen simules
  var_ <- varsim #"FTSW"#"NBI"#"MSA"#"NNI"#"LAI"#
  
  min <- 0
  max <- 1.5*max(simmoy[,var_])
  
  plot_var <- ggplot(data = simmoy, aes(x = STEPS)) +
    geom_line(aes(y = simmoy[,var_]), color=col)+
    geom_ribbon(aes(ymin=simmoy[,var_]-simsd[,var_],ymax=simmoy[,var_]+simsd[,var_]),fill="blue",alpha=0.2)+
    geom_hline(yintercept=0)+
    ylim(min,max)+
    geom_text(x=1.20*min(simmoy$STEPS), y=0.98*max, size=4, label=name)+
    theme(axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))+
    labs(title = "obs",subtitle = "sim",x = "DOY", y = var_)+
    theme(plot.title=element_text(size=10,color = "red"),plot.subtitle = element_text(size=10,color = "blue"))
  
  plot_var
}
#gg_plotsim("LAI", simmoy, simsd, "test")
#titre sera a revoir...



gg_addplotobs <- function(plot_var, var_, obsOK, corresp, colobs="red")
{
  # ajour a un graph simule des points observe pour variable var_
  #obsOK : obs avec tableau meme dimension que les simul (merge)
  #corresp: dataframe de correspondance des nom de variables obs/sim
  nomvarobs <- as.character(corresp[corresp$sim==var_,c("obs")])
  plot_var2 <- plot_var + {if(var_ %in% corresp$sim) geom_point(aes(obsMerge$DOY, obsMerge[,nomvarobs]), fill=colobs,color=colobs , size=2)}
  
  plot_var2
}
#gg_addplotobs(ls_plt[["MSA"]], "MSA", obsMerge, corresp)




gg_plotObsSim <- function(obssim, var_, name="", colpt="red")
{
  #plot obs-sim avec ggplot
  
  #var_ <- "NBI"#"MSA"#
  #nomvarobs <- as.character(corresp[corresp$sim==var_,c("obs")])
  #obssim <- na.omit(data.frame(obs=simmoy[,var_], sim=obsMerge[,nomvarobs]))
  #name <- onglet
  
  
  min <- 0
  max <- 1.5*max(obssim$sim)
  reg   <- lm(obs ~ sim, data = obssim)
  coeff <- coefficients(reg)
  eq    <- paste0("y = ", round(coeff[2],2), "x + ", round(coeff[1],2))
  RMSE_ <- round(rmse(obssim$obs,obssim$sim), 2)
  rmses_ <- rmsesCoucheney(obssim$obs,obssim$sim)
  rmse_ <- rmseuCoucheney(obssim$obs,obssim$sim)
  pRMSEs_ <- round(pRMSEs(rmse_, rmses_), 2)
  rRMSE <- round(rrmseCoucheney(obssim$obs,obssim$sim), 2)
  EF <- round(efficiencyCoucheney(obssim$obs,obssim$sim), 2)
  
  plot_ObsSim <- ggplot(obssim, aes(x = obs, y = sim)) +
    ggtitle(name)+
    geom_abline(intercept = 0, slope = 1, color = "black")+
    geom_point(aes(color = "obs"))+
    geom_smooth(method=lm, se = FALSE, color = colpt)+
    ylim(min,max)+
    xlim(min,max)+
    geom_text(x=0.2*max, y=0.95*max, size=3, label=eq)+
    geom_text(x=0.2*max, y=0.9*max, size=3,label=paste("RMSE: ",RMSE_))+
    geom_text(x=0.2*max, y=0.85*max, size=3,label=paste("rRMSE: ",rRMSE))+
    geom_text(x=0.2*max, y=0.8*max, size=3,label=paste("pRMSEs: ",pRMSEs_))+
    geom_text(x=0.2*max, y=0.75*max, size=3,label=paste("EF: ",EF))+
    labs(x = paste("Obs ", var_), y = paste("Sim ", var_))
  
  
  plot_ObsSim 
}
#plot_ObsSim <- gg_plotObsSim(obssim, "NBI", name=onglet)




mef_dosbssim <- function(var, varsim, obs, simmoy, name='', corobs=1., cutNB=0.)
{
  #mise en forme d'un dictionnaire observe, simule a partir du tableau des obs et du simmoy et des noms de variables var/varsim
  #suppose pas que variable observe et simulee ont le meme nom!
  idvar <- which(names(obs)==var)
  
  obs[obs[,idvar] < 0. & !is.na(obs[,idvar]), idvar] <- NA #retire les -999
  DOYmes <- obs$DOY[!is.na(obs[,idvar])] #recupere les liste de DOY des points mesures sans les NA
  obsvar <- obs[!is.na(obs[,idvar]), idvar]
  
  #recupere les DOY de simulation correspondant
  idDOYsim <- simmoy$STEPS %in% DOYmes
  idvar_sim <- which(names(simmoy)==varsim)
  simvar <- simmoy[idDOYsim, idvar_sim]
  
  #gestion des cas ou pas de simul face aux obs
  if (length(simvar)<length(obsvar))
  {
    #DOYmes qui n'ont pas de sim
    DOYsimOK <-simmoy[idDOYsim, 1]
    DOYsimKO <- DOYmes[! DOYmes %in% DOYsimOK]
    idDOYsimOK <- DOYmes %in% DOYsimOK
    obsvar <- obsvar[idDOYsimOK]
    DOYmes <- DOYsimOK
    #afficher un warnings!
  }
  
  
  if (length(simvar)>0)
  {dobssim <- data.frame(usm= rep(name, length(obsvar)), var=rep(var, length(obsvar)) , DOY=DOYmes, obs=obsvar*corobs, sim=simvar,corobs=rep(corobs, length(obsvar))) 
  names(dobssim)[4] <- "obs" #car bug de changement de nom de colonne
  }
  else
  {dobssim <- NULL
  }
  
  dobssim
}


#var <-'NBI'#ls_var <- c('NBI','nb_phyto_tot','surf_tot','long_pivot','Hmax','MSaerien')
#varsim <- 'NBI'#ls_varsim <- c('NBI','NBphyto','LAI', 'RDepth', 'Hmax','MSA')
#dobssim <- mef_dosbssim(var, varsim, obs, simmoy, corobs=1.)



merge_dobssim <- function(ls_dobssim)
{
  #reuni dans un seul dataframe les differents dobssim d'une liste
  dobssim <- ls_dobssim[[1]]
  if (length(ls_dobssim)>1)
  {
    for (i in 2:length(ls_dobssim)) 
    {
      dobssim <- rbind(dobssim, ls_dobssim[[i]])
    }
  }
  dobssim
}



plot_obssim <- function(dobssim, name='', displayusm=F)
{
  # fonction R obs-sim
  #conversion
  #dobssim$obs <- dobssim$obs*convert #suppose conversion la meme pour tous les obs!!!
  #dobssim deja corrige maintenat!
  
  #plot d'un obs-sim avec les indicateurs de qualite
  #calcul des indicateurs
  res_rmse <- signif(rmse(dobssim$obs,dobssim$sim),3) #RMSE
  res_rrmse <- signif(rrmseCoucheney(dobssim$obs,dobssim$sim),3)#rRMSE
  res_rmses <- rmsesCoucheney(dobssim$obs,dobssim$sim)#RMSEs
  res_rmseu <- rmseuCoucheney(dobssim$obs,dobssim$sim)#RMSEu
  res_prmses <- signif(pRMSEs(rmse(dobssim$obs,dobssim$sim), res_rmses),3)#pRMSEs
  res_prmseu <- signif(pRMSEu(rmse(dobssim$obs,dobssim$sim), res_rmseu),3)#pRMSEu
  res_EF <- signif(efficiencyCoucheney(dobssim$obs,dobssim$sim),3)#Efficiency
  res_r2 <- signif(summary(lm(dobssim$sim~dobssim$obs))$r.squared,3)#r2
  
  
  plot(dobssim$obs, dobssim$sim, xlim=c(0, 1.5*max(dobssim$obs)), ylim=c(0, 1.5*max(dobssim$obs)), xlab='obs', ylab='sim', main=name)
  points(c(0, 1.5*max(dobssim$obs)), c(0, 1.5*max(dobssim$obs)), type='l')
  if (length(dobssim$obs)>2)
  { abline(lm(dobssim$sim~dobssim$obs), col=2)}
  text(0.2*max(dobssim$obs) , 1.4*max(dobssim$obs), cex=0.8, paste("rmse: ", res_rmse))
  text(0.2*max(dobssim$obs) , 1.3*max(dobssim$obs), cex=0.8, paste("rRMSE: ", res_rrmse))
  text(0.2*max(dobssim$obs) , 1.2*max(dobssim$obs), cex=0.8, paste("pRMSEs: ", res_prmses))
  text(0.2*max(dobssim$obs) , 1.1*max(dobssim$obs), cex=0.8, paste("pRMSEu: ", res_prmseu))
  text(0.2*max(dobssim$obs) , 1.*max(dobssim$obs), cex=0.8, paste("EF: ", res_EF))
  text(0.2*max(dobssim$obs) , 0.9*max(dobssim$obs), cex=0.8, paste("R2: ", res_r2))
  
  #ajout des points par usm si option activee (completer dico_col)
  if (displayusm==T)
  {
    splt <- split(dobssim, dobssim$usm)
    #ls_pch <- c(1,16,1,16,1,16)
    #ls_col <- c(1,1,2,2,3,3)
    dico_col <- data.frame(expe=c('morpholeg14', 'morpholeg15','combileg15','combileg16','morpholegRGR15'), col=c(1,1,2,2,3),pch=c(16,1,1,16,1))
    
    for (i in 1:length(splt))
    {
      expe <- strsplit(names(splt)[i], '_')[[1]][1]
      points(splt[[i]]$obs, splt[[i]]$sim, pch=dico_col$pch[dico_col$expe==expe], col=dico_col$col[dico_col$expe==expe])
    }
  }#rq: code couleur change selon les variables, selon nb ds split ; pour garder cst faudrait un dico
  
  
}

#plot_obssim(dobssim, name='test', displayusm=T)





build_ls_dobssim <-function(esp_, ls_expe, ls_var, ls_varsim)
{
  #construit un arbre (liste de liste) de tableau dobssim pour une espece et plusieurs variables / plusieurs usm-expe
  
  #construit un ls_dobssim vide pour une espece et plusieurs expe et variables
  ls_dobssim <- vector('list', length(ls_varsim))
  names(ls_dobssim) <- ls_varsim
  for (i in ls_varsim)
  {
    ls_dobssim[[i]] <- vector('list', length(ls_expe))
    names(ls_dobssim[[i]]) <- ls_expe
  }
  
  #boucle pour toutes les variables et les usm
  for (i in 1:length(ls_varsim))
    try({
      var <- ls_var[i]#'surf_tot'#'nb_phyto_tot'#'NBI'#
      varsim <- ls_varsim[i]#'LAI'#'NBphyto'#'NBI'#
      
      for (key in ls_expe)
        try({
          #calcul simul moy
          ls_toto_paquet <- sp_dtoto[[key]]$name
          
          #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
          ltoto <- read_ltoto(ls_toto_paquet)
          #version locale du paquet de doto
          dtoto <- sp_dtoto[[key]]
          
          mix <- strsplit(ls_toto_paquet[1], '_')[[1]][4] #suppose paquet fait par traitement
          esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
          esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'
          meteo <- strsplit(ls_toto_paquet[1], '_')[[1]][9]
          damier <- strsplit(ls_toto_paquet[1], '_')[[1]][5]
          
          
          #constrcution du dobsim
          #calcul de la moyenne des simuls pour esp_
          simmoy <- build_simmoy(ltoto, lsusm=names(ltoto), esp_)  
          
          #recup des obs correspondant
          namexl <- paste0(meteo, "_obs.xls")#"morpholeg14_obs.xls"
          trait <- if (esp == esp2 & damier=="homogeneous0") "ISO" else "HD-M2"#sera a reprndre pour diffeents traitement
          #trait <- if (esp == esp2 & damier=="homogeneous0" & meteo=="morpholegRGR15") "LD"
          onglet <- paste0(meteo, "_",trait,"_",esp_)#"morpholeg14_ISO_timbale" #marche pour isole; a revoir pour autres
          obs <- read_excel(paste(pathobs,namexl,sep="\\"), sheet = onglet, col_names = TRUE, na = "")
          
          
          #calcul des facteurs de correction selon varsim
          if (varsim=='NBI') {corr=1.}
          if (varsim=='NBphyto') {corr=1./surfsolref}        
          if (varsim=='LAI') {corr=1/(10000*surfsolref)}
          if (varsim=='RDepth') {corr=1.}
          if (varsim=='Hmax') {corr=1.}
          if (varsim=='MSA') {corr=1./surfsolref}
          
          
          ls_dobssim[[varsim]][[key]] <- mef_dosbssim(var, varsim, obs, simmoy, name=onglet, corobs=corr)
        })
    })
  ls_dobssim
}

#esp_ <- 'timbale'#'giga'#'formica'#'sevanskij'#'leo'#'canto'#'kayanne'#
#ls_expe <- c('morpholeg14', 'morpholeg15')
#ls_expe <- names(sp_dtoto)[grepl(esp_, names(sp_dtoto))]#cle comportant le bon geno
#ls_var <- c('NBI','nb_phyto_tot','surf_tot','Hmax','MSaerien')#,'long_pivot')
#ls_varsim <- c('NBI','NBphyto','LAI', 'Hmax','MSA')#, 'RDepth')
#ls_dobssim <- build_ls_dobssim(esp, ls_expe, ls_var, ls_varsim)




build_dtoto <- function(sp_dtoto, key, DOYdeb, DOYScoupe)
{
  ls_toto_paquet <- sp_dtoto[[key]]$name

  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  ltoto <- read_ltoto(ls_toto_paquet)
  #version locale du paquet de doto
  dtoto <- sp_dtoto[[key]]

  #recup du nom des esp
  mix <- strsplit(ls_toto_paquet[1], '_')[[1]][4] #suppose paquet fait par traitement
  esp <- strsplit(mix, '-')[[1]][1] #'Fix2'
  esp2 <- strsplit(mix, '-')[[1]][2] #'nonFixSimTest'

  #visu des rendement moyen m2 / a un DOY
  surfsolref <- NULL
  nbplt <- NULL
  nbplt1 <- NULL
  nbplt2 <- NULL

  #DOYScoupe <- c(165,199,231,271,334)#Avignon
  #DOYScoupe <- c(187,229,282,334)#Lusignan
  #DOYdeb <- 60
  idDOYScoupe <- DOYScoupe - DOYdeb
  Ytot <- NULL
  Ycoupe <- NULL

  YEsp1 <- NULL
  YEsp2 <- NULL

  QNfix <- NULL
  QNupttot <- NULL
  QNuptleg <- NULL

  PARi1 <- NULL
  PARi2 <- NULL
  Surf1 <- NULL
  Surf2 <- NULL
  LRac1 <- NULL
  LRac2 <- NULL
  MRac1 <- NULL
  MRac2 <- NULL
  MGini1 <- NULL
  MGini2 <- NULL
  MAlive1 <- NULL
  MAlive2 <- NULL

  for (i in 1:length(ls_toto_paquet))#(ls_toto))
  {
    name <- ls_toto_paquet[i]
    damier <- strsplit(name, '_')[[1]][5]
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
    #if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    #{
    #  idcols <- as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
    #} else
    #{
    idcols <- grepl(esp, nomcol) & !grepl(esp2, nomcol)#contient esp1 et pas esp2
    #}

    dat1 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
    nb1 <- length(dat1)-2
    nbplt1 <- cbind(nbplt1, nb1)
    if (nb1>0)
    {
      MS1 <- as.matrix(dat1[dat1$V1=='MSaerien' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)
      ProdIaer1 <- rowSums(MS1) / s
      Nuptake_sol_leg <- as.matrix(dat1[dat1$V1=='Nuptake_sol',3:(3+nb1-1)], ncol=nb1)
      Nuptake_sol_leg <- as.numeric(rowSums(Nuptake_sol_leg) / s)
      jPARi1 <- rowSums(as.matrix(dat1[dat1$V1=='PARiPlante' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jSurf1 <- rowSums(as.matrix(dat1[dat1$V1=='SurfPlante' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jLRac1 <- rowSums(as.matrix(dat1[dat1$V1=='RLTot' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jMRac1 <- rowSums(as.matrix(dat1[dat1$V1=='MS_rac_fine' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      jMPiv1 <- rowSums(as.matrix(dat1[dat1$V1=='MS_pivot' & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)) / s
      #gini par date sur MSA (ttes les plantes)
      Gini1 <- NULL
      for (k in 1:length(DOYScoupe))
      { Gini1 <- cbind(Gini1, ineq(MS1[k,], type="Gini"))}
      
      #survie
      matSV1 <- as.matrix(dat1[dat1$V1 == "aliveB" & dat1$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb1)
      matSV1[matSV1>0] <- 1
      aliveP1 <-  as.numeric(nbplt)-as.matrix(rowSums(matSV1))
      aliveDens1 <- t(aliveP1 / s)
      #MortDens1 <-  as.matrix(rowSums(matSV1)) / s
      

    } else
    {
      ProdIaer1 <- 0 #pas de plante de l'esp1
      Nuptake_sol_leg <- 0
      jPARi1 <- 0
      jSurf1 <- 0
      jLRac1 <- 0
      jMRac1 <- 0
      jMPiv1 <- 0
      Gini1 <- c(NA,NA,NA,NA)
      aliveDens1 <- c(0,0,0,0)
    }
    YEsp1 <- cbind(YEsp1, sum(ProdIaer1))#cumul des 5 coupes
    QNuptleg <- cbind(QNuptleg, sum(Nuptake_sol_leg))
    PARi1 <- cbind(PARi1, sum(jPARi1))
    Surf1 <- cbind(Surf1, sum(jSurf1))
    LRac1 <- cbind(LRac1, max(jLRac1))
    MRac1 <- cbind(MRac1, max(jMRac1)+max(jMPiv1))
    MGini1 <- rbind(MGini1, Gini1)
    MAlive1 <- rbind(MAlive1, aliveDens1)

    #YEsp2
    #if (esp==esp2 & grep('damier', damier)==1)#si deux fois le meme nom d'espece, mais mixture damier
    #{
    #  idcols <- !as.logical(didcols[didcols$damier==damier,1:66])#idcols lu dans fichier qui leve les ambiguite
    #  idcols[1:2] <- FALSE #remet a faux les deux premieres colonnes
    #} else
    #{
    idcols <- grepl(esp2, nomcol)#contient esp2
    #}

    dat2 <- cbind(ltoto[[name]][,c(1:2)], ltoto[[name]][,idcols])
    nb2 <- length(dat2)-2
    nbplt2 <- cbind(nbplt2, nb2)
    if (nb2>0)
    {
      MS2 <- as.matrix(dat2[dat2$V1=='MSaerien' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)
      ProdIaer2 <- rowSums(MS2) / s
      jPARi2 <- rowSums(as.matrix(dat2[dat2$V1=='PARiPlante' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jSurf2 <- rowSums(as.matrix(dat2[dat2$V1=='SurfPlante' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jLRac2 <- rowSums(as.matrix(dat2[dat2$V1=='RLTot' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jMRac2 <- rowSums(as.matrix(dat2[dat2$V1=='MS_rac_fine' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      jMPiv2 <- rowSums(as.matrix(dat2[dat2$V1=='MS_pivot' & dat2$steps %in% DOYScoupe,3:(3+nb2-1)], ncol=nb2)) / s
      #gini par date sur MSA (ttes les plantes)
      Gini2 <- NULL
      for (k in 1:length(DOYScoupe))
      { Gini2 <- cbind(Gini2, ineq(MS2[k,], type="Gini"))}
      
      #survie
      matSV2 <- as.matrix(dat2[dat2$V1 == "aliveB" & dat2$steps %in% DOYScoupe,3:(3+nb1-1)], ncol=nb2)
      matSV2[matSV2>0] <- 1
      aliveP2 <-  as.numeric(nbplt)-as.matrix(rowSums(matSV2))
      aliveDens2 <- t(aliveP2 / s)
      #MortDens2 <-  as.matrix(rowSums(matSV2)) / s
    }
    else
    {
      ProdIaer2 <- 0 #pas de plante de l'esp2
      jPARi2 <- 0
      jSurf2 <- 0
      jLRac2 <- 0
      jMRac2 <- 0
      jMPiv2 <- 0
      Gini2 <- c(NA,NA,NA,NA)
      aliveDens2 <- c(0,0,0,0)
    }
    YEsp2 <- cbind(YEsp2, sum(ProdIaer2))#cumul des 5 coupes
    #YEsp2 <- cbind(YEsp2, sum(ProdIaer2))#cumul des 5 coupes
    PARi2 <- cbind(PARi2, sum(jPARi2))
    Surf2 <- cbind(Surf2, sum(jSurf2))
    LRac2 <- cbind(LRac2, max(jLRac2))
    MRac2 <- cbind(MRac2, max(jMRac2)+max(jMPiv2))
    MGini2 <- rbind(MGini2, Gini2)
    MAlive2 <- rbind(MAlive2, aliveDens2)
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

  #new var
  dtoto$Pari1 <- as.numeric(PARi1)
  dtoto$Pari2 <- as.numeric(PARi2)
  dtoto$Surf1 <- as.numeric(Surf1)
  dtoto$Surf2 <- as.numeric(Surf2)
  dtoto$PhiSurf1 <- as.numeric(PARi1) / (as.numeric(Surf1) + 10e-12)#Phi Surf
  dtoto$PhiSurf2 <- as.numeric(PARi2) / (as.numeric(Surf2) + 10e-12)
  dtoto$PhiMass1 <- as.numeric(PARi1) / (as.numeric(YEsp1) + 10e-12)#Phi Mass
  dtoto$PhiMass2 <- as.numeric(PARi2) / (as.numeric(YEsp2) + 10e-12)
  dtoto$LRac1 <- as.numeric(LRac1)
  dtoto$LRac2 <- as.numeric(LRac2)
  dtoto$MRac1 <- as.numeric(MRac1)
  dtoto$MRac2 <- as.numeric(MRac2)
  dtoto$UptNLen1 <- (as.numeric(QNupttot) - as.numeric(QNuptleg)) / (as.numeric(LRac1) + 10e-12)#Uptake par Len
  dtoto$UptNLen2 <- as.numeric(QNuptleg) / (as.numeric(LRac2) + 10e-12)
  dtoto$UptNMass1 <- (as.numeric(QNupttot) - as.numeric(QNuptleg)) / (as.numeric(MRac1) + 10e-12)#Uptake par Mass root
  dtoto$UptNMass2 <- as.numeric(QNuptleg) / (as.numeric(MRac2) + 10e-12)
  dtoto$gini1 <- rowMeans(MGini1) #moyenne des gini de ttes les dates (sans retirer pltes mortes)
  dtoto$gini2 <- rowMeans(MGini2)
  dtoto$alive1 <- as.numeric(MAlive1[,dim(MAlive1)[2]]) #survie derniere date coupe
  dtoto$alive2 <- as.numeric(MAlive2[,dim(MAlive2)[2]]) #survie derniere date coupe
  
  dtoto

}
#fonction a generaliser et a bouger ailleurs




#tab <- ltoto[[id]]
#unique(tab$V1)

#nbplt <- dim(tab)[2] - 2
#matSV <- as.matrix(tab[tab$V1 == "aliveB",3:(nbplt+2)])
#matSV[matSV>0] <- 1
#surfsol <- dtoto$surfsolref[id]

#DOYs <- tab[tab$V1 == "aliveB",2]
#aliveP <- nbplt - rowSums(matSV)
#aliveDens <- aliveP / surfsol






















##old vesions (GL 2)
moysimval1 <- function(ltoto, lsusm, var)
{
  # Fait moyenne de la somme pour toute les plantes d'une variable var pour une liste d'usm simulee
  #utilise pour construire le tableau simmoy
  #version initiale GL (v2)
  
  res <- vector("list",length(lsusm))
  names(res) <- lsusm
  for (usm in lsusm)
  {
    dat <- ltoto[[usm]]
    nbplt <- length(dat)-2
    xplt <- as.matrix(dat[dat$V1==var,3:(3+nbplt-1)], ncol=nbplt)
    xsum <- rowSums(xplt)
    res[[usm]] <- xsum
  }
  xav <- rowSums(as.data.frame(res))/length(lsusm)
  xav
}
#LAI <- moysimval(ltoto, lsusm, var='SurfPlante')/ surfsolref



build_simmoy1 <- function(ltoto, lsusm)
{
  #moy des simul des differentes graines d'un meme usm avec moysimval (pour variables dynamiques)
  
  dat <- ltoto[[lsusm[1]]]
  TT <- dat[dat$V1=='TT',3] #peut changer selon les plantes!
  STEPS <- dat[dat$V1=='TT',2]
  nbplt <- length(dat)-2
  surfsolref <- dat[dat$V1=='pattern',3] #m2
  
  LAI <- moysimval(ltoto, lsusm, var='SurfPlante')/ surfsolref
  MSA <- moysimval(ltoto,lsusm, var='MSaerien')/ surfsolref
  MSpiv <- moysimval(ltoto,lsusm, var='MS_pivot')/ surfsolref
  MSracfine <- moysimval(ltoto,lsusm, var='MS_rac_fine')/ surfsolref
  MSrac <- MSpiv + MSracfine
  NBI <- moysimval(ltoto,lsusm, var='NBI')/ nbplt
  NBI <- pmax(0, NBI - 0.75) #correction des simuls pour les comptages decimaux
  #NBIquart <- quantsimval(ltoto,lsusm, var_='NBI',esp=esp)
  NBphyto <- moysimval(ltoto, lsusm, var='NBphyto')/ surfsolref
  Nbapex <- moysimval(ltoto, lsusm, var='NBapexAct')/ surfsolref
  NBphyto <- pmax(0,NBphyto - 0.5*Nbapex) #correction simuls pour les comptages decimaux
  
  RDepth <- moysimval(ltoto,lsusm, var='RDepth')/ nbplt
  Hmax <- moysimval(ltoto,lsusm, var='Hplante')/ nbplt
  FTSW <- moysimval(ltoto,lsusm, var='FTSW')/ nbplt
  NNI <- moysimval(ltoto,lsusm, var='NNI')/ nbplt
  R_DemandC_Root <- moysimval(ltoto,lsusm, var='R_DemandC_Root')/ nbplt
  cutNB <- moysimval(ltoto,lsusm, var='cutNB')/ nbplt
  
  simmoy <- data.frame(STEPS, TT, NBI, NBphyto, LAI, MSA, MSpiv, MSracfine, MSrac, RDepth, Hmax, FTSW, NNI, R_DemandC_Root, cutNB)
  simmoy
}#version revue par Lucas tient cmpte du nom de l'espece dans les assos

#simmoy <- build_simmoy(ltoto, lsusm=names(ltoto))


#pour gestion des couleur: vecteur 100
col100 <- function(valrel100, lscols)
{
  # fonction pour definir un vecteur de couleur a partir de valeur relative et d'une liste de 101 couleur
  #lscols = vecteur de 100 couleurs
  # valrel100 = position dans ce vecteur (% du max)
  
  #lscols[rdtrel]#pas bon!
  cols_ <- NULL
  for(i in valrel100)
  {
    cols_ <- rbind(cols_, lscols[i+1])
  }
  cols_ <- as.vector(cols_)
  #en faire une fonction
  cols_
}

