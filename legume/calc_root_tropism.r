
rootTropism <- function(alpha0, g, segment=0.3, Long=10.)
{
  #calcul coordonees de n segments sur longueur Long avec incli initiale = alpha0 et tropisme=g
  
  #alpha0 <- 70. #degre (par rapport a verticale)
  #segment <- 0.3 #cm
  #g <- 0.5 #gravitropisme
  #Long <- 10.#cm
  
  #2 1er point (inclinaison initiale)
  ang_actu <- alpha0
  x <- c(0., segment*cos(alpha0*pi/180.))#axe vertical
  y <- c(0., segment*sin(alpha0*pi/180.))#axe horizontal
  reste_angle <- alpha0
  
  #n points pour faire Long
  for (i in 1:as.integer((Long-segment)/segment))
  {
    ang_actu <- ang_actu - reste_angle*g
    x <- c(x, x[length(x)]+segment*cos(ang_actu*pi/180.)) #axe vertical
    y <- c(y, y[length(y)]+segment*sin(ang_actu*pi/180.)) #axe horizontal
    reste_angle = ang_actu
  }
  cumlen <- seq(0,Long,segment)
  data.frame(x,y,cumlen)    
}

#test <- rootTropism(70,0.5)
#plot(test$x,test$y)

idLong <- function(Long,tabTropism)
{
  #id du dataframe immediatement inferieur a Long en longueur cumulee
  id <- max(which(as.numeric(tabTropism$cumlen)<Long))
  id
}

#id <- idLong(4., test)
#test[id,]

interpolate_rootTropism <- function(alpha0, g, x)
{
  #calcule pour une serie de x
  #calcul des segments
  test <- rootTropism(alpha0, g, segment=0.3, Long=300.)
  
  res <- NULL
  for (xinc in x)
  {
    #xinc = x inconnu 
    idinf <- max(which(test$x<xinc)) #id du poiny de segment juste inferieur
    
    #interpolation lineaire sur le segment
    frac_seg <-  (xinc-test$x[idinf])/ (test$x[idinf+1]-test$x[idinf])
    yinc <- test$y[idinf] + frac_seg*(test$y[idinf+1]-test$y[idinf])
    res <- c(res,yinc)
  }
  #renvoi
  res
  
}

#interpolate_rootTropism(70, 0.5, c(0.5,1,2))
