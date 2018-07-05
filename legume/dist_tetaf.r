disttetaf <- function (mf, sdf, n=10000,seed=0)
{
  #calcule une distrib par classe connaisant angle moyen et ecartpe
  set.seed(seed)
  x <- rnorm(n, mean = mf, sd = sdf)
  x[which(x<0)] <- -x[which(x<0)]
  x[which(x>90)] <- 90-(x[which(x>90)]-90)
  x[which(x<0)] <- -x[which(x<0)]
  x[which(x>90)] <- 90-(x[which(x>90)]-90)
  res <- hist(x, breaks=c(0,10,20,30,40,50,60,70,80,90), plot=F)
  res$counts/n
}
#disttetaf (5,10)
