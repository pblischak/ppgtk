
lik <- function(t,r,p){
  
  val <- rep(0,3)
  
  for(a in 0:2){
    gEpsilon <- (a/2) * (1 - 0.005) + (1 - a/2) * 0.005
    val[a] <- exp(dbinom(r, t, gEpsilon, log=T) + dbinom(a, 2, p, log=T))
  }
  
  return(log(sum(val)))
  
}
tt <- rep(200,25)
rr <- rep(20,25)

y <- rep(0, 100)
x <- rep(0, 100)
for(i in 1:100){
  x[i] <- i/101
  for(l in 1:25){
  y[i] <- y[i] + lik(tt[l],rr[l],x[i])
  }
}