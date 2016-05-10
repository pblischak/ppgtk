indLik <- function(r, t, e, p){
  
  vals <- rep(NA, 5)
  
  #cat(r,"\t",t,"\n")
  
  for(a in 0:4){
   gEpsilon = (a / 4) * (1 - e) + (1 - a / 4) * e
   vals[a+1] = dbinom(r, t, gEpsilon, log=T) + dbinom(a, 4, p, log=T)
   #cat(t,"\t",r,"\t",dbinom(r, t, gEpsilon, log=T),"\n")
  }
  #cat("----\n")
  
  #cat(sum(exp(vals)),"\n")
  return(exp(vals))
  
}

rr <- as.matrix(read.table("../example/reference.txt",header=F))
tt <- as.matrix(read.table("../example/total.txt",header=F))

dat <- polyfreqs::sim_reads(0.2, 5, 5, 4, 0.005)
lik <- rep(0,25)
surf <- rep(0,100)

for(f in 1:200){
  x <- f/201
  for(i in 1:25){
    lik[i] <- log(sum(indLik(rr[i,1], tt[i,1], 0.005, x)))
  }
  surf[f] <- sum(lik) + 0.5 * log(x) + 0.5 * log(1-x)
}

plot(1:200/201, surf, type="l", xlab="p", ylab="-lnL")
abline(v = which.max(surf)/201, col="blue", lty="dashed", lwd=2)
abline(v = 0.2, col="red")
cat(which.max(surf)/201)