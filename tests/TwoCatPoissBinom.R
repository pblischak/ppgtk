m = 4
m1 = 1
m2 = m - m1
p1 = 0.25
p2 = 0.7
res = rep(0,m+1)

for(x in 0:m){
  
  if(x <= m1){
    
    for(j in 0:x){
      res[x+1] = res[x+1] + dbinom(j, m1, p1, log=T) + dbinom(x-j, m2, p2, log=T)
    }
    
  } else if(x > m1 && x <= m2){
    
    for(j in 0:m1){
      res[x+1] = res[x+1] + dbinom(j, m1, p1, log=T) + dbinom(x-j, m2, p2, log=T)
    }
    
  } else {
    
    for(j in (x-m2):m1){
      res[x+1] = res[x+1] + dbinom(j, m1, p1, log=T) + dbinom(x-j, m2, p2, log=T)
    }
    
  }
  
}

cat(exp(res))
plot(0:m, exp(res), type="h")
