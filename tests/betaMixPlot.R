pi <- rbeta(1,0.5,0.5)
theta1 <- runif(1, 0.001, 1000)
theta2 <- runif(1, 0.001, 1000)
gam <- runif(1, 0, 1)


x <- 1:1000/1001
p <- gam * dbeta(x, pi * theta1, (1 - pi) * theta1) + (1 - gam) * dbeta(x, pi * theta2, (1 - pi) * theta2)
p1 <- dbeta(x, pi * theta1, (1 - pi) * theta1)
p2 <- dbeta(x, pi * theta2, (1 - pi) * theta2)

if(max(p) > max(p1) && max(p) > max(p2)){
  plot(x,p,type="l", lwd=2, xlab="p", ylab="P(p|gamma, theta1, theta2, pi)")
  lines(x,p1,type="l", col="red", lwd=2, lty="dashed")
  lines(x,p2,type="l", col="blue", lwd=2, lty="dashed")
} else if(max(p1) > max(p) && max(p1) > max(p2)){
  plot(x,p1,type="l", col="red", lwd=2, lty="dashed", xlab="p", ylab="P(p|gamma, theta1, theta2, pi)")
  lines(x,p,type="l", lwd=2)
  lines(x,p2,type="l", col="blue", lwd=2, lty="dashed")
} else {
  plot(x,p2,type="l", col="blue", lwd=2, lty="dashed", xlab="p", ylab="P(p|gamma, theta1, theta2, pi)")
  lines(x,p1,type="l", col="red", lwd=2, lty="dashed")
  lines(x,p,type="l", lwd=2)
}