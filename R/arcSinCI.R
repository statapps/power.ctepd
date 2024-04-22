arcSinCI = function(N, R, alpha = 0.05) {
  z = qnorm(1-alpha/2)
  p = R/N
  px = ifelse(R==N, (R-0.5)/N, p)
  a = asin(sqrt(px))
  b = z/(2*sqrt(N))
  p1 = sin(a-b)^2
  p2 = sin(a+b)^2
  p2 = ifelse(R==N, 1.001, p2) ### useful for forest plot.
  return(cbind(p, p1, p2))
}

normCI = function(N, R, alpha = 0.05) {
  z = qnorm(1-alpha/2)
  p = R/N
  px = ifelse(R==N, (R-1)/N, p)
  a = asin(sqrt(px))
  b = z*sqrt(px*(1-px)/N)
  p1 = p - b
  p2 = p + b
  p2 = ifelse(p2>1, 1.00, p2)
  return(cbind(p, p1, p2))
}
