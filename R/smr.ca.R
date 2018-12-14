###Below is a method using Exponential approximation
#
# Modify the file from SMR in survexp.fr package
#
data(survexp.ca)

smr.ca = function(time, event, age, sex, ratetable = survexp.ca, alpha = 0.05) {
  n = length(age)
  age = ceiling(age)
  fu.y = floor(time) # duration of fullow-up, years.
  fu.p = time - fu.y   # duration of follow-up, partial year.
  
  # Expected death from population
  E = 0
  can.haz = rep(0, 10)
  names(can.haz)=0:9
  for(i in 1:n) {
    #cat('rate = ', survexp.ca[age[i]+1, sex[i], 1], '\n')
    can.haz = can.haz + survexp.ca[age[i]:(age[i]+9), sex[i], 1]
    #cat(can.haz, '\n')
    lambda = 0
    # total risk, full years
    if(fu.y[i]>0) {
      for (j in 1:fu.y[i])
        lambda = lambda + survexp.ca[age[i]+j, sex[i], 1]      
    }
    
    # total risk, partial years
    lambda = lambda + (survexp.ca[age[i]+fu.y[i]+1, sex[i], 1])*fu.p[i]
    
    # expected number of death 
    E = E + (1-exp(-lambda))
    
    # A better method shall be:
    # E = Total death during followup / Number of alive at age[i]
    #
    #cat('E  = ', E)
  }
  can.haz = can.haz/n
  D = sum(event)
  SMR  = D/E
  SMR.lo = D/E * (1 - 1/9/D - qnorm(1 - alpha/2)/3/sqrt(D))^3
  SMR.up = (D + 1)/E * (1 - 1/9/(D + 1) + qnorm(1 - alpha/2)/3/sqrt(D + 1))^3
  if (E >= 10) {
    chisq = ((abs(D - E) - 0.5)^2)/E
  }
  else {
    Mstar = ifelse(D >= E, D, D + 1)
    chisq = 9 * Mstar * (1 - (1/(9 * Mstar)) - ((E/Mstar)^(1/3)))^2
  }
  p.value = 1 - pchisq(chisq, 1)
  SMR.classic = list(SMR = SMR, SMR.lo = SMR.lo, SMR.up = SMR.up, p.value = p.value)
  
  fit = glm(D ~ 1 + offset(log(E)), family = poisson)
  coefs = summary(fit)$coefficients
  SMR = exp(coefs[1, 1])
  SMR.lo = exp(coefs[1, 1] - qnorm(1 - alpha/2) * coefs[1, 2])
  SMR.up = exp(coefs[1, 1] + qnorm(1 - alpha/2) * coefs[1, 2])
  p.value = coefs[1, 4]
  SMR.poisson = list(SMR = SMR, SMR.lo = SMR.lo, SMR.up = SMR.up, p.value = p.value)
  fit = structure(list(death = D, expected = E, SMR.classic = SMR.classic, SMR.poisson = SMR.poisson, 
              can.haz =can.haz), class = "power.htest")
  return(fit)
}
