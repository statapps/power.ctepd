power.surv.test=function(n = NULL, sig.level = 0.05, s0 = NULL, s1 = NULL, 
                         year = 3, followup = 4,
                         power = NULL, tol = .Machine$double.eps^0.25){
  if (sum(sapply(list(n, s1, power), is.null)) != 1) 
    stop("exactly one of 'n', 's1', 'power' must be NULL")

  lambda0 = -log(s0)/year
  p.body <- quote({
    lambda1 = -log(s1)/year
    person.year = n*(1-exp(-lambda0*followup))/lambda0
    person.year1 = n*(1-exp(-lambda1*followup))/lambda1
    
    k = qpois(sig.level, person.year*lambda0)
    alpha = ppois(k, person.year*lambda0)
    while(k>0 & alpha > sig.level){
      k = k - 1
      alpha =  ppois(k, person.year*lambda0)
    }
    sig.level = alpha
      
    pw = ppois(k, person.year1*lambda1)
    ppois(k, n*lambda1*year)
  })

  if(is.null(power))
    power = eval(p.body)
  else if (is.null(n))
    n = uniroot(function(n) eval(p.body) - power, c(25+1e-10, 1e+07), 
                extendInt="upX", tol = tol)$root
  else if (is.null(s1))
    s1 = uniroot(function(s1) eval(p.body) - power, c(s0, 0.98))$root
  n = ceiling(n)
  power = eval(p.body)
  
  alternative = 'less'
  NOTE = 'Based on the assumption of Poisson Process'
  METHOD = 'One sample test for survival distribution'

  structure(list(n = n, s0 = s0, s1 = s1, person.year = person.year, 
     crtl.value = k, sig.level = sig.level, 
     power = power, alternative = alternative, note = NOTE, 
     method = METHOD), class = "power.htest")
}

