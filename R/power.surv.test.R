power.surv.test=function(n = NULL, sig.level = 0.05, s0 = NULL, s1 = NULL, 
                         year = 3, followup = 3, futile = NULL, futile.prob = 0.01,
                         power = NULL, tol = .Machine$double.eps^0.25){
  if (sum(sapply(list(n, s1, power), is.null)) != 1) 
    stop("exactly one of 'n', 's1', 'power', must be NULL")

  lambda0 = -log(s0)/year

  # initial output value for person.year and critical value k
  person.year = 0
  k = 0
  p.body = quote({
    lambda1 = -log(s1)/year
    person.year = n*(1-exp(-lambda0*followup))/lambda0
    
    k = qpois(sig.level, person.year*lambda0)
    alpha = ppois(k, person.year*lambda0)
    while(k>0 & alpha > sig.level){
      k = k - 1
      alpha =  ppois(k, person.year*lambda0)
    }
    sig.level = alpha

    # lambda*t = n*person.year1*lambda1  
    ppois(k, n*(1-exp(-lambda1*followup)))
  })

  if(is.null(power))
    power = eval(p.body)
  else if (is.null(n))
    n = uniroot(function(n) eval(p.body) - power, c(25+1e-10, 1e+07), 
                extendInt="upX", tol = tol)$root
  else if (is.null(s1))
    s1 = uniroot(function(s1) eval(p.body) - power, c(s0, 0.99))$root

  n = ceiling(n)
  power = eval(p.body)
  lambda = k/person.year
  crtl.survf= exp(-lambda*year)

  alternative = 'greater'
  NOTE = 'Based on the assumption of Poisson Process'
  METHOD = 'One sample test for survival distribution'

  fit = structure(list(n = n, s0 = s0, s1 = s1, person.year = person.year, 
     crtl.event = k, crtl.survf = crtl.survf, sig.level = sig.level, 
     power = power, alternative = alternative, 
     note = NOTE, method = METHOD), class = "power.htest")

  if(is.numeric(futile)) {
    lambda1 = -log(s1)/year
    if (futile > 1 | futile < 0) stop("futile shall be between 0 and 1")
    futile.py = person.year*futile
    cat(lambda1, lambda1*futile.py)
    futile.event = qpois(1-futile.prob, lambda1*futile.py)

    fit = structure(list(n = n, s0 = s0, s1 = s1, person.year = person.year,
    crtl.event = k, crtl.survf = crtl.survf, sig.level = sig.level, power = power, 
    futile.pyear = futile.py, futile.event = futile.event, futile.prob = futile.prob, 
    alternative = alternative, note = NOTE, method = METHOD), class = "power.htest")

  }


  return(fit)
}

