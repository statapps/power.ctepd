power.2surv.test = function (n = NULL, sig.level = 0.05, hr = NULL, s0 = 0.7,
          year = 3, accrual = 3, followup = 3, power = NULL, ct.ratio = 1,
          ia.adj = 1.0, alternative = c("two.sided", "one.sided"),
          tol = .Machine$double.eps^0.25) {
  if (sum(sapply(list(n, hr, power), is.null)) != 1) 
    stop("exactly one of 'n', 'hr', 'power', must be NULL")
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.numeric(ia.adj) || ia.adj < 1) 
    stop("'ia.adj' must be numeric and > or = 1")
  
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  #Annual event rate for control arm
  lambda0 = -log(s0)/year
  a = accrual
  f = followup
  td = a + f    #td total duration of the study
  n.event = 0
  za = qnorm(1-sig.level/tside)
  w1 = 1/(1+ct.ratio)
  w0 = 1-w1
  s1 = 0
  
  p.body = quote({
    delta = log(hr)
    lambda1 = lambda0*hr
    s1 = exp(-lambda1*year)
    event_rate0 = 1+(exp(-lambda0*td)-exp(-lambda0*f))/(lambda0*a)
    event_rate1 = 1+(exp(-lambda1*td)-exp(-lambda1*f))/(lambda1*a)
    event_rate = w0*event_rate0 + w1*event_rate1
    #ev0 = floor(w0*n*event_rate0)
    #ev1 = floor(w1*n*event_rate1)
    n.event = n*event_rate 
    1-pnorm(za, -delta/sqrt(1/(w0*w1*n.event)))
  })
  if (is.null(power)) 
    power = eval(p.body)
  else if (is.null(n)) 
    n = uniroot(function(n) eval(p.body) - power, c(2, 1e+07), extendInt = "upX", tol = tol)$root
  else if (is.null(hr)) 
    hr = uniroot(function(hr) eval(p.body) - power, c(0.1, 0.99), tol=tol)$root
  
  n.trt = ceiling(n*w1)
  n.ctl = ceiling(n*w0)
  n = n.ctl + n.trt
  
  n.annual = n/a
  power = eval(p.body)

  NOTE = "Based on the assumption of Exponential distribution"
  METHOD = "Two-sample test power calculation for survival"
  structure(list(n = n, n.trt = n.trt, n.ctl = n.ctl, s0 = s0, s1 = s1, delta = s1-s0, 
                 n.annual= n.annual, n.event = n.event, sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}
