power.fdr.tt = function(n1=NULL, cc.ratio=1, pi0=0.99, fdr=0.05, delta=NULL, power = NULL, 
             sigma = 1, tol = .Machine$double.eps^0.25) {
  
  if (sum(sapply(list(n1, delta, power), is.null)) != 1) 
    stop("exactly one of 'n1', 'delta', 'power', must be NULL")
  
  Lambda = fdr/(1-fdr)*(1-pi0)/pi0
  cg = 1.96
  n0 = n1
  
  p.body = quote({
    n0 = floor(n1*cc.ratio)
    n  = n1 + n0 - 2
    theta = delta/(sigma*sqrt(1/n1+1/n0))
    cg = uniroot(function(x) {
      supressWarings(2*pt(-x, n)/(1-pt(x, n, theta)+pt(-x, n, theta))) - Lambda
    }, c(0.0000001, 100), tol = 0.001)$root
    suppressWarnings(1-pt(cg, n, theta)+pt(-cg, n, theta))
  })
  
  if(is.null(power))
    power = eval(p.body)
  else if (is.null(n1)) {
    n1 = 10
    p0 = try(eval(p.body))
    while(p0 < power | class(p0)=="try-error") {
      n1 = n1*2
      p0 = try(eval(p.body))
    }
    while(p0 > power){
      n1 = n1-1
      p0 = eval(p.body)
    }
    n1 = n1+1
  }
  else if (is.null(delta))
    delta = uniroot(function(delta) eval(p.body) - power, c(0.001, 10), 
                 extendInt="upX", tol = tol)$root
  
  power = eval(p.body)
  
  NOTE = 'Based on two sample t-tst'
  METHOD = 'Power for Micro-array with FDR'
  
  structure(list(case.control = c(n1, n0), n.total = n0+n1, fdr = fdr, pi0 = pi0, c.reject = cg, 
                 Lambda = Lambda, power = power, delta = delta, sigma = sigma,
                 note = NOTE, method = METHOD), class = "power.htest")

}
