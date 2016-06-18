power.cc.hap = function(n1=NULL, cc.ratio=1, p0=0.01, dhf=0.05, n.hap=2, 
                        rr=NULL, rr2=NULL, sig.level=0.05, power = NULL, 
                        genModel=c('co.dom', 'co.dom2', 'dom', 'rec'), 
                        tol = .Machine$double.eps^0.25) {
  if (sum(sapply(list(n1, power, rr), is.null)) != 1) 
    stop("exactly one of 'n', 'power', 'rr', must be NULL")
  
  genModel = match.arg(genModel)
  METHOD = switch(genModel, co.dom = 'Co-dominant model (1 degree of freedom)', 
                  co.dom2 = 'Co-dominant model (2 degrees of freedom)',
                  dom = 'Dominant model', rec = 'Recessive model')
  f1 = dhf
  J = n.hap - 1
  
  ## prob of 0, 1 or 2 of disease haplogype
  Fd = c((1-f1)^2, 2*f1*(1-f1), f1^2)
  n0 = NULL
  Penetrance = NULL
  p.body = quote({
    n0 = n1*cc.ratio
    if(is.null(rr2))  rr2 = rr*rr
    
    RR = switch(genModel, co.dom = c(1, rr, rr^2), co.dom2 = c(1, rr, rr2),
                dom = c(1, rr, rr), rec = c(1, 1, rr))

    #Baseline disease probability = r0
    r0 = p0/sum(Fd*RR)     # p0 is the same as phi in the manuscript
    Penetrance = r0*RR     # Population penetrance
    
    a = r0*(RR[2]*f1+RR[1]*(1-f1))/p0
    b = ((1-r0*RR[2])*f1+(1-r0*RR[1])*(1-f1))/(1-p0)
    d = (n0*a^2+n1*b^2)/(n0*a+n1*b)   #d is the same as c in the manuscript

    delta = n0*n1*(a-b)^2/(n0*a+n1*b) * ((1-f1)+d*(1-f1)^2/(1-d*(1-f1)))
    x2 = qchisq(1-sig.level, J)
    1 - pchisq(x2, J, ncp = delta)
  })
  
  if(is.null(n1)){
    n1 = uniroot(function(n1) eval(p.body)-power, c(2,1e+06), tol=tol)$root
    n1 = ceiling(n1)
  }
  if(is.null(rr)) 
    rr = uniroot(function(rr) eval(p.body)-power, c(1.01, 1000), tol=tol)$root
  
  power = eval(p.body)
  NOTE = 'Power for haplotype case control study'

  structure(list(case.control = c(n1, n0), n.total = n0+n1, sig.level = sig.level, 
                 power = power, rr = rr, rr2 = rr2, Penetrance=Penetrance,
                 note = NOTE, method = METHOD), class = "power.htest")
}
