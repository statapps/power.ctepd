power.cc.snp = function(n1=NULL, cc.ratio=1, p0=0.01, maf=0.05, daf=0.05, Dp=NULL, 
                        R2=NULL, rr=NULL, rr2=NULL, sig.level=0.05, power = NULL, 
                        edf=1, genModel=c('co.dom', 'co.dom2', 'dom', 'rec'), 
                        tol = .Machine$double.eps^0.25) {
#  if (sum(sapply(list(n1, power, rr), is.null)) != 1) 
#    stop("exactly one of 'n', 'power', 'rr', must be NULL")
  genModel = match.arg(genModel)
  pd = daf
  p1 = maf
  
  if(is.null(Dp) & is.null(R2)) R2 = 1
  if(is.null(Dp)) {
    R2.max = min(p1*(1-pd), pd*(1-p1))^2/(p1*(1-p1)*pd*(1-pd))
    D = sqrt(R2*R2.max*p1*(1-p1)*pd*(1-pd))
  } else {
    D.max = min(p1*(1-pd), pd*(1-p1));
    D = Dp*D.max;
  }
  
  # haplotype of marker allele and disease allele
  hw1 = (1-pd)*p1     - D  # +1 haplotype
  hw2 = (1-pd)*(1-p1) + D  # +2 haplotype
  hd1 = pd*p1         + D  # d1 haplotype
  hd2 = pd*(1-p1)     - D  # d2 haplotype
  
  A = matrix(c(hw2^2,     2*hw2*hd2,             hd2^2,      # 0 marker
               2*hw1*hw2, 2*(hw1*hd2 + hw2*hd1), 2*hd1*hd2,  # 1 marker
               hw1^2,     2*hw1*hd1,             hd1^2),     # 2 marker
             byrow = TRUE, nrow = 3)
  #print(A)
  
  alpha = sig.level/edf
  Fd = c((1-pd)^2, 2*pd*(1-pd), pd^2)

  p.body = quote({
    if(is.null(rr2))  rr2 = rr*rr
    RR = switch(genModel, co.dom = c(1, rr, rr^2), co.dom2 = c(1, rr, rr2),
                dom = c(1, rr, rr), rec = c(1, 1, rr))
                                     
    #Baseline disease probability
    f0 = p0/sum(Fd*RR)
    Penetrance = f0*RR     # Population penetrance
    
    # p_i1 = Pr(X=i|D=1), p_i0 = Pr(X=i|D=0)
    # P0 = (1-Penetrance).*Fd./(1-phi)
    # P1 = Penetrance.*Fd./phi
    
    P0 = t(A%*%(1-Penetrance))/(1-p0)
    P1 = t(A%*%(Penetrance))/p0
    
    # expected disease allele frequency in case and control groups
    df0 = P0[2]/2+P0[3]
    df1 = P1[2]/2+P1[3]
    if(genModel == 'co.dom') {
      fullname = 'Co-dominant model (1 degree of freedom)'
      n0 = n1*cc.ratio
      N = rbind(n1*P1, n0*P0)
      Y = c(0, 1, 2)   # score
      A = N[1, ]       # with disease
      B = N[2, ]       # free of disease
      M = A + B        # column total
      ay = A%*%Y
      my = M%*%Y
      my2 = M%*%(Y^2)
      T = n1 + n0
      N1T = n1/T
      chi2 = sum((ay - N1T * my)^2)
      V = n0*N1T *(my2 - (my^2)/T)
      delta = (T-1)*chi2/V
      x = qchisq(1-alpha, 1)
      1 - pchisq(x, 1, ncp=delta)
    } else if (genModel == 'co.dom2') {
      fullname = 'Co-dominant model (2 degrees of freedom)'
      n0 = n1*cc.ratio
      N = rbind(n1*P1, n0*P0)
      dp2 = (P0-P1)^2
      #####Column sum or total sum: Matlab NP = sum(N)
      NP = apply(N, 2, sum)
      delta = n1*n0*sum(dp2/NP)
      x = qchisq(1-alpha, 2)
      1 - pchisq(x, 2, ncp=delta)
    } else if(genModel=='dom'){
      fullname = 'Dominant model'
      rr2 = RR[3]
      n0 = n1*cc.ratio
      N = rbind(n1*P1, n0*P0)
      x = abs(qnorm(alpha/2, 0, 1))
      a = N[1, 1]; b = sum(N[1, 2:3])
      c = N[2, 1]; d = sum(N[2, 2:3])
      sc = log(a)+log(d)-log(b)-log(c)
      sd = sqrt(1/a+1/b+1/c+1/d)
      z = abs(sc)/sd
      1-pnorm(x, z, 1) + pnorm(-x, z, 1);
    } else if(genModel=='rec'){
      fullname = 'Recessive model'
      n0 = n1*cc.ratio
      N = rbind(n1*P1, n0*P0)
      x = abs(qnorm(alpha/2, 0, 1))
      a = sum(N[1, 1:2]); b = N[1, 3]
      c = sum(N[2, 1:2]); d = N[2, 3]
      sc = log(a)+log(d)-log(b)-log(c)
      sd = sqrt(1/a+1/b+1/c+1/d)
      z = abs(sc)/sd
      1-pnorm(x, z, 1) + pnorm(-x, z, 1);
    }
  })  
  
  
  if(is.null(n1)) 
    n1 = uniroot(function(n1) eval(p.body)-power, c(2,1e+06), tol=tol)$root
  n1 = ceiling(n1)
  if(is.null(rr)) 
    rr = uniroot(function(rr) eval(p.body)-power, c(1.01, 1000), tol=tol)$root
  
  power = eval(p.body)
  NOTE = 'Power for genetics case control studies'
  #print(power)
  structure(list(case.control = c(n1, n0), n.total = n0+n1, sig.level = sig.level, 
                 sig.adjust = alpha, power = power, rr = rr, rr2 = rr2, D = D, 
                 note = NOTE, method = fullname), class = "power.htest")
}
