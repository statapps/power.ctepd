conditional.power = function(tk = 0.0, theta=-0.5, sig.level = 0.025, power = 0.8, tm = 0.5) {
  if(theta > 0) {
    cat("\n theta > 0 is transformed to -theta for one-sided test H1: theta=delta (<0)\n")
    tk    = tk - theta
    theta =  0 - theta
  }
  za = -qnorm(sig.level)
  zb =  qnorm(power)
  Im = ((za+zb)/theta)^2        ## Im for full information
  Ik = Im * tm
  zk = tk * sqrt(Ik)

  qtl = -zk*sqrt(Ik) - za*sqrt(Im)-theta*(Im-Ik)
  qtl = qtl/sqrt(Im-Ik)
  cpl = pnorm(qtl)*100
  cat('\n One-sided conditional power at H1: theta=delta(<0):', cpl, '\n')

  qtu = zk*sqrt(Ik) - za*sqrt(Im)+theta*(Im-Ik)
  qtu = qtl/sqrt(Im-Ik)
  cpu = pnorm(qtu)*100
  #cat('\n One-sided conditional power = ', cpu, '\n')
  NOTE = 'cpower is conditional power for one-side futility analysis.' #\n      cpower is for two-sided test.'
  METHOD = 'Conditional Power for sequential tests'

  #cat('\n Two-sided conditional power = ', (cpl+cpu), '\n')
  fit = structure(list(cpower = cpl, #cpower.up = cpu, cpower = cpl + cpu,
            zk = zk, theta = theta, I.max = Im, 
            sig.level = sig.level, power = power, time = tm,
            note = NOTE, method = METHOD), 
            class = "power.htest")
  return(fit)
}

cpower.surv = function(hrk = 1.0, HR = 0.7, sig.level = 0.025, power = 0.8, tm = 0.5) {
  #if(HR > 1) stop("This program only works for superiority trial with HR < 1\n")
  theta = log(HR)
  tk    = log(hrk)
  fit = conditional.power(tk = tk, theta = theta, sig.level = sig.level, power = power, tm = tm)
  fit$HR.k = hrk
  fit$HR   = HR
  fit$events = ceiling(fit$I.max*4)

  return(fit)
}

cpower.tt = function(muk = 0.0, mu = 0.7, sig.level = 0.025, power = 0.8, tm = 0.5, 
                    sigma = 1.0, cc.ratio = 1) {
  theta = mu/sigma
  tk    = muk/sigma

  za = -qnorm(sig.level)
  zb =  qnorm(power)
  n.total = ((za+zb)/theta)^2 *(2 + cc.ratio + 1/cc.ratio)
  n1 = ceiling(n.total*cc.ratio/(1+cc.ratio))
  n0 = ceiling(n.total/(1+cc.ratio))
  n.total = n1 + n0

  fit = conditional.power(tk = tk, theta = theta, sig.level = sig.level, power = power, tm = tm)
  fit$trt.ctl = c(n1, n0)
  fit$n.total = n.total

  return(fit)
}
