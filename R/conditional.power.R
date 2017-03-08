conditional.power = function(tk = 0.0, theta=0.5, sig.level = 0.025, power = 0.8, tm = 0.5) {
  za = -qnorm(sig.level)
  zb =  qnorm(power)
  Im = ((za+zb)/theta)^2        ## Im for full information
  Ik = Im * tm
  zk = tk * sqrt(Ik)

  qtl = -zk*sqrt(Ik) - za*sqrt(Im)-theta*(Im-Ik)
  qtl = qtl/sqrt(Im-Ik)
  cpl = pnorm(qtl)*100
  cat('\n Lower conditional power = ', cpl, '\n')

  qtu = zk*sqrt(Ik) - za*sqrt(Im)+theta*(Im-Ik)
  qtu = qtl/sqrt(Im-Ik)
  cpu = pnorm(qtu)*100
  cat('\n Upper conditional power = ', cpu, '\n')
  NOTE = 'c.power.low is lower conditional power for one-side futility analysis. \n      c.power is for two-sided test.'
  METHOD = 'Conditional Power for sequencial tests'

  cat('\n Two-sided conditional power = ', (cpl+cpu), '\n')
  fit = structure(list(c.power.low = cpl, c.power.up = cpu, c.power = cpl + cpu,
            zk = zk, I.max = Im, 
            sig.level = sig.level, power = power, time = tm,
            note = NOTE, method = METHOD), 
            class = "power.htest")
  return(fit)
}

c.power.surv = function(hrk = 1.0, HR = 0.7, sig.level = 0.025, power = 0.8, tm = 0.5) {
  theta = log(HR)
  tk    = log(hrk)
  fit = conditional.power(tk = tk, theta = theta, sig.level = sig.level, power = power, tm = tm)
  fit$HR.k = hrk
  fit$HR   = HR
  fit$events = ceiling(fit$I.max*4)

  return(fit)
}
