conditional.power = function(tk = 1.3, theta=0.5, sig.level = 0.025, power = 0.8, tm = 0.5) {
  za = -qnorm(sig.level)
  zb =  qnorm(power)
  Im = ((za+zb)/theta)^2  ## Im for full information
  Ik = Im * tm
  zk = Tk/sqrt(Ik)

  qtl = -zk*sqrt(Ik) - za*sqrt(Im)-theta*(Im-Ik)
  qtl = qtl/sqrt(Im-Ik)
  cpl = pnorm(qtl)
  cat('\n Lower conditional power = ', cpl*100, '\n')

  qtu = zk*sqrt(Ik) - za*sqrt(Im)+theta*(Im-Ik)
  qtu = qtl/sqrt(Im-Ik)
  cpu = pnorm(qtu)
  cat('\n Upper conditional power = ', cpu*100, '\n')

  cat('\n Two-sided conditional power = ', (cpl+cpu)*100, '\n')

  return(c(cpu, cpl))
}

conditional.power.surv = function(hrk = 1.1, HR = 0.7, sig.level = 0.025, power = 0.8, tm = 0.5) {
  theta = log(HR)
  tk    = log(hrk)
  cp = conditional.power(tk = tk, theta = theta, sig.level = sig.level, power = power, tm = tm)
  return(cp)
}
