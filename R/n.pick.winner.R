##########################################################################################
# Pick winner design by Richard Simon, Robert E. Wittes and Susan S. Ellenberg 
# Randomized phase II clinical Trials. Cancer Treat Rep 69: 1375-1381, 1985.
##########################################################################################
n.pick.winner = function(n = NULL, k = 2, p0=NULL, p1=NULL, power=NULL, 
                         tol = .Machine$double.eps^0.25) {
  if (sum(sapply(list(n, p1, power), is.null)) != 1) 
    stop("exactly one of 'n', 'p1', 'power', must be NULL")
  if(k<2) stop('Number of arms must > or = 2.\n')
  if(!is.null(p1)) {
    if(p1 >=1) stop('p1 must less than 1')
    if(p0>=p1) stop('p0 must less than p1.\n')
  }
  if(p0 <=0) stop('p0 must greater than 0')
  # probability of picking the winner arm
  p.body = quote({
    n = ceiling(n)
    x = 0:n
    #CDF: cumulativ prob of have x or less responses
    Bi0 = pbinom(x, n, p0)
    Bi1 = pbinom(x, n, p1)
    # PDF: prob of have exactly x responses
    bi0 = dbinom(x, n, p0)
    bi1 = dbinom(x, n, p1)
  
    # Prob that max respose is i in (k-1) null arms
    Bix = c(0, Bi0)[1:(n+1)]
    fx = Bi0^(k-1) - Bix^(k-1)
    
    ## For ties, best arm was selected when being tied with 1 or more null arms.
    gj = 0
    for (j in 1:(k-1)) {
      ckj = choose(k-1, j)
      gj = gj+ckj*bi0^j*Bix^(k-1-j)/(j+1)
    }  
    pb = sum(fx*(1-Bi1)) + sum(bi1*gj)
  })
  
  if(is.null(n)) {
    n = uniroot(function(n) eval(p.body)-power, c(5, 500), tol = tol)$root
    n = n + 1
  }
  if(is.null(p1))
    p1 = uniroot(function(p1) eval(p.body)-power, c(p0-0.001, 0.99), tol = tol)$root
  power = eval(p.body)
  
  NOTE = "Based on method by Richard Simon (1985)."
  METHOD = "Sample size for pick winner design"
  
  structure(list(n = n, arms = k, n.total = n*k, p0 = p0, p1 = p1, power = power,
                 note = NOTE, method = METHOD), class = "power.htest")
}
