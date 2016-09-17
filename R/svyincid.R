## basic prevalence estimator 
Pest <- function(N_neg, N_rec, N_nr){
  return((N_rec+N_nr)/(N_rec+N_nr+N_neg))
}

## gradient of prevalence estimate
dPest <- function(N_neg, N_rec, N_nr){
  c(N_neg, N_neg, N_rec+N_nr, 0, 0) / (N_neg + N_rec + N_nr)^2
}

## prevalence estimator with untested positive samples (N_q)
Pstar <- function(N_neg, N_rec, N_nr, N_q){
  return(setNames((N_rec+N_nr+N_q)/(N_neg+N_rec+N_nr+N_q), "prev"))
}

## gradient of prevalence estimator with untested positive samples (N_q)
dPstar <- function(N_neg, N_rec, N_nr, N_q){
  return(c(-(N_rec+N_nr+N_q), N_neg, N_neg, N_neg, 0, 0) / (N_neg + N_rec + N_nr + N_q)^2)
}


## basic incidence estimator (ignoring untested positive samples N_q)
Iest <- function(N_neg, N_rec, N_nr, mdri = 188 / 365, frr = 0.016, T = 2){
  est <- (N_rec - frr * (N_rec + N_nr)) / (N_neg * (mdri - frr * T))
  return(setNames(est, "incid"))
}
  
  
## gradient if basic incidence estimator Iest (ignoring N_q)
dIest <- function(N_neg, N_rec, N_nr, mdri = 188 / 365, frr = 0.016, T = 2){
  grad <- c(-(N_rec - frr * (N_rec + N_nr)) / (N_neg^2 * (mdri - frr * T)),
            (1-frr)/(N_neg * (mdri - frr * T)),
            -frr/(N_neg * (mdri - frr * T)),
            - (N_rec - frr * (N_rec + N_nr)) / (N_neg * (mdri - frr * T)^2),
            (N_rec * (T - mdri) - N_nr * mdri) / (N_neg * (mdri - frr * T)^2))
  return(grad)
}

## incidence estimator accouting for untested positive samples
Istar <- function(N_neg, N_rec, N_nr, N_q, mdri = 188 / 365, frr = 0.016, T = 2){
  Qstar <- (N_rec + N_nr + N_q) / (N_rec + N_nr)
  return( setNames(Qstar * Iest(N_neg, N_rec, N_nr, mdri, frr), "incid"))
}

## gradient of incidence estimator accouting for untested positive samples
dIstar <- function(N_neg, N_rec, N_nr, N_q, mdri = 188 / 365, frr = 0.016, T = 2){
  Ival <- Iest(N_neg, N_rec, N_nr, mdri, frr)
  dIval <- dIest(N_neg, N_rec, N_nr, mdri, frr)
  Qstar <- (N_rec + N_nr + N_q) / (N_rec + N_nr)
  dQstar <- c(0,
              - N_q / (N_rec + N_nr)^2,
              - N_q / (N_rec + N_nr)^2,
              1 / (N_rec + N_nr),
              0,
              0)
  return(Qstar * c(dIval[1:3], 0, dIval[4:5]) + dQstar * Ival)
}



## transformation from eta (counts) to {prevalence, incidence}
F <- function(pars, T=2.0){
  ## pars[1] = N_neg
  ## pars[2] = N_rec
  ## pars[3] = N_nr
  ## pars[4] = N_q
  ## pars[5] = mdri
  ## pars[6] = frr

  val <- setNames(c(Pstar(pars[1], pars[2], pars[3], pars[4]),
                    Istar(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], T)),
                  c("prev", "incid"))
                    
  return(val)
}

dF <- function(pars, T=2.0){
  cbind(prev=dPstar(pars[1], pars[2], pars[3], pars[4]),
        incid=dIstar(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], T))
}


## transformation of {prev, incid} to {probit(prev), log(incid)}
G <- function(theta){
  ## theta[1] = prevalence
  ## theta[2] = incidence rate
  ## G(theta) = c(probit(prevalence), log(incidence))

  return(c(qnorm(theta[1]), log(theta[2])))
}

dG <- function(theta){
  ## theta[1] = prevalence
  ## theta[2] = incidence rate

  ## g1(x) = qnorm(x)
  ## dg1(x) = 1/dnorm(qnorm(x))  since d/dx finv(x) = 1/f'(finv(x))

  ## g2(x) = log(x)
  ## dg2(x) = 1/x

  return(diag(c(1/dnorm(qnorm(theta[1])), 1/theta[2])))
}



svyprevincid <- function(recent, design, mdri=130/365, frr=0.012,
                         se.mdri=(142-118)/365/(2*qnorm(0.975)), se.frr=(0.025 - 0.0001)/(2*qnorm(0.975)), deff=FALSE){
  tot <- svytotal(recent, design)
  eta <- c(coef(tot), mdri, frr)
  Sigma <- diag(c(rep(0, 4), se.mdri^2, se.frr^2))
  Sigma[1:4,1:4] <- vcov(tot)

  incprev <- F(eta)
  class(incprev) <- "svystat"
  attr(incprev, "var") <- t(dF(eta)) %*% Sigma %*% dF(eta)
  attr(incprev, "statistic") <- "estimate"

  if(is.character(deff) || deff){
    nobs <- sum(weights(design) != 0)
    tot.srs <- nobs*tot/sum(tot)
    eta.srs <- eta
    eta.srs[1:4] <- tot.srs
    Sigma.srs <- Sigma
    Sigma.srs[1:4, 1:4] <- diag(tot.srs) - tot.srs %o% tot.srs / sum(tot.srs) # multinomial variance for simple random sample
    vsrs <-  t(dF(eta.srs)) %*% Sigma.srs %*% dF(eta.srs)
    attr(incprev, "deff") <- attr(incprev, "var") / vsrs
  }
  
  return(incprev)
}

#@' param object
eppformat <- function(object){
  rval <- c(rbind(coef(object), SE(object)), cov2cor(vcov(object))[2,1])
  rval <- data.frame(as.list(rval))
  names(rval) <- c("prev", "prev.se", "incid", "incid.se", "corr")
  return(rval)
}
