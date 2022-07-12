#' Get asymptotic frequencies from observed genotype frequencies.
#'
#' @param yww The observed genotype frequencies.
#' @param niter The number of EM iterations.
#' @param alpha The double reduction rate (assumed known).
#'
#' @return The genotype frequencies at equilibrium.
#'
#' @author Jing Wang and David Gerard
#'
#' @references
#' \itemize{
#'   \item{J. Wang, L. Feng, S. Mu, A. Dong, J. Gan, Z. Wen, J. Meng, M. Li, R. Wu, and L. Sun. Asymptotic tests for Hardy-Weinberg equilibrium in hexaploids. \emph{Horticulture Research}, 9, 05 2022. \doi{10.1093/hr/uhac104}.}
#' }
#'
#' @examples
#' yww <- c(29, 21, 17, 10, 10, 10, 23) / 120
#' hex_recursive(yww = yww, niter = 8, alpha = 0)
#'
#' @export
hex_recursive <- function(yww, niter = 8, alpha = 0) {
  times <- seq(0,niter,1)
  N <- length(times)

  m <- numeric(N)
  n <- numeric(N)
  k <- numeric(N)
  e <- numeric(N)
  h <- numeric(N)
  f <- numeric(N)
  u <- numeric(N)

  m[1] <- as.numeric(yww[1])
  n[1] <- as.numeric(yww[2])
  k[1] <- as.numeric(yww[3])
  e[1] <- as.numeric(yww[4])
  h[1] <- as.numeric(yww[5])
  f[1] <- as.numeric(yww[6])
  u[1] <- as.numeric(yww[7])

  freq <- c()
  for(t in 1:(N-1)){

    m[t+1] <- (m[t]^2*1+2*m[t]*n[t]*(1/2+1/6*alpha)+2*m[t]*k[t]*(1/5+1/5*alpha)+2*m[t]*e[t]*(1/20+3/20*alpha)+2*m[t]*h[t]*(1/15*alpha)+2*m[t]*f[t]*0+2*m[t]*u[t]*0
               +n[t]^2*(1/4+1/36*alpha^2+1/6*alpha)+2*n[t]*k[t]*(1/10+1/30*alpha^2+2/15*alpha)+2*n[t]*e[t]*(1/40+1/40*alpha^2+1/12*alpha)+2*n[t]*h[t]*(1/90*alpha^2+1/30*alpha)+2*n[t]*f[t]*0+2*n[t]*u[t]*0
               +k[t]^2*(1/25+1/25*alpha^2+2/25*alpha)+2*k[t]*e[t]*(1/100+3/100*alpha^2+1/25*alpha)+2*k[t]*h[t]*(1/75*alpha^2+1/75*alpha)+2*k[t]*f[t]*0+2*k[t]*u[t]*0
               +e[t]^2*(1/400+9/400*alpha^2+3/200*alpha)+2*e[t]*h[t]*(1/100*alpha^2+1/300*alpha)+2*e[t]*f[t]*0+2*e[t]*u[t]*0
               +h[t]^2*(1/225*alpha^2)+2*h[t]*f[t]*0+2*h[t]*u[t]*0
               +f[t]^2*0+2*f[t]*u[t]*0
               +u[t]^2*0)

    n[t+1] <- (m[t]^2*0+2*m[t]*n[t]*(1/2-1/3*alpha)+2*m[t]*k[t]*(3/5-1/3*alpha)+2*m[t]*e[t]*(9/20-3/20*alpha)+2*m[t]*h[t]*(1/5+1/15*alpha)+2*m[t]*f[t]*(1/6*alpha)+2*m[t]*u[t]*0
               +n[t]^2*(1/2-1/9*alpha^2-1/6*alpha)+2*n[t]*k[t]*(2/5-11/90*alpha^2-1/30*alpha)+2*n[t]*e[t]*(1/4-3/40*alpha^2+7/120*alpha)+2*n[t]*h[t]*(1/10-1/90*alpha^2+1/10*alpha)+2*n[t]*f[t]*(1/36*alpha^2+1/12*alpha)+2*n[t]*u[t]*0
               +k[t]^2*(6/25-2/15*alpha^2+8/75*alpha)+2*k[t]*e[t]*(3/25-2/25*alpha^2+2/15*alpha)+2*k[t]*h[t]*(1/25-2/225*alpha^2+7/75*alpha)+2*k[t]*f[t]*(1/30*alpha^2+1/30*alpha)+2*k[t]*u[t]*0
               +e[t]^2*(9/200-9/200*alpha^2+3/25*alpha)+2*e[t]*h[t]*(1/100+19/300*alpha)+2*e[t]*f[t]*(1/40*alpha^2+1/120*alpha)+2*e[t]*u[t]*0
               +h[t]^2*(2/225*alpha^2+2/75*alpha)+2*h[t]*f[t]*(1/90*alpha^2)+2*h[t]*u[t]*0
               +f[t]^2*0+2*f[t]*u[t]*0
               +u[t]^2*0)

    k[t+1] <- (m[t]^2*0+2*m[t]*n[t]*(1/6*alpha)+2*m[t]*k[t]*(1/5+1/15*alpha)+2*m[t]*e[t]*(9/20-3/20*alpha)+2*m[t]*h[t]*(3/5-1/3*alpha)+2*m[t]*f[t]*(1/2-1/3*alpha)+2*m[t]*u[t]*0
               +n[t]^2*(1/4+1/6*alpha^2-1/6*alpha)+2*n[t]*k[t]*(2/5+7/45*alpha^2-4/15*alpha)+2*n[t]*e[t]*(9/20+1/20*alpha^2-13/60*alpha)+2*n[t]*h[t]*(2/5-1/15*alpha^2-1/10*alpha)+2*n[t]*f[t]*(1/4-1/9*alpha^2)+2*n[t]*u[t]*0
               +k[t]^2*(11/25+31/225*alpha^2-22/75*alpha)+2*k[t]*e[t]*(37/100+3/100*alpha^2-11/75*alpha)+2*k[t]*h[t]*(6/25-19/225*alpha^2+1/25*alpha)+2*k[t]*f[t]*(1/10-11/90*alpha^2+2/15*alpha)+2*k[t]*u[t]*0
               +e[t]^2*(99/400-9/400*alpha^2-3/200*alpha)+2*e[t]*h[t]*(3/25-7/100*alpha^2+31/300*alpha)+2*e[t]*f[t]*(1/40-3/40*alpha^2+2/15*alpha)+2*e[t]*u[t]*0
               +h[t]^2*(1/25-9/225*alpha^2+8/75*alpha)+2*h[t]*f[t]*((-1/90)*alpha^2+1/15*alpha)+2*h[t]*u[t]*0
               +f[t]^2*(1/36*alpha^2)+2*f[t]*u[t]*0
               +u[t]^2*0)

    e[t+1] <- (m[t]^2*0+2*m[t]*n[t]*0+2*m[t]*k[t]*(1/15*alpha)+2*m[t]*e[t]*(1/20+3/20*alpha)+2*m[t]*h[t]*(1/5+1/5*alpha)+2*m[t]*f[t]*(1/2+1/6*alpha)+2*m[t]*u[t]*1
               +n[t]^2*((-1/9)*alpha^2+1/6*alpha)+2*n[t]*k[t]*(1/10-1/15*alpha^2+1/10*alpha)+2*n[t]*e[t]*(1/4+1/20*alpha^2-1/15*alpha)+2*n[t]*h[t]*(2/5+7/45*alpha^2-1/5*alpha)+2*n[t]*f[t]*(1/2+1/6*alpha^2-1/6*alpha)+2*n[t]*u[t]*(1/2+1/6*alpha)
               +k[t]^2*(6/25-4/225*alpha^2-2/75*alpha)+2*k[t]*e[t]*(37/100+2/25*alpha^2-59/300*alpha)+2*k[t]*h[t]*(11/25+4/25*alpha^2-22/75*alpha)+2*k[t]*f[t]*(2/5+7/45*alpha^2-1/5*alpha)+2*k[t]*u[t]*(1/5+1/5*alpha)
               +e[t]^2*(41/100+9/100*alpha^2-6/25*alpha)+2*e[t]*h[t]*(37/100+2/25*alpha^2-59/300*alpha)+2*e[t]*f[t]*(1/4+1/20*alpha^2-1/15*alpha)+2*e[t]*u[t]*(1/20+3/20*alpha)
               +h[t]^2*(6/25-4/225*alpha^2-2/75*alpha)+2*h[t]*f[t]*(1/10-1/15*alpha^2+1/10*alpha)+2*h[t]*u[t]*(1/15*alpha)
               +f[t]^2*((-1/9)*alpha^2+1/6*alpha)+2*f[t]*u[t]*0
               +u[t]^2*0)

    h[t+1] <- (m[t]^2*0+2*m[t]*n[t]*0+2*m[t]*k[t]*0+2*m[t]*e[t]*0+2*m[t]*h[t]*0+2*m[t]*f[t]*0+2*m[t]*u[t]*0
               +n[t]^2*(1/36*alpha^2)+2*n[t]*k[t]*((-1/90)*alpha^2+1/15*alpha)+2*n[t]*e[t]*(1/40-3/40*alpha^2+2/15*alpha)+2*n[t]*h[t]*(1/10-11/90*alpha^2+2/15*alpha)+2*n[t]*f[t]*(1/4-1/9*alpha^2)+2*n[t]*u[t]*(1/2-1/3*alpha)
               +k[t]^2*(1/25-9/225*alpha^2+8/75*alpha)+2*k[t]*e[t]*(3/25-7/100*alpha^2+31/300*alpha)+2*k[t]*h[t]*(6/25-19/225*alpha^2+1/25*alpha)+2*k[t]*f[t]*(2/5-1/15*alpha^2-1/10*alpha)+2*k[t]*u[t]*(3/5-1/3*alpha)
               +e[t]^2*(99/400-9/400*alpha^2-3/200*alpha)+2*e[t]*h[t]*(37/100+3/100*alpha^2-11/75*alpha)+2*e[t]*f[t]*(9/20+1/20*alpha^2-13/60*alpha)+2*e[t]*u[t]*(9/20-3/20*alpha)
               +h[t]^2*(11/25+31/225*alpha^2-22/75*alpha)+2*h[t]*f[t]*(2/5+7/45*alpha^2-4/15*alpha)+2*h[t]*u[t]*(1/5+1/15*alpha)
               +f[t]^2*(1/4+1/6*alpha^2-1/6*alpha)+2*f[t]*u[t]*(1/6*alpha)
               +u[t]^2*0)

    f[t+1] <- (m[t]^2*0+2*m[t]*n[t]*0+2*m[t]*k[t]*0+2*m[t]*e[t]*0+2*m[t]*h[t]*0+2*m[t]*f[t]*0+2*m[t]*u[t]*0
               +n[t]^2*0+2*n[t]*k[t]*(1/90*alpha^2)+2*n[t]*e[t]*(1/40*alpha^2+1/120*alpha)+2*n[t]*h[t]*(1/30*alpha^2+1/30*alpha)+2*n[t]*f[t]*(1/36*alpha^2+1/12*alpha)+2*n[t]*u[t]*(1/6*alpha)
               +k[t]^2*(2/225*alpha^2+2/75*alpha)+2*k[t]*e[t]*(1/100+19/300*alpha)+2*k[t]*h[t]*(1/25-2/225*alpha^2+7/75*alpha)+2*k[t]*f[t]*(1/10-1/90*alpha^2+1/10*alpha)+2*k[t]*u[t]*(1/5+1/15*alpha)
               +e[t]^2*(9/200-9/200*alpha^2+3/25*alpha)+2*e[t]*h[t]*(3/25-2/25*alpha^2+2/15*alpha)+2*e[t]*f[t]*(1/4-3/40*alpha^2+7/120*alpha)+2*e[t]*u[t]*(9/20-3/20*alpha)
               +h[t]^2*(6/25-2/15*alpha^2+8/75*alpha)+2*h[t]*f[t]*(2/5-11/90*alpha^2-1/30*alpha)+2*h[t]*u[t]*(3/5-1/3*alpha)
               +f[t]^2*(1/2-1/9*alpha^2-1/6*alpha)+2*f[t]*u[t]*(1/2-1/3*alpha)
               +u[t]^2*0)

    u[t+1] <- (m[t]^2*0+2*m[t]*n[t]*0+2*m[t]*k[t]*0+2*m[t]*e[t]*0+2*m[t]*h[t]*0+2*m[t]*f[t]*0+2*m[t]*u[t]*0
               +n[t]^2*0+2*n[t]*k[t]*0+2*n[t]*e[t]*0+2*n[t]*h[t]*0+2*n[t]*f[t]*0+2*n[t]*u[t]*0
               +k[t]^2*(1/225*alpha^2)+2*k[t]*e[t]*(1/100*alpha^2+1/300*alpha)+2*k[t]*h[t]*(1/75*alpha^2+1/75*alpha)+2*k[t]*f[t]*(1/90*alpha^2+1/30*alpha)+2*k[t]*u[t]*(1/15*alpha)
               +e[t]^2*(1/400+9/400*alpha^2+3/200*alpha)+2*e[t]*h[t]*(1/100+3/100*alpha^2+1/25*alpha)+2*e[t]*f[t]*(1/40+1/40*alpha^2+1/12*alpha)+2*e[t]*u[t]*(1/20+3/20*alpha)
               +h[t]^2*(1/25+1/25*alpha^2+2/25*alpha)+2*h[t]*f[t]*(1/10+1/30*alpha^2+2/15*alpha)+2*h[t]*u[t]*(1/5+1/5*alpha)
               +f[t]^2*(1/4+1/36*alpha^2+1/6*alpha)+2*f[t]*u[t]*(1/2+1/6*alpha)
               +u[t]^2*1)
  }
  freq <- rbind(m,n,k,e,h,f,u)
  return(freq[, niter])
}

#' Chi-squared test for HWE using Wang's approach.
#'
#' @inheritParams hex_recursive
#' @param nind The number of individuals in the sample.
#' @param method Should we calculate the chi-squared statistic the correct way
#'     (\code{"correct"}) or the way that Wang et al calculated it
#'     (\code{"incorrect"})?
#'
#' @author Jing Wang and David Gerard
#'
#' @references
#' \itemize{
#'   \item{J. Wang, L. Feng, S. Mu, A. Dong, J. Gan, Z. Wen, J. Meng, M. Li, R. Wu, and L. Sun. Asymptotic tests for Hardy-Weinberg equilibrium in hexaploids. \emph{Horticulture Research}, 9, 05 2022. \doi{10.1093/hr/uhac104}.}
#' }
#'
#' @examples
#' nvec <- c(29, 21, 17, 10, 10, 10, 23)
#' nind <- sum(nvec)
#' yww <- nvec / nind
#' hex_chisq(yww = yww, nind = nind, method = "correct")
#'
#' ## reproduce 6.602 value from page 5 of Wang et al (2022)
#' hex_chisq(yww = yww, nind = 120, method = "incorrect")
#'
#' @export
hex_chisq <- function(yww, nind, niter = 8, alpha = 0, method = c("correct", "incorrect")) {
  method <- match.arg(method)
  rvec <- hex_recursive(yww = yww, niter = niter, alpha = alpha)
  chisq <- sum((yww - rvec)^2 / rvec)
  if (method == "correct") {
    chisq <- nind * chisq
  }
  df <- 6
  p <- pchisq(q = chisq, df = df, lower.tail = FALSE)
  return(list(chisq = chisq, df = df, p = p))
}

#' EM algorithm to estimate gamete frequencies
#'
#' This code was really really weird. I had to edit it slightly to make it
#' generalizable beyond the specific data of (29, 21, 17, 10, 10, 10, 23).
#'
#' @inheritParams hex_recursive
#'
#' @return The estimated gamete frequencies
#'
#' @author Jing Wang and David Gerard
#'
#' @return A list of length 2.
#' \itemize{
#'   \item{\code{p}: The estimated gamete frequencies.}
#'   \item{\code{q}: The estimated genotype frequencies.}
#' }
#'
#' @examples
#' yww <- c(29, 21, 17, 10, 10, 10, 23) / 120
#' hex_em(yww = yww)
#'
#' @export
hex_em <- function(yww, niter = 30) {
  f3 <- c()
  f2 <- c()
  f1 <- c()
  f0 <- c()
  p3= 0.2; p2=0.35; p1=0.15; p0=0.3 # initial values, for some reason

  p222 <- yww[[1]]
  p221 <- yww[[2]]
  p220 <- yww[[3]]
  p210 <- yww[[4]]
  p200 <- yww[[5]]
  p100 <- yww[[6]]
  p000 <- yww[[7]]

  for(i in 1:niter){
    fai6= 2*p3*p1/(p2^2+2*p3*p1)
    fai5= 2*p2^2/(p2^2+2*p3*p1)

    fai4= 2*p3*p0/(2*p3*p0+2*p2*p1)
    fai3= 2*p2*p1/(2*p3*p0+2*p2*p1)

    fai2= 2*p2*p0/(2*p2*p0+p1^2)
    fai1= 2*p1^2/(2*p2*p0+p1^2)

    p3=1/2*(2*p222+p221+fai6*p220+fai4*p210)
    p2=1/2*(p221+fai5*p220+fai3*p210+fai2*p200)
    p1=1/2*(fai6*p220+fai3*p210+fai1*p200+p100)
    p0=1/2*(fai4*p210+fai2*p200+p100+2*p000)

    f3 <- c(f3,p3)
    f2 <- c(f2,p2)
    f1 <- c(f1,p1)
    f0 <- c(f0,p0)

    p3 <- p3
    p2 <- p2
    p1 <- p1
    p0 <- p0
  }

  p3 <- f3[niter];p2 <- f2[niter];p1 <- f1[niter];p0 <- f0[niter]

  p222= (p3)^2
  p221=2*p3*p2
  p220=2*p3*p1+(p2)^2
  p210=2*p3*p0+2*p2*p1
  p200=2*p2*p0+(p1)^2
  p100=2*p1*p0
  p000=(p0)^2

  retlist <- list(p = c(p3, p2, p1, p0),
                  q = c(p222, p221, p220, p210, p200, p100, p000))

  return(retlist)
}

#' Test for random mating in autohexaploids
#'
#' Calls \code{\link{hex_em}()} and runs a chi-squared test on the resulting
#' genotype frequencies.
#'
#' @inheritParams hex_chisq
#'
#' @author Jing Wang and David Gerard
#'
#' @references
#' \itemize{
#'   \item{J. Wang, L. Feng, S. Mu, A. Dong, J. Gan, Z. Wen, J. Meng, M. Li, R. Wu, and L. Sun. Asymptotic tests for Hardy-Weinberg equilibrium in hexaploids. \emph{Horticulture Research}, 9, 05 2022. \doi{10.1093/hr/uhac104}.}
#' }
#'
#' @examples
#' nvec <- c(29, 21, 17, 10, 10, 10, 23)
#' nind <- sum(nvec)
#' yww <- nvec / nind
#' hex_rmtest(yww = yww, nind = nind, method = "correct")
#'
#' ## reproduce 6.649 value from page 5 of Wang et al (2022)
#' hout <- hex_em(yww = yww, niter = 30)
#' rvec <- hex_recursive(yww = hout$q, niter = 8, alpha = 0)
#' sum((yww - rvec)^2 / rvec)
#'
#' @export
hex_rmtest <- function(yww, nind, niter = 30, method = c("correct", "incorrect")) {
  method <- match.arg(method)
  hout <- hex_em(yww = yww, niter = niter)
  rvec <- hout$q
  chisq <- sum((yww - rvec)^2 / rvec)
  if (method == "correct") {
    chisq <- nind * chisq
    df <- 3
  } else {
    df <- 5.5 ## they say "between 5 to 6"
  }
  p <- pchisq(q = chisq, df = df, lower.tail = FALSE)
  return(list(chisq = chisq, df = df, p = p))
}

#' One generation of random mating
#'
#' @inheritParams hex_recursive
#'
#' @author Jing Wang and David Gerard
#'
#' @return The next generation's genotype frequencies
#'
#' @author David Gerard
#'
#' @export
hex_onegen <- function(yww, alpha) {

  m <- yww[1]
  n <- yww[2]
  k <- yww[3]
  e <- yww[4]
  h <- yww[5]
  f <- yww[6]
  u <- yww[7]

  ynext <- rep(NA_real_, 7)
  ynext[[1]] <- (m^2*1+2*m*n*(1/2+1/6*alpha)+2*m*k*(1/5+1/5*alpha)+2*m*e*(1/20+3/20*alpha)+2*m*h*(1/15*alpha)+2*m*f*0+2*m*u*0
             +n^2*(1/4+1/36*alpha^2+1/6*alpha)+2*n*k*(1/10+1/30*alpha^2+2/15*alpha)+2*n*e*(1/40+1/40*alpha^2+1/12*alpha)+2*n*h*(1/90*alpha^2+1/30*alpha)+2*n*f*0+2*n*u*0
             +k^2*(1/25+1/25*alpha^2+2/25*alpha)+2*k*e*(1/100+3/100*alpha^2+1/25*alpha)+2*k*h*(1/75*alpha^2+1/75*alpha)+2*k*f*0+2*k*u*0
             +e^2*(1/400+9/400*alpha^2+3/200*alpha)+2*e*h*(1/100*alpha^2+1/300*alpha)+2*e*f*0+2*e*u*0
             +h^2*(1/225*alpha^2)+2*h*f*0+2*h*u*0
             +f^2*0+2*f*u*0
             +u^2*0)

  ynext[[2]] <- (m^2*0+2*m*n*(1/2-1/3*alpha)+2*m*k*(3/5-1/3*alpha)+2*m*e*(9/20-3/20*alpha)+2*m*h*(1/5+1/15*alpha)+2*m*f*(1/6*alpha)+2*m*u*0
             +n^2*(1/2-1/9*alpha^2-1/6*alpha)+2*n*k*(2/5-11/90*alpha^2-1/30*alpha)+2*n*e*(1/4-3/40*alpha^2+7/120*alpha)+2*n*h*(1/10-1/90*alpha^2+1/10*alpha)+2*n*f*(1/36*alpha^2+1/12*alpha)+2*n*u*0
             +k^2*(6/25-2/15*alpha^2+8/75*alpha)+2*k*e*(3/25-2/25*alpha^2+2/15*alpha)+2*k*h*(1/25-2/225*alpha^2+7/75*alpha)+2*k*f*(1/30*alpha^2+1/30*alpha)+2*k*u*0
             +e^2*(9/200-9/200*alpha^2+3/25*alpha)+2*e*h*(1/100+19/300*alpha)+2*e*f*(1/40*alpha^2+1/120*alpha)+2*e*u*0
             +h^2*(2/225*alpha^2+2/75*alpha)+2*h*f*(1/90*alpha^2)+2*h*u*0
             +f^2*0+2*f*u*0
             +u^2*0)

  ynext[[3]] <- (m^2*0+2*m*n*(1/6*alpha)+2*m*k*(1/5+1/15*alpha)+2*m*e*(9/20-3/20*alpha)+2*m*h*(3/5-1/3*alpha)+2*m*f*(1/2-1/3*alpha)+2*m*u*0
             +n^2*(1/4+1/6*alpha^2-1/6*alpha)+2*n*k*(2/5+7/45*alpha^2-4/15*alpha)+2*n*e*(9/20+1/20*alpha^2-13/60*alpha)+2*n*h*(2/5-1/15*alpha^2-1/10*alpha)+2*n*f*(1/4-1/9*alpha^2)+2*n*u*0
             +k^2*(11/25+31/225*alpha^2-22/75*alpha)+2*k*e*(37/100+3/100*alpha^2-11/75*alpha)+2*k*h*(6/25-19/225*alpha^2+1/25*alpha)+2*k*f*(1/10-11/90*alpha^2+2/15*alpha)+2*k*u*0
             +e^2*(99/400-9/400*alpha^2-3/200*alpha)+2*e*h*(3/25-7/100*alpha^2+31/300*alpha)+2*e*f*(1/40-3/40*alpha^2+2/15*alpha)+2*e*u*0
             +h^2*(1/25-9/225*alpha^2+8/75*alpha)+2*h*f*((-1/90)*alpha^2+1/15*alpha)+2*h*u*0
             +f^2*(1/36*alpha^2)+2*f*u*0
             +u^2*0)

  ynext[[4]] <- (m^2*0+2*m*n*0+2*m*k*(1/15*alpha)+2*m*e*(1/20+3/20*alpha)+2*m*h*(1/5+1/5*alpha)+2*m*f*(1/2+1/6*alpha)+2*m*u*1
             +n^2*((-1/9)*alpha^2+1/6*alpha)+2*n*k*(1/10-1/15*alpha^2+1/10*alpha)+2*n*e*(1/4+1/20*alpha^2-1/15*alpha)+2*n*h*(2/5+7/45*alpha^2-1/5*alpha)+2*n*f*(1/2+1/6*alpha^2-1/6*alpha)+2*n*u*(1/2+1/6*alpha)
             +k^2*(6/25-4/225*alpha^2-2/75*alpha)+2*k*e*(37/100+2/25*alpha^2-59/300*alpha)+2*k*h*(11/25+4/25*alpha^2-22/75*alpha)+2*k*f*(2/5+7/45*alpha^2-1/5*alpha)+2*k*u*(1/5+1/5*alpha)
             +e^2*(41/100+9/100*alpha^2-6/25*alpha)+2*e*h*(37/100+2/25*alpha^2-59/300*alpha)+2*e*f*(1/4+1/20*alpha^2-1/15*alpha)+2*e*u*(1/20+3/20*alpha)
             +h^2*(6/25-4/225*alpha^2-2/75*alpha)+2*h*f*(1/10-1/15*alpha^2+1/10*alpha)+2*h*u*(1/15*alpha)
             +f^2*((-1/9)*alpha^2+1/6*alpha)+2*f*u*0
             +u^2*0)

  ynext[[5]] <- (m^2*0+2*m*n*0+2*m*k*0+2*m*e*0+2*m*h*0+2*m*f*0+2*m*u*0
             +n^2*(1/36*alpha^2)+2*n*k*((-1/90)*alpha^2+1/15*alpha)+2*n*e*(1/40-3/40*alpha^2+2/15*alpha)+2*n*h*(1/10-11/90*alpha^2+2/15*alpha)+2*n*f*(1/4-1/9*alpha^2)+2*n*u*(1/2-1/3*alpha)
             +k^2*(1/25-9/225*alpha^2+8/75*alpha)+2*k*e*(3/25-7/100*alpha^2+31/300*alpha)+2*k*h*(6/25-19/225*alpha^2+1/25*alpha)+2*k*f*(2/5-1/15*alpha^2-1/10*alpha)+2*k*u*(3/5-1/3*alpha)
             +e^2*(99/400-9/400*alpha^2-3/200*alpha)+2*e*h*(37/100+3/100*alpha^2-11/75*alpha)+2*e*f*(9/20+1/20*alpha^2-13/60*alpha)+2*e*u*(9/20-3/20*alpha)
             +h^2*(11/25+31/225*alpha^2-22/75*alpha)+2*h*f*(2/5+7/45*alpha^2-4/15*alpha)+2*h*u*(1/5+1/15*alpha)
             +f^2*(1/4+1/6*alpha^2-1/6*alpha)+2*f*u*(1/6*alpha)
             +u^2*0)

  ynext[[6]] <- (m^2*0+2*m*n*0+2*m*k*0+2*m*e*0+2*m*h*0+2*m*f*0+2*m*u*0
             +n^2*0+2*n*k*(1/90*alpha^2)+2*n*e*(1/40*alpha^2+1/120*alpha)+2*n*h*(1/30*alpha^2+1/30*alpha)+2*n*f*(1/36*alpha^2+1/12*alpha)+2*n*u*(1/6*alpha)
             +k^2*(2/225*alpha^2+2/75*alpha)+2*k*e*(1/100+19/300*alpha)+2*k*h*(1/25-2/225*alpha^2+7/75*alpha)+2*k*f*(1/10-1/90*alpha^2+1/10*alpha)+2*k*u*(1/5+1/15*alpha)
             +e^2*(9/200-9/200*alpha^2+3/25*alpha)+2*e*h*(3/25-2/25*alpha^2+2/15*alpha)+2*e*f*(1/4-3/40*alpha^2+7/120*alpha)+2*e*u*(9/20-3/20*alpha)
             +h^2*(6/25-2/15*alpha^2+8/75*alpha)+2*h*f*(2/5-11/90*alpha^2-1/30*alpha)+2*h*u*(3/5-1/3*alpha)
             +f^2*(1/2-1/9*alpha^2-1/6*alpha)+2*f*u*(1/2-1/3*alpha)
             +u^2*0)

  ynext[[7]] <- (m^2*0+2*m*n*0+2*m*k*0+2*m*e*0+2*m*h*0+2*m*f*0+2*m*u*0
             +n^2*0+2*n*k*0+2*n*e*0+2*n*h*0+2*n*f*0+2*n*u*0
             +k^2*(1/225*alpha^2)+2*k*e*(1/100*alpha^2+1/300*alpha)+2*k*h*(1/75*alpha^2+1/75*alpha)+2*k*f*(1/90*alpha^2+1/30*alpha)+2*k*u*(1/15*alpha)
             +e^2*(1/400+9/400*alpha^2+3/200*alpha)+2*e*h*(1/100+3/100*alpha^2+1/25*alpha)+2*e*f*(1/40+1/40*alpha^2+1/12*alpha)+2*e*u*(1/20+3/20*alpha)
             +h^2*(1/25+1/25*alpha^2+2/25*alpha)+2*h*f*(1/10+1/30*alpha^2+2/15*alpha)+2*h*u*(1/5+1/5*alpha)
             +f^2*(1/4+1/36*alpha^2+1/6*alpha)+2*f*u*(1/2+1/6*alpha)
             +u^2*1)

  return(ynext)
}

#' Test for HWE
#'
#' Their code is super weird, so I try to reproduce it using simpler code.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{J. Wang, L. Feng, S. Mu, A. Dong, J. Gan, Z. Wen, J. Meng, M. Li, R. Wu, and L. Sun. Asymptotic tests for Hardy-Weinberg equilibrium in hexaploids. \emph{Horticulture Research}, 9, 05 2022. \doi{10.1093/hr/uhac104}.}
#' }
#'
#' @examples
#' nvec <- c(29, 21, 17, 10, 10, 10, 23)
#' nind <- sum(nvec)
#' yww <- nvec / nind
#'
#' ## this should be 5.922
#' hex_drtest(yww = yww, nind = nind, niter = 8)
#'
#' @noRd
hex_drtest <- function(yww, nind, niter = 8) {
  hout <- hex_em(yww = yww)
  q1 <- hout$q
}


#' Wang's code to estimate double reduction rate
#'
#' @param NN A vector of genotype counts.
#' @param niter The number of iterations in the EM algorithm.
#' @param tol The stopping criterion for the EM algorithm.
#'
#' @author Jing Wang and David Gerard
#'
#' @examples
#' ## Generate data under assumed mode
#' alpha <- 0.2
#' p <- c(0.3, 0.2, 0.35, 0.15)
#' p <- p / sum(p)
#' q <- stats::convolve(p, rev(p), type = "open")
#' q1 <- hex_onegen(q, alpha = alpha)
#' NN <- c(stats::rmultinom(n = 1, size = 100, prob = q1))
#' hex_estdr(NN = NN)
#' p
#'
#' @export
hex_estdr <- function(NN, niter = 100, tol = 10^-4) {
  n7 <- NN[1]; n6 <- NN[2]; n5 <- NN[3]; n4 <- NN[4]
  n3 <- NN[5]; n2 <- NN[6]; n1 <- NN[7]
  N=n7+n6+n5+n4+n3+n2+n1
  ## Initial values, for some reason
  p3 <- 0.3
  p2 <- 0.2
  p1 <- 0.35
  p0 <- 0.15
  alpha <- 0
  m <- p3^2
  n <- 2*p3*p2
  k <- p2^2+2*p3*p1
  e <- 2*p3*p0+2*p2*p1
  h <- 2*p2*p0+p1^2
  f <- 2*p1*p0
  u <- p0^2

  f4 <- c();f3 <- c();f2 <- c();f1 <- c();f0 <- c()
  for(i in 1:niter){
    alpha_old <- alpha

    Q7 <-expression(((p3^4+p3^2*p2^2+1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+1/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+2*p3^3*p2+2/5*(p3^2*p2^2+2*p3^3*p1)
                      +1/10*(2*p3^3*p0+2*p3^2*p2*p1)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)+1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/50*(2*p3*p2^2*p0
                                                                                                                         +4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2))+alpha*(1/3*2*p3^3*p2+2/5*(p3^2*p2^2+2*p3^3*p1)+3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                                                                                                                                       +2/15*(2*p3^2*p2*p0+p3^2*p1^2)+1/6*(4*p3^2*p2^2)+4/15*(2*p3*p2^3+4*p3^2*p2*p1)+1/6*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)
                                                                                                                                                                       +2/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+2/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                                                                                                                                                                       +3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))+alpha^2*(1/36*4*p3^2*p2^2
                                                                                                                                                                                                                                                                                         +1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/15*(2*p3*p2^3+4*p3^2*p2*p1)
                                                                                                                                                                                                                                                                                         +1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)+3/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                                                                                                                                                                                                                                                         +2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))))

    Q6 <-expression(((2*p3^3*p2+6/5*(p3^2*p2^2+2*p3^3*p1)+9/10*(2*p3^3*p0+2*p3^2*p2*p1)+2/5*(2*p3^2*p2*p0+p3^2*p1^2)+1/2*(4*p3^2*p2^2)+4/5*(2*p3*p2^3+4*p3^2*p2*p1)+1/2*(4*p3^2*p2*p0+4*p3*p2^2*p1)
                      +1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+6/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+6/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                      +9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))+alpha*(-2/3*(2*p3^3*p2)-2/3*(p3^2*p2^2+2*p3^3*p1)-3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                                                                                                     +2/15*(2*p3^2*p2*p0+p3^2*p1^2)+1/3*(2*p3^2*p1*p0)-1/6*(4*p3^2*p2^2)-1/15*(2*p3*p2^3+4*p3^2*p2*p1)+7/60*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/6*(4*p3*p2*p1*p0)
                                                                                                                                     +8/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+4/15*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+14/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                                                                                                                                     +3/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+19/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/60*(4*p3*p1*p0^2+4*p2*p1^2*p0)+2/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4))
                     +alpha^2*(-1/9*(4*p3^2*p2^2)-11/45*(2*p3*p2^3+4*p3^2*p2*p1)-3/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)-1/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)
                               +1/18*(4*p3*p2*p1*p0)-2/15*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-4/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                               -4/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                               +1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+2/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/45*(4*p2*p1*p0^2+2*p1^3*p0))))

    Q5 <-  expression((2/5*(p3^2*p2^2+2*p3^3*p1)+9/10*(2*p3^3*p0+2*p3^2*p2*p1)+6/5*(2*p3^2*p2*p0+p3^2*p1^2)+(2*p3^2*p1*p0)+1/4*(4*p3^2*p2^2)+4/5*(2*p3*p2^3+4*p3^2*p2*p1)
                       +9/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/2*(4*p3*p2*p1*p0)+11/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                       +37/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+12/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                       +99/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+6/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                       +1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4))+alpha*(1/3*(2*p3^3*p2)+2/15*(p3^2*p2^2+2*p3^3*p1)-3/10*(2*p3^3*p0+2*p3^2*p2*p1)-2/3*(2*p3^2*p2*p0+p3^2*p1^2)
                                                                     -2/3*(2*p3^2*p1*p0)-1/6*(4*p3^2*p2^2)-8/15*(2*p3*p2^3+4*p3^2*p2*p1)-13/30*(4*p3^2*p2*p0+4*p3*p2^2*p1)-1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)-22/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                                     -22/75*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                                     +31/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+4/15*(4*p3*p1*p0^2+4*p2*p1^2*p0)+8/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+2/15*(4*p2*p1*p0^2+2*p1^3*p0))
                      +alpha^2*(1/6*(4*p3^2*p2^2)+14/45*(2*p3*p2^3+4*p3^2*p2*p1)+1/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)-2/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/9*(4*p3*p2*p1*p0)+31/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                +3/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)-38/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-11/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                -7/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-3/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)-9/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-1/45*(4*p2*p1*p0^2+2*p1^3*p0)+1/36*(4*p1^2*p0^2)))


    Q4 <- expression((1/10*(2*p3^3*p0+2*p3^2*p2*p1)+2/5*(2*p3^2*p2*p0+p3^2*p1^2)+(2*p3^2*p1*p0)+2*(p3^2*p0^2)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)
                      +1/2*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+(4*p3*p2*p1*p0)+(2*p3*p2*p0^2)
                      +6/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+37/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                      +22/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)+41/100*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                      +37/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/2*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/10*(2*p3*p0^3+2*p2*p1*p0^2)+6/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                      +1/5*(4*p2*p1*p0^2+2*p1^3*p0))+alpha*(2/15*(p3^2*p2^2+2*p3^3*p1)+3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                            +2/5*(2*p3^2*p2*p0+p3^2*p1^2)+1/3*(2*p3^2*p1*p0)+1/6*(4*p3^2*p2^2)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)-2/15*(4*p3^2*p2*p0+4*p3*p2^2*p1)
                                                            -2/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)-1/3*(4*p3*p2*p1*p0)+1/3*(2*p3*p2*p0^2)-2/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-59/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                            -44/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-2/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)
                                                            -6/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)-59/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-2/15*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                                                            +3/10*(2*p3*p0^3+2*p2*p1*p0^2)-2/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/15*(2*p2*p0^3+p1^2*p0^2)
                                                            +1/6*(4*p1^2*p0^2))+alpha^2*(-1/9*(4*p3^2*p2^2)-2/15*(2*p3*p2^3+4*p3^2*p2*p1)+1/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)+14/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/3*(4*p3*p2*p1*p0)
                                                                                         -4/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+4/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+8/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                                                                                         +14/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)+9/100*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+4/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                                                                                         -4/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-2/15*(4*p2*p1*p0^2+2*p1^3*p0)-1/9*(4*p1^2*p0^2)))

    Q3 <- expression(alpha^2*(1/36*(4*p3^2*p2^2)-1/45*(2*p3*p2^3+4*p3^2*p2*p1)-3/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)-11/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/9*(4*p3*p2*p1*p0)-9/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                              -7/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)-38/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-2/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                              +3/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)+31/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+14/45*(4*p2*p1*p0^2+2*p1^3*p0)
                              +1/6*(4*p1^2*p0^2))+alpha*(2/15*(2*p3*p2^3+4*p3^2*p2*p1)+4/15*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/3*(2*p3*p2*p0^2)+8/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                         +31/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)-2/3*(p2^2*p0^2+2*p3*p1*p0^2)
                                                         -3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)-22/75*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-13/30*(4*p3*p1*p0^2+4*p2*p1^2*p0)-3/10*(2*p3*p0^3+2*p2*p1*p0^2)
                                                         -22/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-8/15*(4*p2*p1*p0^2+2*p1^3*p0)+2/15*(2*p2*p0^3+p1^2*p0^2)-1/6*(4*p1^2*p0^2)+1/3*(2*p1*p0^3))+(1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/2*(4*p3*p2*p1*p0)
                                                                                                                                                                                             +(2*p3*p2*p0^2)+1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+6/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+12/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                                                                                                                                                                                             +6/5*(p2^2*p0^2+2*p3*p1*p0^2)+99/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+37/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+9/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)+9/10*(2*p3*p0^3+2*p2*p1*p0^2)
                                                                                                                                                                                             +11/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+4/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/4*(4*p1^2*p0^2)))

    Q2 <- expression(alpha^2*(1/45*(2*p3*p2^3+4*p3^2*p2*p1)+1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/18*(4*p3*p2*p1*p0)
                              +2/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-4/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-1/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                              -4/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-3/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)-2/15*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-11/45*(4*p2*p1*p0^2+2*p1^3*p0)
                              -1/9*(4*p1^2*p0^2))+alpha*(1/60*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/6*(4*p3*p2*p1*p0)+1/3*(2*p3*p2*p0^2)+2/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                         +19/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+14/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/15*(p2^2*p0^2+2*p3*p1*p0^2)+3/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                         +4/15*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+7/60*(4*p3*p1*p0^2+4*p2*p1^2*p0)-3/10*(2*p3*p0^3+2*p2*p1*p0^2)+8/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                                                         -1/15*(4*p2*p1*p0^2+2*p1^3*p0)-2/3*(2*p2*p0^3+p1^2*p0^2)-1/6*(4*p1^2*p0^2)-2/3*(2*p1*p0^3))+(1/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                                                                                                                      +2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)+9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                                                                                                                      +6/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/2*(4*p3*p1*p0^2+4*p2*p1^2*p0)+9/10*(2*p3*p0^3+2*p2*p1*p0^2)+6/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                                                                                                                                                      +4/5*(4*p2*p1*p0^2+2*p1^3*p0)+6/5*(2*p2*p0^3+p1^2*p0^2)+1/2*(4*p1^2*p0^2)+(2*p1*p0^3)))


    Q1 <-expression((1/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/10*(2*p3*p0^3+2*p2*p1*p0^2)
                     +1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/4*(4*p1^2*p0^2)+(2*p1*p0^3)+p0^4)
                    +alpha*(1/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                            +2/15*(p2^2*p0^2+2*p3*p1*p0^2)+3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+2/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/6*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                            +3/10*(2*p3*p0^3+2*p2*p1*p0^2)+2/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+4/15*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/6*(4*p1^2*p0^2)+1/3*(2*p1*p0^3))
                    +alpha^2*(1/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+1/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                              +1/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)+9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+3/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)
                              +1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/15*(4*p2*p1*p0^2+2*p1^3*p0)+1/36*(4*p1^2*p0^2)))



    R7 <-expression(((p3^4+p3^2*p2^2+1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+1/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+2*p3^3*p2+2/5*(p3^2*p2^2+2*p3^3*p1)
                      +1/10*(2*p3^3*p0+2*p3^2*p2*p1)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)+1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/50*(2*p3*p2^2*p0
                                                                                                                         +4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2))+0*(1/3*2*p3^3*p2+2/5*(p3^2*p2^2+2*p3^3*p1)+3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                                                                                                                                   +2/15*(2*p3^2*p2*p0+p3^2*p1^2)+1/6*(4*p3^2*p2^2)+4/15*(2*p3*p2^3+4*p3^2*p2*p1)+1/6*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)
                                                                                                                                                                   +2/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+2/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                                                                                                                                                                   +3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))+0^2*(1/36*4*p3^2*p2^2
                                                                                                                                                                                                                                                                                 +1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/15*(2*p3*p2^3+4*p3^2*p2*p1)
                                                                                                                                                                                                                                                                                 +1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)+3/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                                                                                                                                                                                                                                                 +2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))))

    R6 <-expression(((2*p3^3*p2+6/5*(p3^2*p2^2+2*p3^3*p1)+9/10*(2*p3^3*p0+2*p3^2*p2*p1)+2/5*(2*p3^2*p2*p0+p3^2*p1^2)+1/2*(4*p3^2*p2^2)+4/5*(2*p3*p2^3+4*p3^2*p2*p1)+1/2*(4*p3^2*p2*p0+4*p3*p2^2*p1)
                      +1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+6/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+6/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                      +9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3))+0*(-2/3*(2*p3^3*p2)-2/3*(p3^2*p2^2+2*p3^3*p1)-3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                                                                                                 +2/15*(2*p3^2*p2*p0+p3^2*p1^2)+1/3*(2*p3^2*p1*p0)-1/6*(4*p3^2*p2^2)-1/15*(2*p3*p2^3+4*p3^2*p2*p1)+7/60*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/6*(4*p3*p2*p1*p0)
                                                                                                                                 +8/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+4/15*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+14/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                                                                                                                                 +3/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+19/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/60*(4*p3*p1*p0^2+4*p2*p1^2*p0)+2/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4))
                     +0^2*(-1/9*(4*p3^2*p2^2)-11/45*(2*p3*p2^3+4*p3^2*p2*p1)-3/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)-1/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)
                           +1/18*(4*p3*p2*p1*p0)-2/15*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-4/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                           -4/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                           +1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+2/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/45*(4*p2*p1*p0^2+2*p1^3*p0))))

    R5 <-  expression((2/5*(p3^2*p2^2+2*p3^3*p1)+9/10*(2*p3^3*p0+2*p3^2*p2*p1)+6/5*(2*p3^2*p2*p0+p3^2*p1^2)+(2*p3^2*p1*p0)+1/4*(4*p3^2*p2^2)+4/5*(2*p3*p2^3+4*p3^2*p2*p1)
                       +9/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/2*(4*p3*p2*p1*p0)+11/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                       +37/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+12/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                       +99/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+6/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                       +1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4))+0*(1/3*(2*p3^3*p2)+2/15*(p3^2*p2^2+2*p3^3*p1)-3/10*(2*p3^3*p0+2*p3^2*p2*p1)-2/3*(2*p3^2*p2*p0+p3^2*p1^2)
                                                                 -2/3*(2*p3^2*p1*p0)-1/6*(4*p3^2*p2^2)-8/15*(2*p3*p2^3+4*p3^2*p2*p1)-13/30*(4*p3^2*p2*p0+4*p3*p2^2*p1)-1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)-22/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                                 -22/75*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                                 +31/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+4/15*(4*p3*p1*p0^2+4*p2*p1^2*p0)+8/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+2/15*(4*p2*p1*p0^2+2*p1^3*p0))
                      +0^2*(1/6*(4*p3^2*p2^2)+14/45*(2*p3*p2^3+4*p3^2*p2*p1)+1/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)-2/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/9*(4*p3*p2*p1*p0)+31/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                            +3/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)-38/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-11/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                            -7/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-3/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)-9/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-1/45*(4*p2*p1*p0^2+2*p1^3*p0)+1/36*(4*p1^2*p0^2)))


    R4 <- expression((1/10*(2*p3^3*p0+2*p3^2*p2*p1)+2/5*(2*p3^2*p2*p0+p3^2*p1^2)+(2*p3^2*p1*p0)+2*(p3^2*p0^2)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)
                      +1/2*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+(4*p3*p2*p1*p0)+(2*p3*p2*p0^2)
                      +6/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+37/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                      +22/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)+41/100*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                      +37/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/2*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/10*(2*p3*p0^3+2*p2*p1*p0^2)+6/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                      +1/5*(4*p2*p1*p0^2+2*p1^3*p0))+0*(2/15*(p3^2*p2^2+2*p3^3*p1)+3/10*(2*p3^3*p0+2*p3^2*p2*p1)
                                                        +2/5*(2*p3^2*p2*p0+p3^2*p1^2)+1/3*(2*p3^2*p1*p0)+1/6*(4*p3^2*p2^2)+1/5*(2*p3*p2^3+4*p3^2*p2*p1)-2/15*(4*p3^2*p2*p0+4*p3*p2^2*p1)
                                                        -2/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)-1/3*(4*p3*p2*p1*p0)+1/3*(2*p3*p2*p0^2)-2/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-59/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                        -44/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-2/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)
                                                        -6/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)-59/150*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-2/15*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                                                        +3/10*(2*p3*p0^3+2*p2*p1*p0^2)-2/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/15*(2*p2*p0^3+p1^2*p0^2)
                                                        +1/6*(4*p1^2*p0^2))+0^2*(-1/9*(4*p3^2*p2^2)-2/15*(2*p3*p2^3+4*p3^2*p2*p1)+1/10*(4*p3^2*p2*p0+4*p3*p2^2*p1)+14/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/3*(4*p3*p2*p1*p0)
                                                                                 -4/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+4/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+8/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                                                                                 +14/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)+9/100*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+4/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                                                                                 -4/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-2/15*(4*p2*p1*p0^2+2*p1^3*p0)-1/9*(4*p1^2*p0^2)))

    R3 <- expression(0^2*(1/36*(4*p3^2*p2^2)-1/45*(2*p3*p2^3+4*p3^2*p2*p1)-3/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)-11/45*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/9*(4*p3*p2*p1*p0)-9/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                          -7/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)-38/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-2/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                          +3/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)+31/225*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+14/45*(4*p2*p1*p0^2+2*p1^3*p0)
                          +1/6*(4*p1^2*p0^2))+0*(2/15*(2*p3*p2^3+4*p3^2*p2*p1)+4/15*(4*p3^2*p2*p0+4*p3*p2^2*p1)+4/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)-2/3*(2*p3*p2*p0^2)+8/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                 +31/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)-2/3*(p2^2*p0^2+2*p3*p1*p0^2)
                                                 -3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)-22/75*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-13/30*(4*p3*p1*p0^2+4*p2*p1^2*p0)-3/10*(2*p3*p0^3+2*p2*p1*p0^2)
                                                 -22/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-8/15*(4*p2*p1*p0^2+2*p1^3*p0)+2/15*(2*p2*p0^3+p1^2*p0^2)-1/6*(4*p1^2*p0^2)+1/3*(2*p1*p0^3))+(1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/5*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/2*(4*p3*p2*p1*p0)
                                                                                                                                                                                     +(2*p3*p2*p0^2)+1/25*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+6/25*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+12/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+4/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                                                                                                                                                                                     +6/5*(p2^2*p0^2+2*p3*p1*p0^2)+99/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+37/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+9/10*(4*p3*p1*p0^2+4*p2*p1^2*p0)+9/10*(2*p3*p0^3+2*p2*p1*p0^2)
                                                                                                                                                                                     +11/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+4/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/4*(4*p1^2*p0^2)))

    R2 <- expression(0^2*(1/45*(2*p3*p2^3+4*p3^2*p2*p1)+1/20*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/18*(4*p3*p2*p1*p0)
                          +2/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)-4/225*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)-1/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)-9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                          -4/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)-3/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)-2/15*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)-11/45*(4*p2*p1*p0^2+2*p1^3*p0)
                          -1/9*(4*p1^2*p0^2))+0*(1/60*(4*p3^2*p2*p0+4*p3*p2^2*p1)+1/15*(4*p3*p2^2*p0+2*p3*p2*p1^2)+1/6*(4*p3*p2*p1*p0)+1/3*(2*p3*p2*p0^2)+2/75*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)
                                                 +19/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+14/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/15*(p2^2*p0^2+2*p3*p1*p0^2)+3/25*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                 +4/15*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+7/60*(4*p3*p1*p0^2+4*p2*p1^2*p0)-3/10*(2*p3*p0^3+2*p2*p1*p0^2)+8/75*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                                                 -1/15*(4*p2*p1*p0^2+2*p1^3*p0)-2/3*(2*p2*p0^3+p1^2*p0^2)-1/6*(4*p1^2*p0^2)-2/3*(2*p1*p0^3))+(1/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)
                                                                                                                                              +2/25*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/5*(2*p2^2*p1*p0+4*p3*p1^2*p0)+2/5*(p2^2*p0^2+2*p3*p1*p0^2)+9/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)
                                                                                                                                              +6/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/2*(4*p3*p1*p0^2+4*p2*p1^2*p0)+9/10*(2*p3*p0^3+2*p2*p1*p0^2)+6/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)
                                                                                                                                              +4/5*(4*p2*p1*p0^2+2*p1^3*p0)+6/5*(2*p2*p0^3+p1^2*p0^2)+1/2*(4*p1^2*p0^2)+(2*p1*p0^3)))


    R1 <-expression((1/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+1/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/10*(2*p3*p0^3+2*p2*p1*p0^2)
                     +1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/5*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/4*(4*p1^2*p0^2)+(2*p1*p0^3)+p0^4)
                    +0*(1/150*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)+1/15*(2*p2^2*p1*p0+4*p3*p1^2*p0)
                        +2/15*(p2^2*p0^2+2*p3*p1*p0^2)+3/200*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+2/25*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)+1/6*(4*p3*p1*p0^2+4*p2*p1^2*p0)
                        +3/10*(2*p3*p0^3+2*p2*p1*p0^2)+2/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+4/15*(4*p2*p1*p0^2+2*p1^3*p0)+2/5*(2*p2*p0^3+p1^2*p0^2)+1/6*(4*p1^2*p0^2)+1/3*(2*p1*p0^3))
                    +0^2*(1/225*(p2^4+4*p3*p2^2*p1+4*p3^2*p1^2)+1/50*(2*p3*p2^2*p0+4*p3^2*p1*p0+2*p2^3*p1+4*p3*p2*p1^2)+2/75*(2*p2^3*p0+4*p3*p2*p1*p0+p1^2*p2^2+2*p3*p1^3)
                          +1/45*(2*p2^2*p1*p0+4*p3*p1^2*p0)+9/400*(4*p3^2*p0^2+8*p3*p2*p1*p0+4*p2^2*p1^2)+3/50*(4*p3*p2*p0^2+4*p2^2*p1*p0+2*p3*p1^2*p0+2*p2*p1^3)
                          +1/20*(4*p3*p1*p0^2+4*p2*p1^2*p0)+1/25*(4*p2^2*p0^2+4*p2*p1^2*p0+p1^4)+1/15*(4*p2*p1*p0^2+2*p1^3*p0)+1/36*(4*p1^2*p0^2)))


    Q.fai6A <- expression(

      ((p3^2)^2*(1)+(2*(p3^2)*(2*p3*p2))*(2/3*alpha)+(2*(p3^2)*(p2^2+2*p3*p1))*(2/5*alpha)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(1/5*alpha)+(2*(p3^2)*(2*p2*p0+p1^2))*(1/15*alpha)+(2*(p3^2)*(2*p1*p0))*(0)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*((2*alpha*(1-alpha))/3+(4*(2*alpha^2))/9)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((alpha*(1-alpha))/3+(4*(2*alpha^2))/15)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*((2*alpha*(1-alpha))/15+(2*(2*alpha^2))/15)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/30+(2*(2*alpha^2))/45)+(2*(2*p3*p2)*(2*p1*p0))*(0)+(2*(2*p3*p2)*(p0^2))*(0)
       +((p2^2+2*p3*p1)^2)*((4*alpha*(1-alpha))/25+(4*(2*alpha^2))/25)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((3*alpha*(1-alpha))/50+(2*(2*alpha^2))/25)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/75+(2*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*(0)+(2*(p2^2+2*p3*p1)*(p0^2))*(0)
       +((2*p3*p0+2*p2*p1)^2)*((alpha*(1-alpha))/50+(2*alpha^2)/25)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/300+(2*alpha^2)/75)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*(0)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(0)
       +((2*p2*p0+p1^2)^2)*((2*alpha^2)/225)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*(0)+(2*(2*p2*p0+p1^2)*(p0^2))*(0)
       +((2*p1*p0)^2)*(0)+(2*(2*p1*p0)*(p0^2))*(0)
       +((p0^2)^2)*(0))

    )

    Q.fai5A1a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(1/6*alpha)+(2*(p3^2)*(p2^2+2*p3*p1))*(4/15*alpha)+(2*(p3^2)*e)*(3/10*alpha)+(2*(p3^2)*(2*p2*p0+p1^2))*(4/15*alpha)+(2*(p3^2)*(2*p1*p0))*(1/6*alpha)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*((5*alpha*(1-alpha))/6+(2*(2*alpha^2))/9)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((23*alpha*(1-alpha))/30+(11*(2*alpha^2))/45)+(2*(2*p3*p2)*e)*((67*alpha*(1-alpha))/120+(7*(2*alpha^2))/30)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((3*alpha*(1-alpha))/10+(17*(2*alpha^2))/90)+(2*(2*p3*p2)*(2*p1*p0))*((alpha*(1-alpha))/12+(2*alpha^2)/9)+(2*(2*p3*p2)*(p0^2))*(0)
       +((p2^2+2*p3*p1)^2)*((44*alpha*(1-alpha))/75+(16*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*e)*((28*alpha*(1-alpha))/75+(13*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((13*alpha*(1-alpha))/75+(28*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((alpha*(1-alpha))/30+(2*alpha^2)/15)+(2*(p2^2+2*p3*p1)*(p0^2))*(0)
       +(e^2)*((21*alpha*(1-alpha))/100+(3*(2*alpha^2))/25)+(2*e*(2*p2*p0+p1^2))*((alpha*(1-alpha))/12+(11*(2*alpha^2))/150)+(2*e*(2*p1*p0))*((alpha*(1-alpha))/120+(2*alpha^2)/30)+(2*e*(p0^2))*(0)
       +((2*p2*p0+p1^2)^2)*((2*alpha*(1-alpha))/75+(8*(2*alpha^2))/225)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((2*alpha^2)/90)+(2*(2*p2*p0+p1^2)*(p0^2))*(0)
       +((2*p1*p0)^2)*(0)+(2*(2*p1*p0)*(p0^2))*(0)
       +((p0^2)^2)*(0))

    )

    Q.fai4A2a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(1/6*alpha)+(2*(p3^2)*(p2^2+2*p3*p1))*(4/15*alpha)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(3/10*alpha)+(2*(p3^2)*(2*p2*p0+p1^2))*(4/15*alpha)+(2*(p3^2)*(2*p1*p0))*(1/6*alpha)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*((alpha*(1-alpha))/3+(2*alpha^2)/4)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((8*alpha*(1-alpha))/15+(13*(2*alpha^2))/45)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*((41*alpha*(1-alpha))/60+(17*(2*alpha^2))/60)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((7*alpha*(1-alpha))/10+(7*(2*alpha^2))/30)+(2*(2*p3*p2)*(2*p1*p0))*((alpha*(1-alpha))/2+(5*(2*alpha^2))/36)+(2*(2*p3*p2)*(p0^2))*(0)
       +((p2^2+2*p3*p1)^2)*((44*alpha*(1-alpha))/75+(64*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((89*alpha*(1-alpha))/150+(19*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((13*alpha*(1-alpha))/25+(44*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((alpha*(1-alpha))/3+(2*alpha^2)/9)+(2*(p2^2+2*p3*p1)*(p0^2))*(0)
       +((2*p3*p0+2*p2*p1)^2)*((12*alpha*(1-alpha))/25+(21*(2*alpha^2))/100)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((103*alpha*(1-alpha))/300+(23*(2*alpha^2))/150)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*((11*alpha*(1-alpha))/60+(2*alpha^2)/12)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(0)
       +((2*p2*p0+p1^2)^2)*((14*alpha*(1-alpha))/75+(8*(2*alpha^2))/75)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((alpha*(1-alpha))/15+(2*alpha^2)/18)+(2*(2*p2*p0+p1^2)*(p0^2))*(0)
       +((2*p1*p0)^2)*((2*alpha^2)/36)+(2*(2*p1*p0)*(p0^2))*(0)
       +((p0^2)^2)*(0))
    )

    Q.fai3A3a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(0)+(2*(p3^2)*(p2^2+2*p3*p1))*(1/15*alpha)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(1/5*alpha)+(2*(p3^2)*(2*p2*p0+p1^2))*(2/5*alpha)+(2*(p3^2)*(2*p1*p0))*(2/3*alpha)+(2*(p3^2)*(p0^2))*(1)
       +((2*p3*p2)^2)*((alpha*(1-alpha))/6+(2*alpha^2)/18)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((3*alpha*(1-alpha))/10+(2*(2*alpha^2))/15)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*((13*alpha*(1-alpha))/30+(7*(2*alpha^2))/30)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((3*alpha*(1-alpha))/5+(16*(2*alpha^2))/45)+(2*(2*p3*p2)*(2*p1*p0))*((5*alpha*(1-alpha))/6+(2*alpha^2)/2)+(2*(2*p3*p2)*(p0^2))*(2/3*alpha)
       +((p2^2+2*p3*p1)^2)*((34*alpha*(1-alpha))/75+(44*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((163*alpha*(1-alpha))/300+(19*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((44*alpha*(1-alpha))/75+(23*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((3*alpha*(1-alpha))/5+(16*(2*alpha^2))/45)+(2*(p2^2+2*p3*p1)*(p0^2))*(2/5*alpha)
       +((2*p3*p0+2*p2*p1)^2)*((29*alpha*(1-alpha))/50+(13*(2*alpha^2))/50)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((163*alpha*(1-alpha))/300+(19*(2*alpha^2))/75)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*((13*alpha*(1-alpha))/30+(7*(2*alpha^2))/30)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(1/5*alpha)
       +((2*p2*p0+p1^2)^2)*((34*alpha*(1-alpha))/75+(44*(2*alpha^2))/225)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((3*alpha*(1-alpha))/10+(2*(2*alpha^2))/15)+(2*(2*p2*p0+p1^2)*(p0^2))*(1/15*alpha)
       +((2*p1*p0)^2)*((alpha*(1-alpha))/6+(2*alpha^2)/18)+(2*(2*p1*p0)*(p0^2))*(0)
       +((p0^2)^2)*(0))
    )

    Q.fai2A4a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(0)+(2*(p3^2)*(p2^2+2*p3*p1))*(0)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(0)+(2*(p3^2)*(2*p2*p0+p1^2))*(0)+(2*(p3^2)*(2*p1*p0))*(0)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*((2*alpha^2)/36)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((alpha*(1-alpha))/15+(2*alpha^2)/18)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*((11*alpha*(1-alpha))/60+(2*alpha^2)/12)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/3+(2*alpha^2)/9)+(2*(2*p3*p2)*(2*p1*p0))*((alpha*(1-alpha))/2+(5*(2*alpha^2))/36)+(2*(2*p3*p2)*(p0^2))*(1/6*alpha)
       +((p2^2+2*p3*p1)^2)*((14*alpha*(1-alpha))/75+(8*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((103*alpha*(1-alpha))/300+(23*(2*alpha^2))/150)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((13*alpha*(1-alpha))/25+(44*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((7*alpha*(1-alpha))/10+(7*(2*alpha^2))/30)+(2*(p2^2+2*p3*p1)*(p0^2))*(4/15*alpha)
       +((2*p3*p0+2*p2*p1)^2)*((12*alpha*(1-alpha))/25+(21*(2*alpha^2))/100)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((89*alpha*(1-alpha))/150+(19*(2*alpha^2))/75)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*((41*alpha*(1-alpha))/60+(17*(2*alpha^2))/60)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(3/10*alpha)
       +((2*p2*p0+p1^2)^2)*((44*alpha*(1-alpha))/75+(64*(2*alpha^2))/225)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((8*alpha*(1-alpha))/15+(13*(2*alpha^2))/45)+(2*(2*p2*p0+p1^2)*(p0^2))*(4/15*alpha)
       +((2*p1*p0)^2)*((alpha*(1-alpha))/3+(2*alpha^2)/4)+(2*(2*p1*p0)*(p0^2))*(1/6*alpha)
       +((p0^2)^2)*(0))
    )

    Q.fai1A5a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(0)+(2*(p3^2)*(p2^2+2*p3*p1))*(0)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(0)+(2*(p3^2)*(2*p2*p0+p1^2))*(0)+(2*(p3^2)*(2*p1*p0))*(0)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*(0)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*((2*alpha^2)/90)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*((alpha*(1-alpha))/120+(2*alpha^2)/30)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/30+(2*alpha^2)/15)+(2*(2*p3*p2)*(2*p1*p0))*((alpha*(1-alpha))/12+(2*alpha^2)/9)+(2*(2*p3*p2)*(p0^2))*(1/6*alpha)
       +((p2^2+2*p3*p1)^2)*((2*alpha*(1-alpha))/75+(8*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((alpha*(1-alpha))/12+(11*(2*alpha^2))/150)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((13*alpha*(1-alpha))/75+(28*(2*alpha^2))/225)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((3*alpha*(1-alpha))/10+(17*(2*alpha^2))/90)+(2*(p2^2+2*p3*p1)*(p0^2))*(4/15*alpha)
       +((2*p3*p0+2*p2*p1)^2)*((21*alpha*(1-alpha))/100+(3*(2*alpha^2))/25)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((28*alpha*(1-alpha))/75+(13*(2*alpha^2))/75)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*((67*alpha*(1-alpha))/120+(7*(2*alpha^2))/30)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(3/10*alpha)
       +((2*p2*p0+p1^2)^2)*((44*alpha*(1-alpha))/75+(16*(2*alpha^2))/75)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((23*alpha*(1-alpha))/30+(11*(2*alpha^2))/45)+(2*(2*p2*p0+p1^2)*(p0^2))*(4/15*alpha)
       +((2*p1*p0)^2)*((5*alpha*(1-alpha))/6+(2*(2*alpha^2))/9)+(2*(2*p1*p0)*(p0^2))*(1/6*alpha)
       +((p0^2)^2)*(0))
    )

    Q.fai6a <- expression(

      ((p3^2)^2*(0)+(2*(p3^2)*(2*p3*p2))*(0)+(2*(p3^2)*(p2^2+2*p3*p1))*(0)+(2*(p3^2)*(2*p3*p0+2*p2*p1))*(0)+(2*(p3^2)*(2*p2*p0+p1^2))*(0)+(2*(p3^2)*(2*p1*p0))*(0)+(2*(p3^2)*(p0^2))*(0)
       +((2*p3*p2)^2)*(0)+(2*(2*p3*p2)*(p2^2+2*p3*p1))*(0)+(2*(2*p3*p2)*(2*p3*p0+2*p2*p1))*(0)+(2*(2*p3*p2)*(2*p2*p0+p1^2))*(0)+(2*(2*p3*p2)*(2*p1*p0))*(0)+(2*(2*p3*p2)*(p0^2))*(0)
       +((p2^2+2*p3*p1)^2)*((2*alpha^2)/225)+(2*(p2^2+2*p3*p1)*(2*p3*p0+2*p2*p1))*((alpha*(1-alpha))/300+(2*alpha^2)/75)+(2*(p2^2+2*p3*p1)*(2*p2*p0+p1^2))*((alpha*(1-alpha))/75+(2*(2*alpha^2))/75)+(2*(p2^2+2*p3*p1)*(2*p1*p0))*((alpha*(1-alpha))/30+(2*(2*alpha^2))/45)+(2*(p2^2+2*p3*p1)*(p0^2))*(1/15*alpha)
       +((2*p3*p0+2*p2*p1)^2)*((alpha*(1-alpha))/50+(2*alpha^2)/25)+(2*(2*p3*p0+2*p2*p1)*(2*p2*p0+p1^2))*((3*alpha*(1-alpha))/50+(2*(2*alpha^2))/25)+(2*(2*p3*p0+2*p2*p1)*(2*p1*p0))*((2*alpha*(1-alpha))/15+(2*(2*alpha^2))/15)+(2*(2*p3*p0+2*p2*p1)*(p0^2))*(1/5*alpha)
       +((2*p2*p0+p1^2)^2)*((4*alpha*(1-alpha))/25+(4*(2*alpha^2))/25)+(2*(2*p2*p0+p1^2)*(2*p1*p0))*((alpha*(1-alpha))/3+(4*(2*alpha^2))/15)+(2*(2*p2*p0+p1^2)*(p0^2))*(2/5*alpha)
       +((2*p1*p0)^2)*((2*alpha*(1-alpha))/3+(4*(2*alpha^2))/9)+(2*(2*p1*p0)*(p0^2))*(2/3*alpha)
       +((p0^2)^2)*(1))
    )

    fai6A <- eval(Q.fai6A)/eval(Q7)
    fai5A1a <- eval(Q.fai5A1a)/eval(Q6)
    fai4A2a <- eval(Q.fai4A2a)/eval(Q5)
    fai3A3a <- eval(Q.fai4A2a)/eval(Q4)
    fai2A4a <- eval(Q.fai2A4a)/eval(Q3)
    fai1A5a <- eval(Q.fai1A5a)/eval(Q2)
    fai6a <- eval(Q.fai6a)/eval(Q1)

    fai7.3 <-p3*eval(D(Q7,"p3"))/eval(R7)
    fai7.2 <-p2*eval(D(Q7,"p2"))/eval(R7)
    fai7.1 <-p1*eval(D(Q7,"p1"))/eval(R7)
    fai7.0 <-p0*eval(D(Q7,"p0"))/eval(R7)

    fai6.3 <-p3*eval(D(Q6,"p3"))/eval(R6)
    fai6.2 <-p2*eval(D(Q6,"p2"))/eval(R6)
    fai6.1 <-p1*eval(D(Q6,"p1"))/eval(R6)
    fai6.0 <-p0*eval(D(Q6,"p0"))/eval(R6)

    fai5.3 <-p3*eval(D(Q5,"p3"))/eval(R5)
    fai5.2 <-p2*eval(D(Q5,"p2"))/eval(R5)
    fai5.1 <-p1*eval(D(Q5,"p1"))/eval(R5)
    fai5.0 <-p0*eval(D(Q5,"p0"))/eval(R5)

    fai4.3 <-p3*eval(D(Q4,"p3"))/eval(R4)
    fai4.2 <-p2*eval(D(Q4,"p2"))/eval(R4)
    fai4.1 <-p1*eval(D(Q4,"p1"))/eval(R4)
    fai4.0 <-p0*eval(D(Q4,"p0"))/eval(R4)

    fai3.3 <-p3*eval(D(Q3,"p3"))/eval(R3)
    fai3.2 <-p2*eval(D(Q3,"p2"))/eval(R3)
    fai3.1 <-p1*eval(D(Q3,"p1"))/eval(R3)
    fai3.0 <-p0*eval(D(Q3,"p0"))/eval(R3)

    fai2.3 <-p3*eval(D(Q2,"p3"))/eval(R2)
    fai2.2 <-p2*eval(D(Q2,"p2"))/eval(R2)
    fai2.1 <-p1*eval(D(Q2,"p1"))/eval(R2)
    fai2.0 <-p0*eval(D(Q2,"p0"))/eval(R2)

    fai1.3 <-p3*eval(D(Q1,"p3"))/eval(R1)
    fai1.2 <-p2*eval(D(Q1,"p2"))/eval(R1)
    fai1.1 <-p1*eval(D(Q1,"p1"))/eval(R1)
    fai1.0 <-p0*eval(D(Q1,"p0"))/eval(R1)

    alpha = (fai6A*n7+fai5A1a*n6+fai4A2a*n5+fai3A3a*n4+fai2A4a*n3+fai1A5a*n2+fai6a*n1)/(2*N)
    PAAA = (fai7.3*n7+fai6.3*n6+fai5.3*n5+fai4.3*n4+fai3.3*n3+fai2.3*n2+fai1.3*n1)/(4*N)
    PAAa = (fai7.2*n7+fai6.2*n6+fai5.2*n5+fai4.2*n4+fai3.2*n3+fai2.2*n2+fai1.2*n1)/(4*N)
    PAaa = (fai7.1*n7+fai6.1*n6+fai5.1*n5+fai4.1*n4+fai3.1*n3+fai2.1*n2+fai1.1*n1)/(4*N)
    Paaa = (fai7.0*n7+fai6.0*n6+fai5.0*n5+fai4.0*n4+fai3.0*n3+fai2.0*n2+fai1.0*n1)/(4*N)

    f4 <-c(f4,alpha);f3 <-c(f3,PAAA);f2 <-c(f2,PAAa);f1 <-c(f1,PAaa);f0 <-c(f0,Paaa)
    p4 <- alpha ;p3 <- PAAA ;p2 <-  PAAa;p1 <- PAaa;p0 <- Paaa

    # cat(alpha, "\n")

    if (abs(alpha_old - alpha) < tol) {
      break
    }

  }

  retlist <- list(alpha = p4,
                  p = c(p3, p2, p1, p0))

  return(retlist)
}
