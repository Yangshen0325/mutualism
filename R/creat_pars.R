# creat pars
# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,M0,qgain,qloss,lambda1)
mutualism_pars <- create_mutualism_pars(lac_par = c(0.1,0.1),
                                      mu_par = c(0.2,0.2,0.5,0.5),
                                      K_par = c(Inf,Inf,0.5,0.5),
                                      gam_par = c(0.3,0.3),
                                      laa_par = c(0.4,0.4,0.5,0.5),
M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},
                                      qgain = 0.5,
                                      qloss = 0.5,
                                      lambda1 = 0.5)


create_mutualism_pars <- function(lac_par,
                              mu_par,
                              K_par,
                              gam_par,
                              laa_par,
                              M0,
                              qgain,
                              qloss,
                              lambda1) {
  testit::assert(is.numeric(lac_par))
  testit::assert(is.numeric(mu_par))
  testit::assert(is.numeric(K_par))
  testit::assert(is.numeric(gam_par))
  testit::assert(is.numeric(laa_par))
  testit::assert(is.matrix(M0))
  testit::assert(is.numeric(qgain))
  testit::assert(is.numeric(qloss))
  testit::assert(is.numeric(lambda1))
  testit::assert(lac_par >= 0.0)
  testit::assert(mu_par >= 0.0)
  testit::assert(K_par >= 0.0)
  testit::assert(gam_par >= 0.0)
  testit::assert(laa_par >= 0.0)
  testit::assert(qgain >=0.0)
  testit::assert(qloss >=0)
  testit::assert(lambda1 >=0)
  list(lac_par = lac_par,
       mu_par = mu_par,
       K_par = K_par,
       gam_par = gam_par,
       laa_par = laa_par,
       M0 = M0,
       qgain = qgain,
       qloss = qloss,
       lambda1 = lambda1)
}