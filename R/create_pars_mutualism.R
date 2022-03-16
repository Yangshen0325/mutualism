# creat pars
# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,qgain,qloss,lambda1, M0, pro)
# lac_par <- c(lac_plant, lac_animal)
# mu_par <- c(mu_P0, mu_A0, mu_P1, mu_A1)
# K_par <- c(K_P0, K_A0, K_P1, K_A1)
# gam_par <- c(gam_plant, gam_animal)
# laa_par <- c(laa_P0, laa_A0, laa_P1, laa_A1)

create_mutualism_pars <- function(lac_par,
                              mu_par,
                              K_par,
                              gam_par,
                              laa_par,
                              qgain,
                              qloss,
                              lambda1,
                              M0,
                              pro) {
  testit::assert(is.numeric(lac_par))
  testit::assert(is.numeric(mu_par))
  testit::assert(is.numeric(K_par))
  testit::assert(is.numeric(gam_par))
  testit::assert(is.numeric(laa_par))
  testit::assert(is.matrix(M0))
  testit::assert(is.numeric(qgain))
  testit::assert(is.numeric(qloss))
  testit::assert(is.numeric(lambda1))
  testit::assert(is.numeric(pro))
  testit::assert(lac_par >= 0.0)
  testit::assert(mu_par >= 0.0)
  testit::assert(K_par >= 0.0)
  testit::assert(gam_par >= 0.0)
  testit::assert(laa_par >= 0.0)
  testit::assert(qgain >=0.0)
  testit::assert(qloss >=0)
  testit::assert(lambda1 >=0)
  testit::assert(pro >=0 & pro <=1 )
  list(lac_par = lac_par,
       mu_par = mu_par,
       K_par = K_par,
       gam_par = gam_par,
       laa_par = laa_par,
       qgain = qgain,
       qloss = qloss,
       lambda1 = lambda1,
       M0 = M0,
      pro = pro)}

mutualism_pars <- create_mutualism_pars(lac_par = c(0,0),
                                      mu_par = c(0,0,0,0),
                                      K_par = c(Inf,Inf,0.5,0.5),
                                      gam_par = c(0, 1),
                                      laa_par = c(0,0,2,0),
                                      qgain = 0,
                                      qloss = 0,
                                      lambda1 = 0,
                                      M0 = {set.seed(1);matrix(sample(c(0,1),40,replace = TRUE),ncol=5,nrow=8)},
                                      pro = 1)


