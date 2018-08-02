library(devtools)
devtools::load_all()
# library(rgl)

N <- 100
T <- 1000
d <- 2

message(" computing LAIS ...")
set.seed(1)
lais_chain <- lais(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10), sigma = rep(sqrt(3), d), lposterior_6, sigma_mh = rep(1, d), compute_denom = compute_denomtable_byrow)
estim_lais <- compute_expectation(lais_chain$x, lais_chain$w)

message(" computing indep PMC ...")
set.seed(1)
pmc_idchain <- indep_pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
estim_idpmc <- compute_expectation(pmc_idchain$x, pmc_idchain$w)

message(" computing dep PMC ...")
set.seed(1)
pmc_dchain <- pmc(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(13), d), lposterior_6, compute_denom = compute_denomtable_byrow, reuse_weights = FALSE)
estim_dpmc <- compute_expectation(pmc_dchain$x, pmc_dchain$w)

message(" computing indep APIS ...")
set.seed(1)
apis_idchain <- indep_apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow)
estim_idapis <- compute_expectation(apis_idchain$x, apis_idchain$w)

message(" computing dep APIS ...")
set.seed(1)
apis_dchain <- apis(d = d, N = N, T = T, M = 2, mu = rnorm(N * d, sd = 10) , sigma = rep(sqrt(3), d), lposterior_6, compute_denom = compute_denomtable_byrow)
estim_dapis <- compute_expectation(apis_dchain$x, apis_dchain$w)


message(" Summerizing results ...")
res_table <- matrix(c(estim_lais, estim_idpmc, estim_dpmc, estim_idapis, estim_dapis),
                    nrow = d, ncol = 5)
colnames(res_table) <- c("lais", "idpmc", "dpmc", "idapis", "dapis")


#theoretical mu = c(0, -1.09)
theo <- c(-1.09, 0)
diff_table <- res_table - theo
rmse <- apply(diff_table, 2, function(x) sqrt(mean(x^2)))

library(rgl)
rgl::plot3d(lais_chain$mu[, 1,], lais_chain$mu[, 2,], 0)
rgl::plot3d(lais_chain$x[, 1,,], lais_chain$x[, 2,,], exp(lais_chain$pi))
rgl::plot3d(lais_chain$x[, 1,,], lais_chain$x[, 2,,], lais_chain$w)
