## Build stochastic simulator
reactions <- c("S [beta*S*I] -> I", "I [gamma*I] -> R")
sir_simu <- suppressMessages(build_simulator(reactions))
detfile <- system.file("extdata", "sir_det.rda", package = "TiPS")
load(detfile)

test_that("80% enveloppe of simulated trajectories by GDA include deterministic trajectory",{
  ## Parameters
  initialStates <- c(I = 1, S = 9999, R = 0);time <- c(2000, 2020);theta <- c(gamma = 1, beta = 2e-4)

  cuts <- seq(from=2000.2,to=2020,by=0.200)
  detI <- vector(mode = "numeric", length = length(cuts))
  for(i in 1:length(cuts)){
    indx <- which(round(dDet001$Time,1) <= cuts[i])
    indx <- indx[length(indx)]
    detI[i] <- dDet001$I[indx]
  }
  
  ## Simulations
  nb <- 1
  direct_stats <- vector("list", length = length(cuts)) ; names(direct_stats) <- cuts
  while(nb < 1001){
    suppressMessages(traj_dm <- sir_simu(paramValues = theta,
                                         initialStates = initialStates,
                                         times = time,
                                         verbose = FALSE,
                                         nTrials = 3)
    )
    if(length(traj_dm$traj) > 0){
      for (cut in cuts) {
        indx <- which(round(traj_dm$traj$Time,1) <= cut)
        indx <- indx[length(indx)]
        tmp <- traj_dm$traj$I[indx]
        direct_stats[[as.character(cut)]] <- append(direct_stats[[as.character(cut)]], tmp)
      }
      rm(traj_dm)
      nb <- nb + 1
    }
  }
  
  ## Get 80% HDI
  range <- 0.80
  prob1=(1-range)/2
  prob2=1-prob1
  low <- as.vector(unlist(lapply(direct_stats, quantile, prob=prob1)))
  upp <- as.vector(unlist(lapply(direct_stats, quantile, prob=prob2)))
  means <- as.vector(unlist(lapply(direct_stats, mean)))
  medians <- as.vector(unlist(lapply(direct_stats, median)))
  vartraj <- data.frame(Time=seq(from = 2000.2,to = 2020,by = 0.2),
                        Upper=upp,Lower=low, Mean=means, Median=medians, Det=detI)
  
  indx <- which(vartraj$Det > vartraj$Upper | vartraj$Det < vartraj$Lower)
  
  ## Check that deterministic trajectory is included between lower and upper range of 80% of stochatic trajectories
  expect_equal(length(indx), 0)

})

test_that("80% enveloppe of simulated trajectories by GTA include deterministic trajectory",{
  ## Parameters
  initialStates <- c(I = 1, S = 9999, R = 0);time <- c(2000, 2020);theta <- c(gamma = 1, beta = 2e-4);dT <- 0.001
  
  cuts <- seq(from=2000.2,to=2020,by=0.200)
  detI <- vector(mode = "numeric", length = length(cuts))
  for(i in 1:length(cuts)){
    indx <- which(round(dDet001$Time,1) <= cuts[i])
    indx <- indx[length(indx)]
    detI[i] <- dDet001$I[indx]
  }
  
  ## Simulations
  nb <- 1
  tau_stats <- vector("list", length = length(cuts)) ; names(tau_stats) <- cuts
  while(nb < 1001){
    suppressMessages(traj_tl<- sir_simu(paramValues = theta,
                                         initialStates = initialStates,
                                         times = time,
                                         verbose = FALSE,
                                         nTrials = 3,
                                         method = "approximate",
                                         tau = dT)
    )
    if(length(traj_tl$traj) > 0){
      for (cut in cuts) {
        indx <- which(round(traj_tl$traj$Time,1) <= cut)
        indx <- indx[length(indx)]
        tmp <- traj_tl$traj$I[indx]
        tau_stats[[as.character(cut)]] <- append(tau_stats[[as.character(cut)]], tmp)
      }
      rm(traj_tl)
      nb <- nb + 1
    }
  }
  
  ## Get 80% HDI
  range <- 0.80
  prob1=(1-range)/2
  prob2=1-prob1
  low <- as.vector(unlist(lapply(tau_stats, quantile, prob=prob1)))
  upp <- as.vector(unlist(lapply(tau_stats, quantile, prob=prob2)))
  means <- as.vector(unlist(lapply(tau_stats, mean)))
  medians <- as.vector(unlist(lapply(tau_stats, median)))
  vartraj <- data.frame(Time=seq(from = 2000.2,to = 2020,by = 0.2),
                        Upper=upp,Lower=low, Mean=means, Median=medians, Det=detI)
  
  indx <- which(vartraj$Det > vartraj$Upper | vartraj$Det < vartraj$Lower)
  
  ## Check that deterministic trajectory is included between lower and upper range of 80% of stochatic trajectories
  expect_equal(length(indx), 0)
  
})

test_that("80% enveloppe of simulated trajectories by MSA include deterministic trajectory",{
  ## Parameters
  initialStates <- c(I = 1, S = 9999, R = 0);time <- c(2000, 2020);theta <- c(gamma = 1, beta = 2e-4);dT <- 0.001

  cuts <- seq(from=2000.2,to=2020,by=0.200)
  detI <- vector(mode = "numeric", length = length(cuts))
  for(i in 1:length(cuts)){
    indx <- which(round(dDet001$Time,1) <= cuts[i])
    indx <- indx[length(indx)]
    detI[i] <- dDet001$I[indx]
  }
  
  ## Simulations
  nb <- 1
  mix_stats <- vector("list", length = length(cuts)) ; names(mix_stats) <- cuts
  while(nb < 1001){
    suppressMessages(traj_mm<- sir_simu(paramValues = theta,
                                        initialStates = initialStates,
                                        times = time,
                                        verbose = FALSE,
                                        nTrials = 3,
                                        method = "mixed",
                                        tau = dT)
    )
    if(length(traj_mm$traj) > 0){
      for (cut in cuts) {
        indx <- which(round(traj_mm$traj$Time,1) <= cut)
        indx <- indx[length(indx)]
        tmp <- traj_mm$traj$I[indx]
        mix_stats[[as.character(cut)]] <- append(mix_stats[[as.character(cut)]], tmp)
      }
      rm(traj_mm)
      nb <- nb + 1
    }
  }
  
  ## Get 80% HDI
  range <- 0.80
  prob1=(1-range)/2
  prob2=1-prob1
  low <- as.vector(unlist(lapply(mix_stats, quantile, prob=prob1)))
  upp <- as.vector(unlist(lapply(mix_stats, quantile, prob=prob2)))
  means <- as.vector(unlist(lapply(mix_stats, mean)))
  medians <- as.vector(unlist(lapply(mix_stats, median)))
  vartraj <- data.frame(Time=seq(from = 2000.2,to = 2020,by = 0.2),
                        Upper=upp,Lower=low, Mean=means, Median=medians,Det=detI)
  
  indx <- which(vartraj$Det > vartraj$Upper | vartraj$Det < vartraj$Lower)
  ## Check that deterministic trajectory is included between lower and upper range of 80% of stochatic trajectories
  expect_equal(length(indx), 0)
  
})
