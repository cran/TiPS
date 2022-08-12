test_that("structure of phylogeny (number of leaves, nodes, heights, ...) is as expected" , {

	reactions <- c("S [beta*S*I] -> I", "I [gamma*I] -> R")
	sir_simu <- suppressMessages(build_simulator(reactions))
	initialStates <- c(I = 1, S = 9999, R = 0)
	time <- c(0, 20)
	theta <- list(gamma = 1, beta = 2e-4)
	dT <- 0.001
	suppressMessages(traj_dm <- suppressMessages(sir_simu(paramValues = theta,
  						initialStates = initialStates,
  						times = time,
  						seed = 166,
  						verbose=FALSE,
						method="exact"))
	)
	dates <- system.file("extdata", "SIR-dates.txt", package = "TiPS")
	dates_data <- read.table(dates,header=T,stringsAsFactors=F)

	sir_tree <- simulate_tree(simuResults = traj_dm,
  							  dates = dates,
  							  deme = c("I"),
  							  sampled = c(I = 1), 
  							  root = "I",
  							  isFullTrajectory = FALSE, 
  							  nTrials = 5,
  							  addInfos = TRUE,
  							  verbose=FALSE)

	## Checks the number of tips
	expect_equal(length(sir_tree$tip.label), nrow(dates_data))

	## Checks the tips height : should be identical to maximum/min dates value
	expect_equal(max(sir_tree$tip.height), max(dates_data$Dates))
	expect_equal(min(sir_tree$tip.height), min(dates_data$Dates))

	## Checks the number of nodes (ntip -1)
	expect_equal(sir_tree$Nnode, (nrow(dates_data)-1))

})
