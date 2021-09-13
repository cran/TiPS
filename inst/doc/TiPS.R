## ---- include = FALSE---------------------------------------------------------
knitr::knit_hooks$set(
  margin = function(before, options, envir) {
    if (before) par(mgp = c(1.5, .5, 0), bty = "n", plt = c(.105, .97, .13, .97))
    else NULL
  },
  prompt = function(before, options, envir) {
    options(prompt = if (options$engine %in% c("sh", "bash")) "$ " else "> ")
  })

knitr::opts_chunk$set(
  margin = TRUE,
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  message = FALSE,
  warning = FALSE,
  dev.args = list(pointsize = 11),
  fig.height = 3.5,
  fig.width = 6,
  fig.retina = 2,
  fig.align = "center"
)

#options(width = 137)

## ----setup--------------------------------------------------------------------
library(TiPS)
library(ape)

## -----------------------------------------------------------------------------
reactions <- c("S [beta*S*I] -> I",
               "I [gamma*I] -> R")

## ----results = "hide"---------------------------------------------------------
sir_simu <- build_simulator(reactions)

## -----------------------------------------------------------------------------
initialStates <- c(I = 1, S = 9999, R = 0)

## -----------------------------------------------------------------------------
time <- c(0, 20)

## -----------------------------------------------------------------------------
theta <- list(gamma = 1, beta = 2e-4)

## -----------------------------------------------------------------------------
dT <- 0.001

## -----------------------------------------------------------------------------
safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}

## -----------------------------------------------------------------------------
safe_sir_simu <- function(...) safe_run(sir_simu, ...)

## ----results = "hide"---------------------------------------------------------
traj_dm <- safe_sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time)

## -----------------------------------------------------------------------------
names(traj_dm)

## -----------------------------------------------------------------------------
head(traj_dm$traj)

## -----------------------------------------------------------------------------
plot(traj_dm)

## ----results = "hide"---------------------------------------------------------
traj_tl <- safe_sir_simu(
  paramValues = theta,
	initialStates = initialStates,
	times = time,
  method = "approximate",
	tau = 0.009)

## -----------------------------------------------------------------------------
head(traj_tl$traj)

## -----------------------------------------------------------------------------
plot(traj_tl)

## ----results = "hide"---------------------------------------------------------
traj_mm <- safe_sir_simu(
	paramValues = theta,
	initialStates = initialStates,
	times = time,
	method = "mixed",
	tau = 0.009)

## -----------------------------------------------------------------------------
names(traj_mm)
head(traj_mm$traj)
plot(traj_mm)

## -----------------------------------------------------------------------------
dates <- system.file("extdata", "SIR-dates.txt", package = "TiPS")

## ---- results = "hide"--------------------------------------------------------
sir_tree <- simulate_tree(
  simuResults = traj_dm,
  dates = dates,
  deme = c("I"), # the type of individuals that contribute to the phylogeny
  sampled = c(I = 1), # the type of individuals that are sampled and their proportion of sampling
  root = "I", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 5,
  addInfos = FALSE) # additional info for each node

## ----fig.height = 10----------------------------------------------------------
ape::plot.phylo(sir_tree, cex = .5)

## -----------------------------------------------------------------------------
reactions <- c("0 [beta1 * I1] -> I1",
               "I1 [gamma1 * I1] -> 0",
               "I1 [mu1 * I1] -> I2",
               "0 [beta2 * I2] -> I2",
               "I2 [gamma2 * I2] -> 0",
               "I2 [mu2 * I2] -> I1")

## ----results = "hide"---------------------------------------------------------
bd_simu <- build_simulator(reactions)

## -----------------------------------------------------------------------------
initialStates <- c(I1 = 0, I2 = 1)

## -----------------------------------------------------------------------------
time <- c(1975, 1998, 2018)

## -----------------------------------------------------------------------------
theta <- list(gamma1 = c(0.2, 0.09), gamma2 = 0.1, mu1 = 0.025, mu2 = 0.021, beta1 = c(0.26,0.37), beta2 = 0.414)

## -----------------------------------------------------------------------------
dT <- 0.01

## -----------------------------------------------------------------------------
safe_bd_simu <- function(...) safe_run(bd_simu, ...)

## ----results = "hide"---------------------------------------------------------
trajbd_tl <- safe_bd_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "approximate",
  tau = 0.001)

## -----------------------------------------------------------------------------
head(trajbd_tl$traj)

## -----------------------------------------------------------------------------
plot(trajbd_tl)

## ----results = "hide"---------------------------------------------------------
dates_bd <- seq(from=2015, to=2018, length.out=100)

## ----results = "hide"---------------------------------------------------------
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  deme = c("I1", "I2"),
  sampled = c(I1 = 0.2, I2 = 0.8), # the type of individuals that are sampled and their proportion of sampling
  root = "I2", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 3,
  addInfos = TRUE) # additional info for each node

## ----results = "hide"---------------------------------------------------------
inode_cols <- ifelse(grepl(x=bd_tree$node.label,pattern="I2"),"blue","red")

## ----fig.height = 10----------------------------------------------------------
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, align.tip.label = T)
nodelabels(pch=20,col=inode_cols)

## -----------------------------------------------------------------------------
dates_bd <- seq(from=2015, to=2018, length.out=100)
dates_bd <- data.frame(Date=sample(dates_bd),Comp=c(rep("I1",20),rep("I2",80)))
head(dates_bd)

## ----results = "hide"---------------------------------------------------------
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 3,
  addInfos = TRUE) # additional info for each node

## ----results = "hide"---------------------------------------------------------
tips_cols <- ifelse(grepl(x=bd_tree$tip.label,pattern="I2"),"blue","red")

## ----fig.height = 10----------------------------------------------------------
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)
tiplabels(pch=20,col=tips_cols)

## -----------------------------------------------------------------------------
dates_bd <- seq(from=2015, to=2018, length.out=100)

## ----results = "hide"---------------------------------------------------------
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  sampled = c("I1","I2"),
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 10,
  addInfos = TRUE) # additional info for each node

## ----fig.height = 10----------------------------------------------------------
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)

## ----results = "hide"---------------------------------------------------------
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 10,
  isFullTrajectory = TRUE, # deads generate leaves
  addInfos = TRUE) # additional info for each node

## ----fig.height = 10----------------------------------------------------------
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)

