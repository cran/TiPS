# TiPS : Trajectories and Phylogenies Simulator

TiPS is an R package to generate trajectories and phylogenetic trees associated with a compartmental model. Trajectories are simulated using Gillespie's exact or approximate stochastic simulation algorithm, or a newly proposed mixed version of the two. Phylogenetic trees are simulated from a trajectory under a backwards-in-time approach (i.e. coalescent).


``` r
library(TiPS)
library(ape)
```

Simulating trajectories
=======================

Building the simulator
----------------------

We use the classical SIR epidemiological model to illustrate the functionning of `TiPS`.

This SIR model can be described by a system of reactions such as

```math
S+I \xrightarrow[]{\beta S I} I+I
```

and

```math
 I \xrightarrow[]{\gamma I} R
```

or by a system of differential equations such as

```math
\frac{dS}{dt} = -\beta S I
```

and

```math
\frac{dI}{dt} = \beta S I - \gamma I
```


In R, with `TiPS`, this model can either be written as

``` r
reactions <- c("S [beta*S*I] -> I",
               "I [gamma*I] -> R")
```

Let's now build the simulator:

We then build the simulator that will allow us to run multiple trajectories:

``` r
sir_simu <- build_simulator(reactions)
```

This typically takes 10-15' as it involves compilation.

Defining simulations parameters
-------------------------------

To run numerical simulations, we define the initial values of the state variables,

``` r
initialStates <- c(I = 1, S = 9999, R = 0)
```

the time range of the simulations,

``` r
time <- c(0, 20)
```

and the parameters values

``` r
theta <- list(gamma = 1, beta = 2e-4)
```

For the *τ*-leap and mixed algorithms, a time step is also required:

``` r
dT <- 0.001
```

Running simulations
-------------------

In some simulations, the population size of a deme compartment may be zero before the upper time limit is reached, because of stochasticity or parameter values. In this case, the simulation is considered to have failed and is halted.

To bypass these failures, we can define the following wrapper:

``` r
safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}
```

A safe version of our simulator `sir_simu()` is then:

``` r
safe_sir_simu <- function(...) safe_run(sir_simu, ...)
```

### Direct method

A trajectory using Gillespie's direct method is obtained by

``` r
traj_dm <- safe_sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "exact")
```

The output consists of a named list containing the reactions of the model (with `$reactions`), the parameter values (with `$values`), the time range (with `$times`), the algorithm used to simulate (with `$algo`), the time-step in case the algorithm is *τ*-leap or the mixed algorithm (with `$dT`) and finally the simulated trajectory (with `$traj`) :

``` r
names(traj_dm)
#> [1] "reactions" "values"    "times"     "method"    "tau"       "seed"     
#> [7] "traj"
```

The simulated trajectory is also a named list, where each simulated reaction event is recorded `$Reaction`, along with the time at which it occured `$Time`, the number of times it occured `$Nrep` (if *τ*-leap or mixed algorithm chosen), and the size of each compartment througt time, here `$I` `$R` `$S`.

``` r
head(traj_dm$traj)
#>         Time          Reaction Nrep    S I R
#> 1 0.00000000              init    1 9999 1 0
#> 2 0.09966524 S [beta*S*I] -> I    1 9998 2 0
#> 3 0.13240687  I [gamma*I] -> R    1 9998 1 1
#> 4 0.63636802 S [beta*S*I] -> I    1 9997 2 1
#> 5 0.76161520 S [beta*S*I] -> I    1 9996 3 1
#> 6 0.76990828 S [beta*S*I] -> I    1 9995 4 1
```

The trajectory can readily be plotted using the `plot()` function:

``` r
plot(traj_dm)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-13-1.png" width="576" style="display: block; margin: auto;" />

You can also specify the seed to generate reproducible results with the parameter seed. By default, its value is null (set to 0), in which case the seed is randomly generated. Let's fix the seed to 166:

``` r
traj_dm <- sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "exact",
  seed = 166)
```

When running this multiple times, you will always get the exact same results:

``` r
head(traj_dm$traj)
#>         Time          Reaction Nrep    S I R
#> 1 0.00000000              init    1 9999 1 0
#> 2 0.06279117 S [beta*S*I] -> I    1 9998 2 0
#> 3 0.10048223 S [beta*S*I] -> I    1 9997 3 0
#> 4 0.13123844  I [gamma*I] -> R    1 9997 2 1
#> 5 0.22783501  I [gamma*I] -> R    1 9997 1 2
#> 6 0.35976195 S [beta*S*I] -> I    1 9996 2 2
```

### *τ*-leap method

The default mode of the simulator is the *τ*-leap method. The `method` argument must be specified and fixed to `àpproximate`. The time-step `tau` is by default set to 0.05. Let's fix its value to 0.009:

``` r
traj_tl <- safe_sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "approximate",
  tau = 0.009)
```

We obtain the same type of output as with the direct method:

``` r
head(traj_tl$traj)
#>    Time          Reaction Nrep    S I R
#> 1 0.000              init    1 9999 1 0
#> 2 0.342 S [beta*S*I] -> I    1 9998 2 0
#> 3 0.648 S [beta*S*I] -> I    1 9997 3 0
#> 4 0.963 S [beta*S*I] -> I    1 9996 4 0
#> 5 0.990  I [gamma*I] -> R    1 9996 3 1
#> 6 1.008 S [beta*S*I] -> I    1 9995 4 1
```

The trajectory can also be plotted:

``` r
plot(traj_tl)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-16-1.png" width="576" style="display: block; margin: auto;" />

### Mixed method

To run simulations with the `mixed` algorithm (basically switching between the direct method and the *τ*-leap method depending on the number of reactions occurring per unit of time), the `method` argument must be specified and fixed to `mixed`. The time-step `tau` is by default set to 0.05. Let's fix its value to 0.009 :

``` r
traj_mm <- safe_sir_simu(
    paramValues = theta,
    initialStates = initialStates,
    times = time,
    method = "mixed",
    tau = 0.009)
```

Outputs are similar to the other methods and the trajectory can also be plotted :

``` r
names(traj_mm)
#> [1] "reactions" "values"    "times"     "method"    "tau"       "seed"     
#> [7] "traj"
head(traj_mm$traj)
#>        Time          Reaction Nrep    S I R
#> 1 0.0000000              init    1 9999 1 0
#> 2 0.2947102 S [beta*S*I] -> I    1 9998 2 0
#> 3 0.5586645 S [beta*S*I] -> I    1 9997 3 0
#> 4 0.6503104  I [gamma*I] -> R    1 9997 2 1
#> 5 0.7677044 S [beta*S*I] -> I    1 9996 3 1
#> 6 0.8610195 S [beta*S*I] -> I    1 9995 4 1
plot(traj_mm)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-18-1.png" width="576" style="display: block; margin: auto;" />

The MSA is a new algorithm we introduce that switches from the direct method (GDA) to the tau-leap algorithm (GTA) if over n iterations (option `msaIt`, by default `msaIt = 10`) the time until the next event *δt* is below a user-defined threshold (option `msaTau`, by default `msaTau = tau/10`), and from GTA to GDA if the total number of realised events is lower than the number of possible events. When simulating trajectories with the mixed algorithm, it is possible to modify these parameters such as :

```r
traj_mm <- safe_sir_simu(
    paramValues = theta,
    initialStates = initialStates,
    times = time,
    method = "mixed",
    tau = 0.009,
    msaTau = 0.0001,
    msaIt = 20)
```

### Trajectory in output file

It is possible to specify to write the trajectory directly in a tab delimited output file with the option `outFile`. By default, the trajectory is as presented previously and `outFile` is empty.

```r
traj_dm <- sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "exact",
  seed = 166,
  outFile = "sir_traj.txt")
```

In that case, `traj_dm$traj` is `NULL` and the output file name is indicated in `traj_dm$outFile`. It is possible to visualize the simulated trajectory by reading with the `read.table` R function :

```r
trajectory <- read.table(file = "sir_traj.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
```

Simulating phylogenies
======================

A great advantage of `TiPS`, besides its compulational efficiency, is that it can generate phylogenies from the population dynamics trajectories using a coalescent approach.

For this, we may need a vector of sampling dates.

For the SIR example, these are stored here:

``` r
dates <- system.file("extdata", "SIR-dates.txt", package = "TiPS")
```

The `simulate_tree` function simulates a phylogeny from a trajectory object and using a set of sampling dates:

``` r
sir_tree <- simulate_tree(
  simuResults = traj_dm,
  dates = dates,
  deme = c("I"), # the type of individuals that contribute to the phylogeny
  sampled = c(I = 1), # the type of individuals that are sampled and their proportion of sampling
  root = "I", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 5,
  addInfos = FALSE) # additional info for each node
```

The `sampled` option can be used for labelled phylogenies (i.e. with multiple host types) but it requires specifying the proportion of each label type. The `root` option indicates the state of the individual initiating the dynnamics. See the next section for details.

The full phylogeny can be obtained (therefore neglecting the sampling dates) with the option `isFullTrajectory`.

Finally, some runs may fail to simulate a phylogeny from a trajectory for stochastic reasons. The `nTrials` parameter indicates the number of unsuccessful trials allowed before giving up.

The simulated phylogeny can be visualised using:

``` r
ape::plot.phylo(sir_tree, cex = .5)
```

It is possible to write the tree in a output file with the option `outFile`. If a file name is given as input, by default the tree is written in a Newick format. To write the simulated tree in a Nexus format, the option `format` must be specified and fixed to nexus (option `format = nexus`).

```r
sir_tree <- simulate_tree(
  simuResults = traj_dm,
  dates = dates,
  deme = c("I"),
  sampled = c(I = 1),
  root = "I",
  isFullTrajectory = FALSE,
  nTrials = 5,
  addInfos = FALSE,
  outFile = "sir_tree.nexus",
  format = "nexus"
)
```

Multiple demes
==============

We sometimes have multiple demes, i.e. different types of individuals that contribute to the pylogeny or that can be sampled (e.g. juveniles vs. adults or acute vs. chronic infections).

We illustrate this example using an SIR model with two patches (labelled 1 and 2) and migration between these patches (at a rate *μ*).

Initialising the system
-----------------------

$$
\\frac{dI\_1}{dt} = \\beta\_1 I\_1 - \\gamma\_1 I\_1 - \\mu\_1 I\_1 + \\mu\_2 I\_2 \\\\
\\frac{dI\_2}{dt} = \\beta\_2 I\_2 - \\gamma\_2 I\_2 - \\mu\_2 I\_2 + \\mu\_1 I\_1 \\\\
$$

The associated reactions are:

``` r
reactions <- c("0 [beta1 * I1] -> I1",
               "I1 [gamma1 * I1] -> 0",
               "I1 [mu1 * I1] -> I2",
               "0 [beta2 * I2] -> I2",
               "I2 [gamma2 * I2] -> 0",
               "I2 [mu2 * I2] -> I1")
```

We then build the simulator:

``` r
bd_simu <- build_simulator(reactions)
```

The initial state variables values are

``` r
initialStates <- c(I1 = 0, I2 = 1)
```

The time range of the simulation is between 1975 and 2018:

``` r
time <- c(1975, 1998, 2018)
```

the parameters values are

``` r
theta <- list(gamma1 = c(0.2, 0.09), gamma2 = 0.1, mu1 = 0.025, mu2 = 0.021, beta1 = c(0.26,0.37), beta2 = 0.414)
```

and the time step (for the *τ*-leap and mixed algorithms) is:

``` r
dT <- 0.01
```

A safe version of the simulator `bd_simu()` is:

``` r
safe_bd_simu <- function(...) safe_run(bd_simu, ...)
```

Tau-leap trajectory simulation
------------------------------

We perform the simulations using:

``` r
trajbd_tl <- safe_bd_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "approximate",
  tau = 0.001)
```

We obtain:

``` r
head(trajbd_tl$traj)
#>       Time              Reaction Nrep I1 I2
#> 1 1975.000                  init    1  0  1
#> 2 1975.249  0 [beta2 * I2] -> I2    1  0  2
#> 3 1977.053  0 [beta2 * I2] -> I2    1  0  3
#> 4 1978.664 I2 [gamma2 * I2] -> 0    1  0  2
#> 5 1979.283 I2 [gamma2 * I2] -> 0    1  0  1
#> 6 1981.790  0 [beta2 * I2] -> I2    1  0  2
```

Graphically, we get:

``` r
plot(trajbd_tl)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-31-1.png" width="576" style="display: block; margin: auto;" />

Phylogeny simulation
--------------------

### With known sampling dates and known proportion of sampling

Instead of loading a vector, we assume we have 100 samples at 100 sampling dates between 2015 and 2018. We can generate the dates vector as:

``` r
dates_bd <- seq(from=2015, to=2018, length.out=100)
```

We then simulate a phylogeny where 20% of the sampling dates correspond to the I1 compartment, and 80% to the I2 compartment:

``` r
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  deme = c("I1", "I2"),
  sampled = c(I1 = 0.2, I2 = 0.8), # the type of individuals that are sampled and their proportion of sampling
  root = "I2", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 3,
  addInfos = TRUE) # additional info for each node
```

This is done using a coaslescence process informed by the trajectory. Therefore, each internal node of the phylogeny corresponds to a coalescence event and is associated with a label stoed in `$node.label`.

In our two-patches example, there are two types of coalesence: I2 individuals creating a new I2 individual, and I1 individuals creating a new I1 individual.

We can plot the phylogeny and color the internal nodes based on the type of coalescence.

First we generate a vector of colors for the nodes (if we find `I2` in the node label we color it in blue, otherwise in red):

``` r
inode_cols <- ifelse(grepl(x=bd_tree$node.label,pattern="I2"),"blue","red")
```

Then we plot the phylogeny:

``` r
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, align.tip.label = T)
nodelabels(pch=20,col=inode_cols)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-35-1.png" width="576" style="display: block; margin: auto;" />

### With known sampling dates, each assigned to a compartment by the user

One can give as input, sampling dates assigned to a compartment, in which case the option `sampled` is not required.

``` r
dates_bd <- seq(from=2015, to=2018, length.out=100)
dates_bd <- data.frame(Date=sample(dates_bd),Comp=c(rep("I1",20),rep("I2",80)))
head(dates_bd)
#>       Date Comp
#> 1 2017.212   I1
#> 2 2015.515   I1
#> 3 2017.606   I1
#> 4 2015.576   I1
#> 5 2017.515   I1
#> 6 2015.030   I1
```

Now let's simulate a phylogeny with sampling dates assigned to a compartment by the user.

``` r
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 3,
  addInfos = TRUE) # additional info for each node
```

We can plot the phylogeny and color the external nodes given the compartment.

``` r
tips_cols <- ifelse(grepl(x=bd_tree$tip.label,pattern="I2"),"blue","red")
```

``` r
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)
tiplabels(pch=20,col=tips_cols)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-39-1.png" width="576" style="display: block; margin: auto;" />

### With only known sampling dates

In the case where the user has no information on the sampling proportions or the assignment of sampling dates on any compartment, the algorithm will randomly assign each sampling date to a compartment. The user gives as input sampling dates:

``` r
dates_bd <- seq(from=2015, to=2018, length.out=100)
```

Now let's simulate a phylogeny with sampling dates and no information about the sampling schemes :

``` r
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  sampled = c("I1","I2"),
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 10,
  addInfos = TRUE) # additional info for each node
```

``` r
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-42-1.png" width="576" style="display: block; margin: auto;" />

### Without given sampling dates

Let's simuate a phylogeny where simulated deaths in the trajectory generate leaves. This can be done with the `isFullTrajectory` option.

``` r
bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  deme = c("I1", "I2"),
  root = "I2", # type of individual at the root of the tree
  nTrials = 10,
  isFullTrajectory = TRUE, # deads generate leaves
  addInfos = TRUE) # additional info for each node
```

``` r
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, show.tip.label = F)
```

<img src="vignettes/TiPS_files/figure-markdown_github/unnamed-chunk-44-1.png" width="576" style="display: block; margin: auto;" />

How to cite this work
==========

Danesh G, Saulnier E, Gascuel O, Choisy M, Alizon S. 2022. TiPS: Rapidly simulating trajectories and phylogenies from compartmental models. *Methods in Ecology and Evolution*, [10.1111/2041-210X.14038](https://doi.org/10.1111/2041-210X.14038).
