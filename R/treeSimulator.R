#' Simulates a phylogeny using a beackward-in-time process using sampling dates and a trajectory
#'
#' @param simuResults Object of class \code{simutraj} resulting from running a simulator of trajectories built using the \code{build_simulator} function.
#' @param dates Contains the sampling dates. Can be a vector (for example using \code{seq} function), can be a named list or a file with header.
#' @param deme A vector with the compartments that contribute to the simulation of phylogeny.
#' @param sampled A named vector with the proportions of sampling for each compartment. This is used in case there are multiple deme compartments where the sampling dates will be randomly associated to a compartment to sample. Sum of \code{sampled} must be 1.
#' @param root Name of the compartment from which the phylogeny is rooted.
#' @param isFullTrajectory Boolean to define if death events generate or not leaves. By default, \code{isFullTrajectory=FALSE}.
#' @param addInfos Boolean to define if each internal node' name will be the reaction. By default, \code{addInfos=FALSE}.
#' @param resampling Boolean to allow a sampled individual to transmit the pathogen once again. By default, \code{resampling=FALSE}.
#' @param seed Seed to initialize the random generator, for better reproducibility. By default, \code{seed=0} and the seed value is randomly generated.
#' @param outFile Output file name to write tree. By default, tree is not written in output file.
#' @param format Output tree format if output file is given. Values are either \code{format = "newick"} ou \code{format = "nexus"}. By default, \code{format = "newick"}
#' @param nTrials Integer indicating the number of unsuccessful simulation trials allowed before giving up.
#' @param verbose Boolean to allow printing time execution of simulation
#'
#' @returns An object of class \code{ape::phylo}.
#'
#' @author Gonche Danesh
#' @export
#'
#' @examples
#' \dontrun{
#' # A multi-type birth-death model with birth rate beta,
#' # death rate gamma, mutation rates m1 and m2
#' # and I1 and I2 the number of infected individuals of each type.
#' # With parameter beta varying over 2 time intervals.
#' reactions <- c("0 [beta1 * I1] -> I1","I1 [gamma1 * I1] -> 0",
#' "I1 [mu1 * I1] -> I2","0 [beta2 * I2] -> I2",
#' "I2 [gamma2 * I2] -> 0","I2 [mu2 * I2] -> I1")
#'
#' BD_simu <- build_simulator(reactions)
#' initialStates <- c(I1 = 0, I2 = 1)
#' times <- c(1975, 1998, 2018)
#' theta <- list(gamma1 = c(0.2, 0.09), gamma2 = 0.1, mu1 = 0.025,
#' mu2 = 0.021, beta1 = c(0.26,0.37), beta2 = 0.414)
#' BDres <- BD_simu(paramValues = theta,
#'                   initialStates = initialStates,
#'                   times = times,
#'                   tau = 0.08,
#'                   method = "approximate",
#'									 seed = 994543)
#' # Let's generate 100 sampling dates from 2015 and 2018
#' dates_bd <- seq(from=2015, to=2018, length.out=100)
#' BD_tree <- simulate_tree(simuResults = BDres,
#'                          dates = dates,
#'                          deme = c("I"),
#'                          sampled = c(I=1),
#'                          root = "I",
#'                          isFullTrajectory = FALSE,
#'                          seed = 973360)
#' BD_tree$seed
#' # Plot the simulated phylogeny using the \code{ape::plot.phylo} function.
#' ape::plot.phylo(BD_tree)
#' }
simulate_tree <- function(simuResults, dates, deme, sampled, root, isFullTrajectory = FALSE, nTrials = 1, addInfos = FALSE, resampling = FALSE, verbose=FALSE, seed=0, outFile="", format="newick"){

	checkArgs<-TRUE
	if(missing(simuResults)){
		warning("Simulation results object necessary to simulate tree.")
		checkArgs<-FALSE
	} else if(length(simuResults$traj) == 0 & length(simuResults$outFile) == 0){
		warning("The trajectory in simuResults argument is an empty list.")
		checkArgs<-FALSE
	}
	if(missing(deme)){
		warning("Vector of deme(s) not given.")
		checkArgs<-FALSE
	}
	if(missing(dates) & missing(isFullTrajectory)){
		warning("No dates given, and no information if simulation of full phylogeny.")
		checkArgs<-FALSE
	}
	if(missing(root)){
	    warning("Name of the individual at the origin of the rooting of the tree not given.")
	    checkArgs<-FALSE
	}
	if(missing(dates) & isFullTrajectory){
	    isFullTrajectory = TRUE
	    dates = NULL
	}

	if(checkArgs){

		reactions <- simuResults$reactions
		warning(simuResults$outFile)
		if(length(simuResults$outFile) != 0 ){#& utils::file_test(op="-f", x=simuResults$outFile)){
			trajectory <- utils::read.table(file = simuResults$outFile, header = T, stringsAsFactors = F, sep = "\t")
		} else{
			trajectory <- simuResults$traj
		}

		options <- list("seed" = seed, "nTrials" = nTrials, "deme" = deme)

		rtrajectory = .merge_trajectory_dates(reactions, trajectory, dates, sampled, root, deme, isFullTrajectory)
		Phylo=Phyloepid$new(rtrajectory$reactions,rtrajectory$trajectory,isFullTrajectory,resampling,rtrajectory$nbdates, verbose, options)
		ok=Phylo$simulationTree()

		tree<-list()

	    if(ok){
	    	message("Success !")
	      	tree<-list()
	      	tree=Phylo$createTreeObject()
	      	edges=cbind(tree$from,tree$to)
	      	tree[["edge"]]=edges
	      	tree$from=NULL
	      	tree$to=NULL
	      	if(!addInfos){
	      		tree$node.label=NULL
	      		tree$tip.height=NULL
	      	}
					if(outFile != ""){
						if(format == "newick"){
							tree_str <- Phylo$getNewickTree(addInfos)
						}
						else if(format == "nexus"){
							tree_str <- Phylo$getNexusTree(addInfos)
						}
						write(x = tree_str, file = outFile)
					}
	      	class(tree)<-"phylo"

	    }else{
	    	message("Failure.")
	    }
	    rm(trajectory)
	    rm(Phylo)
	    rm(ok)
	    return(tree)
    }

}

`%notin%` <- Negate(`%in%`)

#' Merges simulated trajectory with the sampling dates chronologically
#' @param reactions a vector of reactions
#' @param trajectory simulated trajectory
#' @param dates a vector of sampling dates
#' @param sampled list with the proportion of sampling for each deme compartment
#' @param root compartment rooting the tree
#' @param deme vector of deme compartments
#' @param isFullTrajectory simulate a complete tree
#' @NoRd
.merge_trajectory_dates <- function(reactions, trajectory, dates, sampled, root, deme, isFullTrajectory){

	## From traj simulation formalism to tree simulation formalism
	lnR<-length(reactions)
	rates<-unlist(stringr::str_extract_all(string = reactions, pattern = "(?<=\\[).*(?=\\])")) # récupérer ce qui est entre crochets
	froms<-vector()
	tos<-vector()
	for(i in 1:lnR){
	  first<-strsplit(x=reactions[[i]],split="\\[")[[1]]
	  second<-strsplit(x=reactions[[i]],split="->")[[1]]
	  froms[i]<-stringr::str_extract_all(first[1], stringr::boundary("word"))[[1]]
	  tos[i]<-stringr::str_extract_all(second[2], stringr::boundary("word"))[[1]]
	}
	indivNames<-unique(c(froms,tos))
	indivNames<-indivNames[!grepl("^[[:digit:]]+$",indivNames)] # retirer le 0 dans le cas où naissance sans individu
	non_deme<-indivNames[indivNames %notin% deme]

	rates_sep<-sapply(rates, stringr::str_extract_all, stringr::boundary("word"))
	matches <- lapply(rates_sep , match  , indivNames )
	matches <- lapply(matches, function(x) x[!is.na(x)])

	tree_reactions<-vector(length = lnR)
	for(i in 1:lnR){
	  tmp<-matches[[i]]
	  indx<-which(indivNames[tmp] %notin% froms[i])
	  if(length(indx)>0){
			## if birth reaction
			if(length(indx)>1){
				## if for example rate = S1*I1/(I1+I2+S1...)
				indiv_add<-indivNames[tmp][indx][1]
			} else{
				indiv_add<-indivNames[tmp][indx]
			}
			if(froms[i] == "0"){
				reac<-paste0(indiv_add,"+=[",rates[i],"]",tos[i],"+",indiv_add)
			}else{
				reac<-paste0(froms[i],"+",indiv_add,"+=[",rates[i],"]",tos[i],"+",indiv_add)
			}
	  } else{
	    if(tos[i] == "0"){
	      ## if death
	      reac<-paste0(froms[i],"-=[",rates[i],"]")
	    } else{
	      ## else, migration
	      reac<-paste0(froms[i],"=[",rates[i],"]",tos[i])
	    }
	  }
	  tree_reactions[i]<-reac
	}
	#for(i in 1:lnR){
	#  trajectory$Reaction[trajectory$Reaction == reactions[i]] <- tree_reactions[i]
	#}
	#trajectory[non_deme]<-NULL
	trajectory$tmp <- NA
	for(i in 1:lnR){
	  trajectory$tmp[which(trajectory$Reaction == reactions[i])] <- tree_reactions[i]
	}
	trajectory$tmp[which(trajectory$Reaction == "init")] <- "init"
	trajectory$Reaction <- trajectory$tmp
	trajectory$tmp <- NULL

	if(sum(!indivNames %in% deme) > 0){
		non_deme_reactions<-vector()
		deme_indx_to_search<-which(indivNames %in% deme)
		indx_non_deme_reactions<-which(names(rates_sep) %notin% names(which(sapply(matches, function(r) any(r %in% deme_indx_to_search)))))
		non_deme_reactions<-tree_reactions[indx_non_deme_reactions]
		trajectory[trajectory$Reaction %notin% non_deme_reactions,]
	}
	reactions<-tree_reactions
	############################
	if(isFullTrajectory){
		 deme_death_reactions<-tree_reactions[stringr::str_detect(tree_reactions,pattern="-=")]
		 tmp<-trajectory[which(trajectory$Reaction %in% deme_death_reactions),]
		 lndates<-sum(tmp$Nrep)
	}

	reactions<-c(reactions,paste('+=',root,sep=""))
	trajectory$Reaction[1] <- paste('+=',root,sep="") ## add rooting reaction

	#### Sampled issue
	## Get dates
	# If file is given as input, it must have a header
	if(!is.null(dates)){
		cols<-c("Time","Reaction")
		if(is.list(dates)){
		 	len<-length(colnames(dates))
		  	colnames(dates)[1:len]<-cols[1:len]
		  	sampled_nouns<-as.vector(unique(dates$Reaction))
		  	n<-length(sampled_nouns)
		  	sampled<-vector(mode = "list", length = n)
		  	names(sampled)<-sampled_nouns
		} else if(length(dates) == 1){
			if(file.exists(dates)){
		    	tmp<-read.table(file = dates, header = T, stringsAsFactors = F)
		      	len<-length(colnames(tmp))
		      	colnames(tmp)[1:len]<-cols[1:len]
		      	dates<-tmp
		    } else{
		      	ok = FALSE
		      	warning("Error : input of sampling dates is wrong.")
		    }
		} else if(is.vector(dates)){
		    if(length(dates) == 1){
		    	ok = FALSE
		    	warning("Error : input of sampling dates is wrong.")
		  	} else{
		    	dates<-data.frame(Time=dates, stringsAsFactors = F)
		  	}
		} else {
		  	ok = FALSE
		  	warning("Error : input of sampling dates is wrong.")
		}

		## Get sampled : if no proportion given (ex: sampled<-c('I1','I2'))
		if(is.null(names(sampled))){
		  	tmp<-as.list(gtools::rdirichlet(1,c(rep(1,length(sampled)))))
		  	names(tmp)<-sampled
		  	sampled<-tmp
		}

		## Dates + Reaction
		if(length(colnames(dates)) == 1){
		  	# Only dates -> generate randomly
		  	dates$Reaction<-sample(names(sampled), nrow(dates), prob=sampled, replace=T)
		}
		if(length(colnames(dates)) == 0){
			# Only dates -> generate randomly
			dates$Reaction<-sample(names(sampled), length(dates), prob=sampled, replace=T)
		}

		for(i in 1:length(sampled)){
			reactions<-c(reactions,names(sampled)[i])
		}
		compNames<-colnames(trajectory)[!(colnames(trajectory) %in% c("Time", "Reaction", "Nrep"))]

		dates$Nrep<-1

		lndates<-length(dates$Nrep)

		dates<-stats::aggregate(Nrep~Time+Reaction, dates, sum)

		if(sum(!(colnames(trajectory) %in% c("Time", "Reaction", "Nrep")))!=0){
			dates<-cbind(dates, as.data.frame(matrix(0,ncol=ncol(trajectory)-length(c("Time", "Reaction", "Nrep")),nrow=nrow(dates), dimnames=list(NULL,colnames(trajectory)[!(colnames(trajectory) %in% c("Time", "Reaction", "Nrep"))]))))
		}

		trajectory<-trajectory[trajectory$Time<max(dates$Time),]

		fullTraj<-rbind(trajectory, dates[,colnames(trajectory)])


	}else{
		fullTraj<-trajectory
	}


	fullTraj$Reaction<-factor(fullTraj$Reaction, levels=reactions, ordered = TRUE)

	fullTraj<-fullTraj[order(fullTraj$Time, fullTraj$Reaction, decreasing=T),colnames(trajectory)]
	fullTraj$Reaction=as.character(fullTraj$Reaction)


	results=list("trajectory" = fullTraj, "reactions" = reactions, "nbdates" = lndates)
	results

}
