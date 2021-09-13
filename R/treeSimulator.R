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
#' 
#' @return An object of class \code{ape::phylo}. 
#' 
#' @author Gonche Danesh
#' @export
#' 
#' @example 
#' # A birth-death model with birth rate beta, death rate gamma and I the number of infected individuals.
#' # With parameter beta varying over 10 time intervals.
#' reactions <- c("0 [beta*I] -> I", "I [gamma*I] -> 0")
#' BD.simu <- build_simulator(reactions = reactions)
#' theta <- tibble::tribble(
#'      ~beta,
#'     .8827628,
#'     .714422,
#'     .7390572,
#'     .5871399,
#'     .5542929,
#'     .6323045,
#'     .8334923,
#'     .5707164,
#'     .51734,
#'     .7308455) %$%
#' list(beta = beta, gamma = .410587333412841)
#' times <- seq(1926.17, by = 8.683, le = 11)
#' BDres <- BD.simu(paramValues = theta,
#'                   initialStates = c(I = 1),
#'                   times = times,
#'                   dT = 0.08,
#'                   nTrials = 10)
#' # File with 157 sampling dates from 1987 to 2012.
#' dates <- system.file("extdata", "BD-10_dates.txt", package = "TiPS")
#' BD_tree <- simulate_tree(simuResults = BDres,
#'                          dates = dates,
#'                          deme = c("I"),
#'                          sampled = c(I=1),
#'                          root = "I",
#'                          isFullTrajectory = FALSE,
#'                          nTrials = 5)
#' # Plot the simulated phylogeny using the \code{ape::plot.phylo} function. 
#' ape::plot.phylo(BD_tree)
simulate_tree <- function(simuResults, dates, deme, sampled, root, isFullTrajectory = FALSE, nTrials = 1, addInfos = FALSE, resampling = FALSE, verbose=FALSE){

	checkArgs<-TRUE
	if(missing(simuResults)){
		warning("Simulation results object necessary to simulate tree.")
		checkArgs<-FALSE
	} else if(length(simuResults$traj) == 0){
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
		trajectory <- simuResults$traj
	    
		rtrajectory = .merge_trajectory_dates(reactions, trajectory, dates, sampled, root, deme, isFullTrajectory)
		Phylo=Phyloepid$new(rtrajectory$reactions,rtrajectory$trajectory,isFullTrajectory,nTrials,resampling,rtrajectory$nbdates, verbose)
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
	      	class(tree)<-"phylo"

	    }else{
	    	message("Failure.")
	    }
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
	rates<-unlist(stringr::str_extract_all(string = reactions, pattern = "(?<=\\[).*(?=\\])"))
	froms<-vector()
	tos<-vector()
	for(i in 1:lnR){
	  first<-strsplit(x=reactions[[i]],split="\\[")[[1]]
	  second<-strsplit(x=reactions[[i]],split="->")[[1]]
	  froms[i]<-stringr::str_extract_all(first[1], stringr::boundary("word"))[[1]]
	  tos[i]<-stringr::str_extract_all(second[2], stringr::boundary("word"))[[1]]
	}
	indivNames<-unique(c(froms,tos))
	indivNames<-indivNames[!grepl("^[[:digit:]]+$",indivNames)]

	non_deme<-indivNames[indivNames %notin% deme]
	rates_sep<-sapply(rates, stringr::str_extract_all, stringr::boundary("word"))
	matches <- lapply(rates_sep , match  , indivNames )
	matches <- lapply(matches, function(x) x[!is.na(x)])
	
	tree_reactions<-vector()
	for(i in 1:lnR){
	  tmp<-matches[[i]]
	  indx<-which(indivNames[tmp] %notin% froms[i])
	  if(length(indx)>0){ 
	    ## if birth reaction 
	    indiv_add<-indivNames[tmp][indx]
	    if(froms[i] == "0"){
	      reac<-paste0(indiv_add,"+=[",rates[i],"]",tos[i],"+",indiv_add)
	    }else{
	      reac<-paste0(froms[1],"+",indiv_add,"+=[",rates[i],"]",tos[i],"+",indiv_add)
	    }
	  } 
	  else{
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
	for(i in 1:lnR){
	  trajectory$Reaction[trajectory$Reaction == reactions[i]] <- tree_reactions[i]
	}
	#trajectory[non_deme]<-NULL
	
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
