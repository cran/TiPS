
#' Checks that input vectors match
#' @param params a vector of parameter values of the model
#' @param times a vetor of time limits for the simulation
#' @param indiv a vector of compartment names
#' @NoRd
.checkParamValues<-function(params,times,indiv){
  ## Check if NA values in the parameter values
  i=1
  j=1
  k=1
  ok=TRUE
  while(i < length(params) && j < length(times) && k < length(indiv)){
    if( sum(is.na(params[[i]])) || sum(is.na(times[j])) || sum(is.na(indiv[[k]])) ){
      ok=FALSE
      i=length(params)+1
    } else{
      i=i+1
      j=j+1
      k=k+1
    }
  }
  return(ok)
}

`%notin%` <- Negate(`%in%`)

#' Generates a simulator in string
#' @param reactions a character vector of reactions
#' @param functions a character vector of functions
#' @NoRd
simulator_generator<-function(reactions, functions=NULL){

  message("Building the simulator ...")


    ############### CONSTRUCTEUR ###########
  constructor_gillespie<-'
  reactions(List params, NumericVector initialStates, NumericVector times, double deltaT, int nTrials, bool isSwitchingMode, bool verbose, long seed, double msaTau, double msaIt, string outFile) : nTrials_(nTrials), deltaT_(deltaT), isSwitchingMode_(isSwitchingMode), seed_(seed), msaTau_(msaTau), msaIt_(msaIt), outFile_(outFile){
    verbose_ = verbose;
    for (unsigned i = 0; i < params.size(); i++){
      parameters_.push_back(params[i]);
    }
    initialStates_ = initialStates;
    vTimes_ = times;
    initStructures();
    if(seed_ == 0){
      initRandomSeed();
    }
    randomGenerator_.seed(seed_);
    time_ = vTimes_[0];
    isOutFile_ = !outFile_.empty();
  }
  '

  ############ RUN SWITCHING MODE ##############
  runSwitch<-'
  bool runWithWitchingMode(){
    unsigned timePart = 0;
    time_ = vTimes_[timePart];
    int timePartLim = (vTimes_.size()-1);
    initParams(timePart);
    double jump = 1;
    double totalRate = 1;
    unsigned switchingThreshold = 0;
    int nEvents = 0;
    ofstream outTrajectory;
    if(!isOutFile_) {
      updateDataFrameOutput();
    }
    else {
      outTrajectory.open(outFile_);
      outTrajectory << trajectoryHeaderToString() << endl;
      outTrajectory << "" << setprecision(10) << time_ << "\\t" ;
      outTrajectory << compartmentStatesToString() << "\\t" ;
      outTrajectory << "init\\t1" << endl;;
    }
    for ( ; timePart<timePartLim && computeTotalRate()!=0 && nEvents>=0 ; timePart++) {
      initParams(timePart);
      totalRate = computeTotalRate();
      jump = rexp(totalRate);
      while(continueCondition(timePart) && nEvents>=0 ){
        if(jump < msaTau_){
          switchingThreshold ++;
        }
        else{
          switchingThreshold = 0;
        }
        if(switchingThreshold > msaIt_ || ((double)nEvents/(double)vreactions.size()) > 1){
          nEvents = tauLeapAlgo(outTrajectory);
          jump = deltaT_;
          switchingThreshold = 0;
        }
        else{
          jump = directAlgo(outTrajectory);
          nEvents =1;
        }
      }
    }
    if (isOutFile_) {
      outTrajectory.close();
    }
    return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
  }
  '

  ######### RUN DIRECT ##########
  runDirect_gillespie<-'
  bool runDirect(){
    unsigned timePart = 0;
    double jump = 1;
    time_ = vTimes_[timePart];
    int timePartLim = (vTimes_.size()-1);
    initParams(timePart);
    ofstream outTrajectory;
    if(!isOutFile_) {
      updateDataFrameOutput();
    }
    else {
      outTrajectory.open(outFile_);
      outTrajectory << trajectoryHeaderToString() << endl;
      outTrajectory << "" << setprecision(10) << time_ << "\\t" ;
      outTrajectory << compartmentStatesToString() << "\\t" ;
      outTrajectory << "\\tinit\\t1" << endl;;
    }
    for ( ; timePart<timePartLim && computeTotalRate()!=0 ; timePart++) {
      initParams(timePart);
      while(continueCondition(timePart)){
        jump = directAlgo(outTrajectory);
      }
    }
    if (isOutFile_) {
      outTrajectory.close();
    }
    return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
  }
  '

  ############# DIRECT ALGORITHM ###########
  directAlgo_gilespie<-'
  double directAlgo(ofstream& outTrajectory){
    double totalRate = 0.0;
    double jump = 0.0;
    int reactIndx = 0;
    totalRate = computeTotalRate();
    jump = rexp(totalRate);
    time_ += jump;
    reactIndx = whichReaction(totalRate);
    performReaction(reactIndx,1);
    updateCompartments();
    outTrajectory << "" << setprecision(10) << time_ << "\\t" ;
    updateOutput(reactIndx, 1, outTrajectory);
    updateRates();
    return jump;
  }

  '

  ############ GILLESPIE SIMULATION ###########
  gillespieSimulation_gil<-'
  List GillespieSimulation(){
    if(verbose_){
      Rcout << "Running simulation of the trajectory..." << endl;
    }
    bool ok = false;
    auto start_ = std::chrono::high_resolution_clock::now();
    auto stop_ = std::chrono::high_resolution_clock::now();
    if (deltaT_ > 0){
      if(isSwitchingMode_){
        unsigned i = 0;
        if (nTrials_ > 1){
          if(verbose_){
            Rcout << "- Trial " << (i+1) << "..." << endl;
          }
        }
        start_ = std::chrono::high_resolution_clock::now();
        ok = runWithWitchingMode();
        stop_ = std::chrono::high_resolution_clock::now();
        for (i = 1; i < nTrials_ && !ok; i++){
          reInitializeDataFrame();
          if(verbose_){
            Rcout << "- Trial " << (i+1) << "..." << endl;
          }
          start_ = std::chrono::high_resolution_clock::now();
          ok = runWithWitchingMode();
          stop_ = std::chrono::high_resolution_clock::now();
        }
      }
      else{
        unsigned i = 0;
        if (nTrials_ > 1){
          if(verbose_){
            Rcout << "- Trial " << (i+1) << "..." << endl;
          }
        }
        start_ = std::chrono::high_resolution_clock::now();
        ok = runTau();
        stop_ = std::chrono::high_resolution_clock::now();
        for (i = 1; i < nTrials_ && !ok; i++){
          reInitializeDataFrame();
          if(verbose_){
            Rcout << "- Trial " << (i+1) << "..." << endl;
          }
          start_ = std::chrono::high_resolution_clock::now();
          ok = runTau();
          stop_ = std::chrono::high_resolution_clock::now();
        }
      }
    }
    else{
      unsigned i = 0;
      if (nTrials_ > 1){
        if(verbose_){
          Rcout << "- Trial " << (i+1) << "..." << endl;
        }
      }
      start_ = std::chrono::high_resolution_clock::now();
      ok = runDirect();
      stop_ = std::chrono::high_resolution_clock::now();
      for (i = 1; i < nTrials_ && !ok; i++){
        reInitializeDataFrame();
          if(verbose_){
            Rcout << "- Trial " << (i+1) << "..." << endl;
          }
        start_ = std::chrono::high_resolution_clock::now();
        ok = runDirect();
        stop_ = std::chrono::high_resolution_clock::now();
      }
    }
    if (ok){
      List res;
      res["seed"] = seed_;
  '

  ############# RUN TAU #############
  runTau<-'
  bool runTau(){
  unsigned timePart = 0;
  time_ = vTimes_[timePart];
  int timePartLim = (vTimes_.size()-1);
  initParams(timePart);
  int nEvents = 0;
  ofstream outTrajectory;
  if(!isOutFile_) {
    updateDataFrameOutput();
  }
  else {
    outTrajectory.open(outFile_);
    outTrajectory << trajectoryHeaderToString() << endl;
    outTrajectory << "" << setprecision(10) << time_ << "\\t" ;
    outTrajectory << compartmentStatesToString() << "\\t" ;
    outTrajectory << "\\tinit\\t1" << endl;;
  }
  for ( ; timePart<timePartLim && computeTotalRate()!=0 && nEvents>=0 ; timePart++) {
    initParams(timePart);
    while(continueCondition(timePart) && nEvents>=0 ){
      nEvents = tauLeapAlgo(outTrajectory);
    }
  }
  if (isOutFile_) {
    outTrajectory.close();
  }
  return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
}

  '

  ############ TAU LEAP ALGO #########
  tauLeapAlgo<-'
  int tauLeapAlgo(ofstream& outTrajectory){
  unsigned nEvents = 0;
    vector<pair<int,int> > reactionsToPerform;
    int indx = 0 ;
    vector<double>::iterator it;
    bool ok = TRUE;
    for (it = vRates_.begin(); it!=vRates_.end() && ok ; ++it,indx++) {
      unsigned nTime = 0 ;
      double rate = *it * deltaT_;
      if(rate > 0 && rate < 1000000000){
        nTime = rpois(rate) ;
        unsigned limit = 20;
        while(limit > 0 && !isFromSizeEnough(indx, nTime)){
          nTime = rpois(rate);
          limit--;
        }
        if(limit == 0 && !isFromSizeEnough(indx, nTime)){
          ok = FALSE;
          warning("Error : Stuck in Poisson distribution.You should try a smaller value for the time step tau.");
        }
      }
      if (rate > 1000000000){
        ok = FALSE;
        warning("Error: total rate too high. You should try other parameter or deltaT values.");
      }
      if (nTime >= numeric_limits<int>::max()){
        ok = FALSE;
        warning("Error: Number of events to perform has reached the maximal INT value");
      }
      if (nTime > 0){
        reactionsToPerform.push_back(make_pair(indx, nTime));
        nEvents += nTime;
      }
    }
    if(ok){
      time_ += deltaT_ ;
      for(unsigned i = 0; i<reactionsToPerform.size() ; i++){
        outTrajectory << "" << setprecision(10) << time_ << "\\t" ;
        int reac = reactionsToPerform[i].first;
        int nbTime = reactionsToPerform[i].second;
        performReaction(reac, nbTime);
        updateCompartments();
        updateOutput(reac, nbTime, outTrajectory);
        updateRates();
        if (!checkCompartmentSize()) {
          i = reactionsToPerform.size(); // to end the loop
          warning("Error : Population size negative. You should try a smaller value for the time step tau.");
          nEvents = -1;
        }
      }
    }
    else{
      nEvents = -1;
    }
    return nEvents;
  }
  '

  ################## RCPP MODULE ###############
  rcppModule_gillespie<-'
  RCPP_MODULE(reactionsmodule){
  Rcpp::class_<reactions>( "reactions" )
  .constructor<List, NumericVector, NumericVector, double, long, int, bool, int, int, string>("documentation for constructor")
  .method( "GillespieSimulation" , &reactions::GillespieSimulation, "final simulation")
  ;
  }
  '

  ############ HEADER ############
  header_gillespie<-'
  // [[Rcpp::plugins(cpp11)]]
  #include <Rcpp.h>
  #include <map>
  #include <utility>
  #include <algorithm>
  #include <vector>
  #include <cstdlib>
  #include <chrono>
  #include <sys/time.h>
  #include <random>
  #include <iostream>
  #include <fstream>
  #include <stdio.h>


  using namespace Rcpp;
  using namespace std;

  class reactions {
  '

  ############ INIT RANDOM GENERATOR ##############
  initRandomGenerator<-'
  void initRandomSeed(){
  struct timeval start;
  gettimeofday(&start,NULL);
  seed_ = start.tv_usec;
  randomGenerator_.seed(seed_);
  }
  '

  ############## COMPUTE TOTAL RATE ###########
  computeTotRate<-'
  double computeTotalRate(){
    double sum = 0;
    for(vector<double>::iterator it=vRates_.begin();it!=vRates_.end();it++){
      sum += *it;
    }
    return sum;
  }
  '

  ############## WHICH REACTION ##############
  whichReac<-'
  int whichReaction(const double& totalRate){
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double r = uniform(randomGenerator_);
    double sum = 0.0 ;
    int indx = 0;
    vector<double>::iterator it;
    for(it=vRates_.begin() ; it!=vRates_.end() && sum<(r*totalRate) ; ++it, indx++){
      sum += *it;
    }
    indx--;
    return indx;
  }
  '

  ############### PERFORM REACTION ###############
  performReac<-'
  void performReaction(int reac, int nbTime){
    string indiv="";
    for(unsigned i=0;i<vFrom[reac].size();i++){
      indiv = vFrom[reac][i];
      compartments[indiv]-=nbTime;
    }
    for(unsigned i=0;i<vTo[reac].size();i++){
      indiv = vTo[reac][i];
      compartments[indiv]+=nbTime;
    }
  }
  '

  ############## EXP JUMP ################
  expJump<-"
  double expJump(const double& rate){
    return -(1.0 / rate) * log(1.0 - ((rand() % 100) / 100.0));
  }
  "

  ################### UPDATE reactions OUTPUT #######################
  updateReacOutput<-"
  void updatereactionsOutput(int indxReac, int nTime){
    dfreactions.push_back(vreactions[indxReac]);
    dfNrep.push_back(nTime);
  }
  "

  ################# REXP METHOD ################
  rexp_gil<-'
  double rexp(const double& rate){
    double rand = 0.0;
    try{
      std::exponential_distribution<double> exponential(rate);
      rand = exponential(randomGenerator_);
    }
    catch (int e){
      Rcpp::stop("Error in Gillespie::rexp(rate)");
    }
    return rand;
  }
  '

  ################# RPOIS METHOD ##############
  rpois_gil<-'
  unsigned rpois(const double& rate){
    unsigned rand = 0;
    try{
      std::poisson_distribution<unsigned> poisson(rate);
      rand = poisson(randomGenerator_);
    }
    catch (int e){
      Rcpp::stop("Error in Gillespie::rpois(rate)");
    }
    return rand;
  }
  '

  ############## CHECK POPULATION SIZE ##############
  checkCompartmentSize <- '
  bool checkCompartmentSize(){
    // returns a boolean : True if no negative population size, False otherwise
    double value = 0;
    auto result = find_if(compartments.begin(), compartments.end(), [&value](pair<const string, double> &mo) {return mo.second < 0; });
    return (result == std::end(compartments));
  }
  '
  ################## IS FROM SIZE ENOUGH ##########
  isFromSizeEnough<-'
  bool isFromSizeEnough(int reac, int nbReac){
    vector<bool> isEnough;
    int nbFrom = 0;
    string indiv;
    for (unsigned i = 0; i < vFrom[reac].size(); i++){
      if(find(vTo[reac].begin(),vTo[reac].end(),vFrom[reac][i]) == vTo[reac].end()){
        indiv = vFrom[reac][i];
        nbFrom =compartments[indiv];
        isEnough.push_back(nbFrom>=nbReac);
      }
    }
    return (is_true(all(as<LogicalVector>( wrap(isEnough)) == TRUE))) ;
  }
  '

  ##### UPDATE DATA OUTPUT #####
  updateDataOutputAll <- '
  string trajectoryHeaderToString(){
    	stringstream ss;
      map<string,double>::iterator it;
    	ss << "Time";
    	for (it = compartments.begin(); it != compartments.end() ; ++it){
    			ss << "\\t";
    			ss << it->first;
    	}
    	ss << "\\tReaction\\tNrep";
    	return ss.str();
  }

  string compartmentStatesToString(){
  	stringstream ss;
  	map<string,double>::iterator it = compartments.begin();
  	ss << it->second;
  	++it;
  	for (; it != compartments.end() ; ++it){
  		ss << "\t" << it->second;
  	}
  	return ss.str();
  }

  string whichReactionsToString(int indxReac, int nTime){
    stringstream ss;
    ss << vreactions[indxReac];
    ss << "\t";
    ss << nTime;
    return ss.str();
  }

  void updateOutput(int indxReac, int nTime, ofstream& outTrajectory) {
    if (!isOutFile_) {
      updateDataFrameOutput();
      updatereactionsOutput(indxReac, nTime);
    }
    else{
      outTrajectory << compartmentStatesToString() << "\t";
      outTrajectory << whichReactionsToString(indxReac, nTime) << endl;
    }
  }

  '

  ########### ADD FUNCTIONS TO COMPUTE VALUES ##########
  RtoCfunctions = ""
  updateParam=""
  function_params=""
  if(!is.null(functions)){
    n=length(functions)
    function_params = sapply(strsplit(functions,"="),"[[",1)
    for (i in 1:n) {
      if(grepl('[=]',functions[i])){
        name=strsplit(x=functions[i],split='=')[[1]][1]
        fonction=strsplit(x=functions[i],split='=')[[1]][2]
      } else if(grepl('[<-]',functions[i])){
        name=strsplit(x=functions[i],split='<-')[[1]][1]
        fonction=strsplit(x=functions[i],split='<-')[[1]][2]
      }
      RtoCfunctions=paste(RtoCfunctions,"double get",name,"(){\n\treturn (",fonction,") ;\n}\n",sep="")

      updateParam=paste(updateParam,"\t",name," = get",name,"() ; \n",sep="")
    }
  }

  ############ reactions #############
  lnR<-length(reactions)
  lignes<-reactions

  rates<-unlist(stringr::str_extract_all(string = lignes, pattern = "(?<=\\[).*(?=\\])"))
  vReac<-paste("vreactions.assign(",lnR,',"");',sep="")
  updateRates<-""
  for(i in 1:lnR){
    rates[i]<-gsub(".*\\[(.*)\\].*", "\\1", lignes[i])
    rate <- rates[i] # copy since the rates vector is modified
    if(grepl(pattern = "\\^",x = rates[i])){ ## in case of square -> convert to Rcpp pow(m,n) function
      extracted <- stringr::str_extract_all(rates[i], "[^\\^]+")[[1]]
      pos <- stringr::str_locate_all(rates[i],"\\^")[[1]]
      npos <- nrow(pos)
      for(j in 1:npos){
        tmp <- extracted[j] # what's before the jth '^'
        tmp <- unlist(str_extract_all(tmp, stringr::boundary("word"))) # get only the words before the '^'
        tosquare <- tmp[length(tmp)] # get the last word, the one preceeding the '^'
        if(j == npos){
          tmp <- stringr::str_sub(rate, start = (pos[j,1]+1)) # what's after the last '^'
          tmp <- regmatches(tmp,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",tmp, perl=TRUE)) # get only numerical values from that
          byN <-  as.numeric(unlist(tmp[1]))[1] # get the first one (in case there are other numbers after)
        } else{
          tmp <- stringr::str_sub(string = rate, start = (pos[j,1]+1), end=(pos[(j+1),1]-1)) # what's after the jth '^'
          tmp <- regmatches(tmp,gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*",tmp, perl=TRUE)) # get only numerical values from that
          byN <-  as.numeric(unlist(tmp[1]))[1] # get the first one (in case there are other numbers after)
        }
        old <- paste0(tosquare,"\\^",byN)
        new <- paste0("pow(",tosquare,",",byN,")")
        tmp <- gsub(pattern = old, x = rates[i], replacement = new)
        rates[i] <- tmp
      }
    }
    vReac<-paste(vReac,paste("vreactions[",i-1,']="',lignes[i],'";',sep=""),sep="\n\t")
    #### UPDATE RATES FUNCTION ####
    updateRates<-paste(updateRates,paste("vRates_[",i-1,"]=",rates[i],";",sep=""),sep="\n\t")
  }
  ## Definir les individus et les parametres
  indivNames<-c()
  for (i in 1:lnR) {
    tmp<-unlist(stringr::str_extract_all(strsplit(reactions[i],"[[]")[[1]][1],stringr::boundary("word")))
    indivNames[length(indivNames)+1]<-tmp
    tmp<-unlist(stringr::str_extract_all(strsplit(reactions[i],"[]]")[[1]][2],stringr::boundary("word")))
    if(!NA %in% tmp){
      indivNames[length(indivNames)+1]<-tmp
    }
  }
  indivNames=indivNames[!grepl("^[[:digit:]]+$",indivNames)]
  indivNames=unique(indivNames)

  all=unique(unlist(stringr::str_extract_all(unlist(reactions),stringr::boundary("word"))))
  paramsNames=all[which(all %notin% indivNames)]
  paramNames=paramsNames[!grepl("^[[:digit:]]+$",paramsNames)]
  lnP<-length(paramNames)
  lnI<-length(indivNames)

  ### PRIVATE ###
  doubles<-paste("double time_","deltaT_", sep=" , ")
  output<-paste("vector<double> ", "dfTimes;",sep="")
  getParam<-""
  compartments<-""
  for(i in 1:lnP){
    doubles<-paste(doubles,paramNames[i],sep=' , ')
    ##### GET PARAM FOR PUBLIC CONSTRUCTOR #####
    if(paramNames[i] %notin% function_params){
      getParam<-paste(getParam,paste(paramNames[i],' = parameters_[timePart]["',paramNames[i],'"];',sep=""),sep="\n\t")
    }
  }
  updateCompartments<-""
  updateOutput<-""
  nameVec<-'namevec[0]="Time";'
  nameVec<-paste(nameVec,"\n\tnamevec[",1,']="Reaction";',"\n\tnamevec[",2,']="Nrep";',sep="")
  longList<-"long_list[0]=dfTimes;"
  longList<-paste(longList,"\n\tlong_list[",1,']=dfreactions;',"\n\tlong_list[",2,']=dfNrep;',sep="")
  dfClear<-paste('dfreactions.clear();','dfNrep.clear();','dfTimes.clear();',sep='\n')
  for(i in 1:lnI){
    doubles<-paste(doubles,indivNames[i],sep=' , ')
    output<-paste(output,paste("vector<double> df",indivNames[i],";",sep=""),sep="\n\t ")
    #### GET INDIV FOR PUBLIC CONSTRUCTOR #####
    compartments<-paste(compartments,paste('compartments["',indivNames[i],'"] = initialStates_["',indivNames[i],'"];',sep=""),sep="\n\t")
    ############################################
    updateCompartments<-paste(updateCompartments,paste(indivNames[i],' = compartments["',indivNames[i],'"];',sep=""),sep="\n\t")
    updateOutput<-paste(updateOutput,paste("df",indivNames[i],".push_back(",indivNames[i],");",sep=""),sep="\n\t")
    longList<-paste(longList,paste("long_list[",i+2,"]=df",indivNames[i],";",sep=""),sep="\n\t")
    dfClear<-paste(dfClear,paste('df',indivNames[i],'.clear();',sep=""),sep='\n')
    nameVec<-paste(nameVec,paste("namevec[",i+2,']="',indivNames[i],'";',sep=""),sep="\n\t")
  }
  ############# vTo & vFrom #############
  To<-""
  From<-""
  for(i in 1:lnR){
    first<-strsplit(x=lignes[[i]],split="\\[")[[1]]
    second<-strsplit(x=lignes[[i]],split="->")[[1]]
    froms<-stringr::str_extract_all(first[1], stringr::boundary("word"))[[1]]
    tos<-stringr::str_extract_all(second[2], stringr::boundary("word"))[[1]]
    for (j in 1:length(froms)) {
      if(froms[j]!="0"){
        From<-paste(From,paste("vFrom[",i-1,'].push_back("',froms[j],'");',sep=""),sep="\n\t")
      }
    }
    for(k in 1:length(tos)){
      if(tos[k]!="0"){
        To<-paste(To,paste("vTo[",i-1,'].push_back("',tos[k],'");',sep=""),sep="\n\t")
      }
    }
  }


  constructor<-constructor_gillespie
  initRandomGenerator<-initRandomGenerator

  ############# INIT STRUCTURES ###########
  initStructures<-paste('void initStructures(){',compartments,"vector<string> indiv(0);",paste("vRates_.assign(",lnR,",0);",sep=""),paste("vFrom.assign(",lnR,",indiv);",sep=""),paste("vTo.assign(",lnR,",indiv);",sep=""),From,To,vReac,"updateCompartments();",'dfreactions.assign(1,"init");','dfNrep.assign(1,1);',sep="\n\t")
  initStructures<-paste(initStructures,'}',sep='\n')
  ############## UPDATE RATES ###########
  updateRates=paste(updateParam,updateRates,sep="")
  updateRates<-paste("void updateRates(){",updateRates,"}",sep="\n")

  ############## UPDATE COMPARTMENTS ###########
  updateCompartments<-paste("void updateCompartments(){",updateCompartments,"}",sep="\n")

  ############### CONTINUE CONDITION ################
  continueCondition<-paste("bool continueCondition(unsigned timePart){",paste("return ( computeTotalRate()!=0 && ((time_ < vTimes_[timePart+1]) && !(abs(vTimes_[timePart+1] - time_)<(deltaT_/10)))) ;",sep=""),sep="\n\t")
  continueCondition<-paste(continueCondition,"}",sep="\n")

  computeTotRate<-computeTotRate
  whichReac<-whichReac
  performReac<-performReac
  expJump<-expJump
  runSwitch<-runSwitch
  runDirect<-runDirect_gillespie
  directAlgo<-directAlgo_gilespie
  gillespieSimulation<-gillespieSimulation_gil

  listL<-paste("if(!isOutFile_){\n\t","int numCol=",lnI+3,";",sep="")
  listL<-paste(listL,"Rcpp::List long_list(numCol);","Rcpp::CharacterVector namevec(numCol);",sep="\n\t\t")
  listL<-paste(listL,longList,nameVec,'long_list.attr("names") = namevec;','long_list.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, dfTimes.size());','long_list.attr("class") = "data.frame";','res["traj"] = long_list;',sep="\n\t\t")
  listL<-paste(listL,'}','auto elapsed = stop_ - start_;','if(verbose_){Rcout << "Operation took: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() << " microseconds." << endl;}','return res;',sep='\n\t')

  gillespieSimulation<-paste(gillespieSimulation,listL,'}',sep="\n\t")
  # gillespieSimulation<-paste(gillespieSimulation,paste('else{','Rcpp::stop("Failure. You should try other parameter values.");','}',sep="\n"),'}',sep='\n')

  gillespieSimulation<-paste(gillespieSimulation,paste('else{','remove((outFile_.c_str()));','return NULL;','}',sep="\n"),'}',sep='\n')

  ################### UPDATE DATAFRAME OUTPUT ####################
  updateOutput<-paste("void updateDataFrameOutput(){","dfTimes.push_back(time_);",updateOutput,sep="\n\t")
  updateOutput<-paste(updateOutput,"}",sep="\n")

  updateReacOutput<-updateReacOutput

  ################### RE INITIALIZE DATA FRAME ##############
  reInitialize<-paste(dfClear,'time_ = vTimes_[0];','dfreactions.assign(1,"init");','dfNrep.assign(1,1);',compartments,'updateCompartments();','initParams(0);',sep='\n')
  reInitialize<-paste('void reInitializeDataFrame(){',reInitialize,sep='\n\t')
  reInitialize<-paste(reInitialize,'}',sep='\n')

  rexp<-rexp_gil
  rpois<-rpois_gil
  isFromSizeEnough<-isFromSizeEnough
  checkCompartmentSize <- checkCompartmentSize
  runTau<-runTau
  tauLeapAlgo<-tauLeapAlgo
  updateDataOutputAll <- updateDataOutputAll

  destructor<-'~reactions(){}'

  ########## INIT PARAMS #########
  initParams<-paste('void initParams(unsigned timePart){',getParam,'updateRates();',sep="\n\t")
  initParams<-paste(initParams,"}",sep="\n")


  ################### PUBLIC PRIVATE ###################
  ######################################################
  doubles<-paste(doubles,"", ";")
  private<-paste("private:","vector<double> vRates_;","long seed_;","int nTrials_;",doubles,"NumericVector initialStates_;","NumericVector vTimes_;","map<string,double> compartments;","vector<vector<string> > vFrom;","vector<vector<string> > vTo;","vector<string> vreactions;",output,"vector<string> dfreactions;","vector<int> dfNrep;","std::mt19937 randomGenerator_;","vector<NumericVector> parameters_;","bool isSwitchingMode_;", "bool verbose_;", "double msaTau_;", "int msaIt_;","string outFile_;","bool isOutFile_;",sep='\n\t')

  public<-paste(constructor,destructor,initRandomGenerator,initStructures,initParams,reInitialize,updateRates,updateCompartments,computeTotRate,whichReac,isFromSizeEnough,checkCompartmentSize,performReac,rexp,rpois,gillespieSimulation,runSwitch,runDirect,directAlgo,runTau,tauLeapAlgo,continueCondition,updateOutput,updateReacOutput,updateDataOutputAll,RtoCfunctions,"};",sep="\n")
  public<-paste("public:",public,sep="\n\t")

  rcppModule<-rcppModule_gillespie
  header<-header_gillespie
  ############### CODE ##############
  # codeSimulation<-paste(header,private,public,rcppModule,sep="\n")
  codeSimulation<-paste(header,private,public,sep="\n")

}

  ###########################################
  ######### BUILD COMPILE SIMULATOR #########
  ###########################################

  #' Build a simulator of dynamics of population-model
  #'
  #' A simulator is built by supplying reactions of the model described by our formalism or described by differential equations
  #' The returned function will be used to simulate trajectories, that can further be used to simulate phylogenies.
  #'
  #' @param reactions A character vector of reactions describing the input model.
  #' @param functions A named vector where functions are defined.
  #'
  #' @return An object of class \code{simulation}, which is a function that can be used to simulate trajectories from the model.
  #' @author Gonche Danesh
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' # Build a simulator for an SIR model
  #' reactions <- c('S [beta * S * I] -> I',
  #'                'I [gamma * I] -> R')
  #'
  #' sir.simu <- build_simulator(reactions = reactions)
  #'
  #' # Run a simulation of a trajectory
  #' sir_traj <- sir.simu(paramValues = c(gamma = 1, beta = 2e-4),
  #'                      initialStates = c(I = 1, S = 9999, R = 0),
  #'                      times = c(0, 20), ,
  #'                      method = "exact",
  #'                      seed = 166)
  #'
  #' # The output is a named list containing the trajectory, the algorithm used,
  #' # the parameter values and the reactions of the model.
  #' names(sir_traj)
  #'
  #' # Print head of the simulated trajectory
  #' head(sir_traj$traj)
  #'
  #' # Plot the trajectory
  #' plot(sir_traj)
  #' }
  build_simulator<-function(reactions, functions=NULL){

    ### Check for redonduncy
    # tmp = as.vector(unlist(compartmentNames))
    # if(length(str_replace(tmp,"[*]","")) != length(compartmentNames)){
    #   stop('List of compartment names has redundancy. Names must be different.')
    # } else if(length(unique(paramNames)) != length(paramNames)){
    #   stop('List of parameter names has redundancy. Names must be different.')
    # }

    src='
    // [[Rcpp::plugins(cpp11)]]
    List theta_(paramValues);
    NumericVector times_(times);
    NumericVector init_(initialStates);
    double dT_ = as<double>(dT);
    int nTrial_ = as<int>(nTrials);
    long seed_ = as<long>(seed);

    bool switching_ = LOGICAL(x)[0];
    bool verbose_ = LOGICAL(x)[1];

    double msaTau_ = as<double>(msaTau);
    int msaIt_ = as<int>(msaIt);

    string outFile_ = as<string>(outFile);

    reactions simu(theta_,init_,times_,dT_,nTrial_,switching_,verbose_,seed_, msaTau_, msaIt_, outFile_);
    List traj = simu.GillespieSimulation();
    return traj;
    '
    codeSimu <- simulator_generator(reactions, functions)
    message("Compiling the simulator ...")

    # ptm <- proc.time()

    ## For c++ compilation
    myplugin <- getPlugin("Rcpp")
    myplugin$env$PKG_CXXFLAGS <- "-std=c++11"

    simu.cpp<-cxxfunction(
      signature(paramValues="List",times="NumericVector",initialStates="NumericVector",
                dT="numeric",nTrials="integer",x="logical",seed="integer", msaTau="numeric", msaIt="numeric", outFile="character"
      ),
      plugin="Rcpp",body=src,includes=codeSimu, settings=myplugin
    )

    lnR = length(reactions)
    simulation<-function(paramValues,initialStates,times,method="approximate",tau=0.01,nTrials=1, verbose=FALSE, seed=0, msaTau=tau/10, msaIt=10, outFile=""){
      ok = TRUE
      traj = list()

      if(method %notin% c("exact","approximate","mixed")){
        ok=FALSE
        message('The algorithm given in the method argument given is not correct.\nThe different methods are "exact", "approximate" or "mixed."')
      }

      # if( length(times) != (length(paramValues)+1) ){
      #   print("Error : number of vectors of parameter values doesn't correspond to the number of time intervals.")
      #   ok = FALSE
      # } else
      if(! .checkParamValues(paramValues,times,initialStates)){
        warning("Error : parameter or individual or times values contain NA value.")
        ok = FALSE
      } else{
        ok = TRUE

        #### Transformer theta sous forme de liste de vecteurs
        ## Combien d'intervalles de temps
        nb<-length(paramValues[[1]])
        for (i in 2:length(paramValues)) {
          nb<-max(nb,length(paramValues[[i]]))
        }
        parameters<-names(paramValues)
        theta<-list()
        for (j in 1:nb) {
          tmp<-c()
          for (i in 1:length(parameters)) {
            len<-length(paramValues[[i]])
            if(j > len){
              tmp[length(tmp)+1]=paramValues[[i]][len]
            } else{
              tmp[length(tmp)+1]=paramValues[[i]][j]
            }
          }
          names(tmp)<-parameters
          theta[[j]]<-tmp
        }


      }

      if(method == "exact"){
        switchingMode = FALSE
        algo = "Exact Method"
        dT = 0
      } else if(method == "approximate"){
        switchingMode = FALSE
        algo = "Tau-Leap Method"
        dT = tau
      } else if(method == "mixed"){
        switchingMode = TRUE
        algo = "Mixed Method"
        dT = tau
        ## ajouter les seuils
      }

      results<-list()

      if(ok){
        x=c(switchingMode, verbose)
        traj = simu.cpp(paramValues=theta,initialStates=initialStates,times=times,dT=dT,nTrials=nTrials,x=x,seed=seed,msaTau=msaTau,msaIt=msaIt,outFile=outFile)
        if(length(traj) > 0){

          message("Success!")

          reacs <- unique(traj$Reaction)

          results[["reactions"]] <- reactions
          results[["values"]] <- c(paramValues,initialStates)
          results[["times"]] <- times
          results[["method"]] <- algo
          results[["tau"]] <- dT
          results[["seed"]] <- traj$seed
          if(outFile == ""){
            results[["traj"]] <- traj$traj
          }else{
            results[["outFile"]] <- outFile
          }
          class(results) <- "simutraj"
        }else{
          message("Failure. You should try other parameter values.")
        }

      }

      results
      # traj

    }

    message("Process complete.")

    class(simulation) <- c('simulation', 'function')
    simulation

  }

  ## Plot a trajectory ; input is the resulted output from the simulation

  #' Plot an object of class \code{simutraj}.
  #'
  #' @param x An object of \code{simutraj} resulting from running a simulator of trajectories built using the \code{build_simulator} function.
  #' @param ... Arguments to be passed to methods, such as graphical parameters
  #' @author Gonche Danesh
  #' @rdname plot.simutraj
  #' @export
  plot.simutraj <- function(x, ...){
    indx_indiv <- which(names(x$traj) %notin% c("Time","Reaction","Nrep"))
    len_indx <- length(indx_indiv)
    ymin <- min(unlist(lapply(indx_indiv,function(i) min(x$traj[[i]]))))
    ymax <- max(unlist(lapply(indx_indiv,function(i) max(x$traj[[i]]))))
    ylim <- c(ymin, ymax)
    xlim <- c(min(x$times),max(x$times))
    cl <- grDevices::rainbow(len_indx)
    graphics::plot(0,0,xlim=xlim,ylim=ylim,type = "n",xlab="Time",ylab="Number of individuals")
    invisible(lapply(1:len_indx, function(i) graphics::lines(y=x$traj[[i+3]], x=x$traj$Time, col=cl[i])))
    my.legend.size <- graphics::legend("topright",names(x$traj)[indx_indiv],col=cl[1:len_indx],lty=1)
  }
