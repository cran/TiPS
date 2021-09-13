
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
#' @param reactions a vector of reactions
#' @param functions a list of functions
#' @NoRd
simulator_generator<-function(reactions, functions=NULL){
  
  message("Building the simulator ...")


    ############### CONSTRUCTEUR ###########
  constructor_gillespie<-'
  reactions(List params, NumericVector initialStates, NumericVector times, double deltaT, int nTrials, bool isSwitchingMode, bool verbose) : nTrials_(nTrials), deltaT_(deltaT), isSwitchingMode_(isSwitchingMode){
    verbose_ = verbose;
  for (unsigned i = 0; i < params.size(); i++){
  parameters_.push_back(params[i]);
  }
  initialStates_ = initialStates;
  vTimes_ = times;
  initStructures();
  initRandomSeed();
  randomGenerator_.seed(seed_);
  time_ = vTimes_[0];
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
  for ( ; timePart<timePartLim && computeTotalRate()!=0 && nEvents>=0 ; timePart++) {
  initParams(timePart);
  totalRate = computeTotalRate();
  jump = rexp(totalRate);
  while(continueCondition(timePart) && nEvents>=0 ){
  if(jump < deltaT_/10.0){
  switchingThreshold ++;
  }
  else{
  switchingThreshold = 0;
  }
  if(switchingThreshold > 10 || ((double)nEvents/(double)vreactions.size()) > 1){
  nEvents = tauLeapAlgo();
  jump = deltaT_;
  switchingThreshold = 0;
  }
  else{
  jump = directAlgo();
  nEvents =1;
  }
  }
  }
  return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
  }
  '

  ######### RUN DIRECT ##########
  runDirect_gillespie<-'
  bool runDirect(){
  unsigned timePart = 0;
  time_ = vTimes_[timePart];
  int timePartLim = (vTimes_.size()-1);
  initParams(timePart);
  for ( ; timePart<timePartLim && computeTotalRate()!=0 ; timePart++) {
  initParams(timePart);
  while(continueCondition(timePart)){
  double totalRate = 0.0;
  double jump = 0.0;
  int reactIndx = 0;
  totalRate = computeTotalRate();
  jump = rexp(totalRate);
  time_ += jump;
  reactIndx = whichReaction(totalRate);
  performReaction(reactIndx,1);
  updateCompartments();
  updateDataFrameOutput();
  updatereactionsOutput(reactIndx,1);
  updateRates();
  }
  }
  return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
  }
  '

  ############# DIRECT ALGORITHM ###########
  directAlgo_gilespie<-'
  double directAlgo(){
  double totalRate = 0.0;
  double jump = 0.0;
  int reactIndx = 0;
  totalRate = computeTotalRate();
  jump = rexp(totalRate);
  time_ += jump;
  reactIndx = whichReaction(totalRate);
  performReaction(reactIndx,1);
  updateCompartments();
  updateDataFrameOutput(); 
  updatereactionsOutput(reactIndx,1);
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
  updateDataFrameOutput();
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
  updateDataFrameOutput();
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
  updateDataFrameOutput();
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
  '

  ############# RUN TAU #############
  runTau<-'
  bool runTau(){
    unsigned timePart = 0;
    time_ = vTimes_[timePart];
    int timePartLim = (vTimes_.size()-1);
    initParams(timePart);
    int nEvents = 0;
    for ( ; timePart<timePartLim && computeTotalRate()!=0 && nEvents>=0 ; timePart++) {
      initParams(timePart);
      while(continueCondition(timePart) && nEvents>=0 ){
        nEvents = tauLeapAlgo();
      }
    }
    return (time_>vTimes_[vTimes_.size()-1] || (abs(vTimes_[vTimes_.size()-1] - time_)<(deltaT_/10)));
  }
  '

  ############ TAU LEAP ALGO #########
  tauLeapAlgo<-'
  int tauLeapAlgo(){
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
        int reac = reactionsToPerform[i].first;
        int nbTime = reactionsToPerform[i].second;
        performReaction(reac, nbTime);
        updateCompartments();
        updateDataFrameOutput();
        updatereactionsOutput(reac,nbTime);
        updateRates();
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
  .constructor<List, NumericVector, NumericVector, double, long, int, bool>("documentation for constructor")
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

  using namespace Rcpp;
  using namespace std;

  class reactions {
  '

  ########## RANDOM GENERATOR ############
  randGenerator<-'
  void initRandomSeed(){
  struct timeval start;
  gettimeofday(&start,NULL);
  seed_ = start.tv_usec;
  randomGenerator_.seed(seed_); 
  }
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
  whichReac<-"
  int whichReaction(const double& totalRate){
    double r = (double)rand() / RAND_MAX ;
    double sum = 0.0 ;
    int indx = 0;
    vector<double>::iterator it;
    for(it=vRates_.begin() ; it!=vRates_.end() && sum<(r*totalRate) ; ++it, indx++){
      sum += *it;
    }
    indx--;
    return indx;
  }
  "

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
  
  ############ reactions #############
  lnR<-length(reactions)
  lignes<-reactions
  
  rates<-unlist(stringr::str_extract_all(string = lignes, pattern = "(?<=\\[).*(?=\\])"))
  vReac<-paste("vreactions.assign(",lnR,',"");',sep="")
  updateRates<-""
  for(i in 1:lnR){
    rates[i]<-gsub(".*\\[(.*)\\].*", "\\1", lignes[i])
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
  
  all=unique(unlist(stringr::str_extract_all(unlist(rates),stringr::boundary("word"))))
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
    getParam<-paste(getParam,paste(paramNames[i],' = parameters_[timePart]["',paramNames[i],'"];',sep=""),sep="\n\t")
    ############################################
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
  
  
  randGenerator<-randGenerator
  constructor<-constructor_gillespie
  initRandomGenerator<-initRandomGenerator
  
  ########### ADD FUNCTIONS TO COMPUTE VALUES ##########
  RtoCfunctions = ""
  updateParam=""
  if(!is.null(functions)){
    n=length(functions)
    
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
  
  listL<-paste("int numCol=",lnI+3,";",sep="")
  listL<-paste(listL,"Rcpp::List long_list(numCol);","Rcpp::CharacterVector namevec(numCol);",sep="\n\t")
  listL<-paste(listL,longList,nameVec,'long_list.attr("names") = namevec;','long_list.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, dfTimes.size());','long_list.attr("class") = "data.frame";','auto elapsed = stop_ - start_;','     if(verbose_){Rcout << "Operation took: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() << " microseconds." << endl;}','return long_list;',sep="\n\t")
  
  gillespieSimulation<-paste(gillespieSimulation,listL,'}',sep="\n\t")
  # gillespieSimulation<-paste(gillespieSimulation,paste('else{','Rcpp::stop("Failure. You should try other parameter values.");','}',sep="\n"),'}',sep='\n')
  
  gillespieSimulation<-paste(gillespieSimulation,paste('else{','return NULL;','}',sep="\n"),'}',sep='\n')
  
  ################### UPDATE DATAFRAME OUTPUT ####################
  updateOutput<-paste("void updateDataFrameOutput(){","dfTimes.push_back(time_);",updateOutput,sep="\n\t")
  updateOutput<-paste(updateOutput,"}",sep="\n")
  
  updateReacOutput<-updateReacOutput
  
  ################### RE INITIALIZE DATA FRAME ##############
  reInitialize<-paste(dfClear,'time_ = vTimes_[0];','dfreactions.assign(1,"init");','dfNrep.assign(1,1);',compartments,'updateCompartments();','initParams(0);','updateDataFrameOutput();',sep='\n')
  reInitialize<-paste('void reInitializeDataFrame(){',reInitialize,sep='\n\t')
  reInitialize<-paste(reInitialize,'}',sep='\n')
  
  rexp<-rexp_gil
  rpois<-rpois_gil
  isFromSizeEnough<-isFromSizeEnough
  runTau<-runTau
  tauLeapAlgo<-tauLeapAlgo
  
  destructor<-'~reactions(){}'
  
  ########## INIT PARAMS #########
  initParams<-paste('void initParams(unsigned timePart){',getParam,'updateRates();',sep="\n\t")
  initParams<-paste(initParams,"}",sep="\n")
  
  
  ################### PUBLIC PRIVATE ###################
  ######################################################
  doubles<-paste(doubles,"", ";")
  private<-paste("private:","vector<double> vRates_;","long seed_;","int nTrials_;",doubles,"NumericVector initialStates_;","NumericVector vTimes_;","map<string,double> compartments;","vector<vector<string> > vFrom;","vector<vector<string> > vTo;","vector<string> vreactions;",output,"vector<string> dfreactions;","vector<int> dfNrep;","std::mt19937 randomGenerator_;","vector<NumericVector> parameters_;","bool isSwitchingMode_;", "bool verbose_;",sep='\n\t')
  
  public<-paste(constructor,destructor,initRandomGenerator,initStructures,initParams,reInitialize,updateRates,updateCompartments,computeTotRate,whichReac,isFromSizeEnough,performReac,rexp,rpois,gillespieSimulation,runSwitch,runDirect,directAlgo,runTau,tauLeapAlgo,continueCondition,updateOutput,updateReacOutput,RtoCfunctions,"};",sep="\n")
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
  #' # Build a simulator for an SIR model
  #' reactions <- c('I [beta * S * I] -> I',
  #'                'I [gamma * I] -> R')
  #'                
  #' sir.simu <- build_simulator(reactions = reactions)
  #' 
  #' # Run a simulation of a trajectory 
  #' sir_traj <- sir.simu(paramValues = c(gamma = 1, beta = 2e-4),
  #'                      initialStates = c(I = 1, S = 9999, R = 0),
  #'                      times = c(0, 20))
  #'                  
  #' # The output is a named list containing the trajectory, the algorithm, the parameter values and the reactions of the model.
  #' names(sir_traj)
  #' 
  #' # Print head of the simulated trajectory
  #' head(sir_traj$traj) 
  #' 
  #' # Plot the trajectory
  #' plot(sir_traj)
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
    
    bool switching_ = LOGICAL(x)[0];
    bool verbose_ = LOGICAL(x)[1];
    
    reactions simu(theta_,init_,times_,dT_,nTrial_,switching_,verbose_);
    
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
                dT="numeric",nTrials="integer",x="logical"
      ),
      plugin="Rcpp",body=src,includes=codeSimu, settings=myplugin
    )
    
    simulation<-function(paramValues,initialStates,times,method="mixed",tau=0.01,nTrials=1, verbose=FALSE){
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
        traj = simu.cpp(paramValues=theta,initialStates=initialStates,times=times,dT=dT,nTrials=nTrials,x=x)
        
        if(length(traj) > 0){

          message("Success!")

          reacs <- unique(traj$Reaction)
          
          results[["reactions"]] <- reactions
          results[["values"]] <- c(paramValues,initialStates)
          results[["times"]] <- times
          results[["method"]] <- algo
          results[["tau"]] <- dT
          results[["traj"]] <- traj
          
          class(results) <- "simutraj"
        } else{
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
  #' @param simuResults An object of \code{simutraj} resulting from running a simulator of trajectories built using the \code{build_simulator} function.
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
