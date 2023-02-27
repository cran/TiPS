#include "Phyloepid.h"
#include <typeinfo>
#include <chrono>
#include <sys/time.h>
#include <random>

// static boost::posix_time::ptime start_, stop_;

using namespace std;
using namespace Rcpp;


Phyloepid::Phyloepid(List reactions, List traj, bool fullTraj, bool isReSampling,  unsigned int nbdates, bool verbose, NumericVector options) :
		compartments_(),
		reactions_(),
		roots_(),
		leafcount_(0),
		fullTraj_(fullTraj),
		traj_(traj),
		outTree_(""),
		nTrials_(options["nTrials"]),
		strReactions_(),
		reSampling_(isReSampling),
		nbdates_(nbdates),
		verbose_(verbose),
		seed_(options["seed"])
		//reacFind_(false)
{
	initCompartments();
	readReactions(reactions);
	if (seed_ == 0){
		initRandomSeed();
	}
	randomGenerator_.seed(seed_);
}

// void Phyloepid::initRandomSeed(){
// 	boost::posix_time::ptime epoch = boost::posix_time::second_clock::local_time();
// 	boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
// 	boost::posix_time::time_duration tzero = now - epoch;
// 	seed_ = tzero.total_microseconds();
// }

void Phyloepid::initRandomSeed(){
	struct timeval start;
	gettimeofday(&start,NULL);
	seed_ = start.tv_usec;
	randomGenerator_.seed(seed_);
}


bool Phyloepid::simulationTree(){
	if (verbose_){
		Rcout << "Running simulation of the tree based on the trajectory..." << endl;
	}
	bool ok = false;
	int i = 0;

	auto start_ = std::chrono::high_resolution_clock::now();
	auto stop_ = std::chrono::high_resolution_clock::now();

	ok = run();
	if (ok){
		roots_[0]->clean();
		while(roots_[0]->getNbSons() == 1 && !(roots_[0]->getSons()[0]->isLeaf())){
			Node* tmpSon = roots_[0]->getSons()[0];
			roots_[0]->removeSon(tmpSon);
			roots_[0] = tmpSon;
		}
		roots_[0]->initializeDistances();
		stop_ = std::chrono::high_resolution_clock::now();
		if (roots_[0]->getNbLeaves() != nbdates_){
			ok = false;
		}
		// else{
		// 	// Rcout << "Success !" << endl;
		// 	auto elapsed = stop_ - start_;
		// 	// Rcout << "Operation took: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() << " microseconds." << endl;
		// }
	}
	stop_ = std::chrono::high_resolution_clock::now();
	for (i=1 ; i < nTrials_ && !ok ; i++){
		if (verbose_){
			Rcout << "- Trial " << (i+1) << "..." << endl;
		}
		start_ = std::chrono::high_resolution_clock::now();
		ok = run();
		if (ok){
			roots_[0]->clean();
			while(roots_[0]->getNbSons() == 1 && !(roots_[0]->getSons()[0]->isLeaf())){
				Node* tmpSon = roots_[0]->getSons()[0];
				roots_[0]->removeSon(tmpSon);
				roots_[0] = tmpSon;
			}
			roots_[0]->initializeDistances();
			stop_ = std::chrono::high_resolution_clock::now();
			if (roots_[0]->getNbLeaves() != nbdates_){
				ok = false;
			}
			// else{
			// 	Rcout << "Success !" << endl;
			// 	auto elapsed = stop_ - start_;
			// 	Rcout << "Operation took: " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() << " microseconds." << endl;
			// }
		}
	}

	// if(ok){
	// 	// Rcout << "Success !" << endl;
	// 	// Rcout << "#Random seed = " << endl << seed_ << endl;
	// 	//Rcout << "#Number of trees: " << endl;
	// 	//Rcout << roots_.size() << endl;
	// 	//Rcout << "#Number of leaves:" << endl;
	// 	// for (unsigned i=0 ; i<roots_.size() ; i++){
	// 		roots_[0]->clean();
	// 		while(roots_[0]->getNbSons() == 1 && !(roots_[0]->getSons()[0]->isLeaf())){
	// 			Node* tmpSon = roots_[0]->getSons()[0];
	// 			roots_[0]->removeSon(tmpSon);
	// 			roots_[0] = tmpSon;
	// 		}
	// 		roots_[0]->initializeDistances();
	// 		if (roots_[0]->getNbLeaves() != nbdates_){
	// 			ok = false;
	// 			Rcout << "Failure." << endl;
	// 		} else{
	// 			boost::posix_time::time_duration duration = stop_ - start_;
	// 			Rcout << "Operation took: " << duration.total_microseconds() << " microseconds." << endl;
	// 		}
	// 		// Rcout << "##Tree " << i+1 << ": " << roots_[i]->getNbLeaves() << endl;
	// 	// }
	// 	// boost::posix_time::time_duration duration = stop_ - start_;
	// 	// Rcout << "Operation took: " << duration.total_microseconds() << " microseconds." << endl;
	// }
	// else{
	// 	Rcout << "Failure." << endl;
	// }
	// if(!ok){
	// 	Rcout << "Failure." << endl;
	// }
	return ok;
}

string Phyloepid::getNexusTree(const bool& withInfos){
	stringstream tree;
	tree << "#NEXUS" << endl;
	tree << "begin taxa;" << endl;;
	tree << "\t" << "dimensions ntax=" << (leafcount_-1) << ";" << endl;
	tree << "\t" << "taxlabels" << endl;
	for(int i=1 ; i<leafcount_ ; i++){
		tree << "\t\"I_" << i << "\"" << endl ;
	}
	tree << ";" << endl;
	tree << "end;" << endl << endl;
	tree << "begin trees;" << endl;
	for (unsigned i=0 ; i<roots_.size() ; i++){
		tree << "\t" << "tree TREE" << (i+1) << " = [&R] " << roots_[i]->getNewick(withInfos) << endl;
	}
	tree << "end;" << endl;
	return tree.str();
}

string Phyloepid::getNewickTree(const bool& withInfos){
	stringstream tree;
	for (unsigned i=0 ; i<roots_.size() ; i++){
		tree << roots_[i]->getNewick(withInfos) << endl;
	}
	return tree.str();
}

void Phyloepid::readReactions(List reactions) {
	Reaction* tmpReaction;
	string tmpTo;
	string tmpFrom;
	for (unsigned i = 0; i < reactions.size(); ++i){
		// Rcpp::print(reactions[i]);

		strReactions_.push_back(reactions[i]);

		tmpReaction = new Reaction();
	  	string strReaction = reactions[i];
	  	size_t equalPos = strReaction.find('=');
	  	if (equalPos == string::npos){
	  		/* Si la réaction est 'I' alors il s'agit d'un échantillonnage de I */
	  		tmpReaction->addTo(compartments_[strReaction]);
	  		tmpReaction->setIsSampling();
	  	}
	  	else{
	  		char sign = strReaction.at(equalPos-1);
	  		size_t endTo = equalPos-1;
	  		switch(sign){
	  			case '+':
	  				tmpReaction->setIsBirth();
	  				break;
	  			case '-':
	  				tmpReaction->setIsDeath(fullTraj_);
	  				break;
	  			default:
	  				sign = ' ';
	  				endTo++;
	  				tmpReaction->setIsMigration();
	  				break;
	  		}
	  		tmpTo = strReaction.substr(0,endTo);
	  		if(tmpTo.size() != 0){
	  			parseTo(tmpTo, tmpReaction);
	  		}
	  		if(sign != '-'){
	  			size_t endRate = strReaction.find(']');
	  			if(endRate == string::npos){
	  				endRate = equalPos ;
	  			}
	  			tmpFrom = strReaction.substr(endRate+1);
	  			parseFrom(tmpFrom, tmpReaction);
	  		}
	  	}
	  	reactions_[strReaction] = tmpReaction ;
	}
}


void Phyloepid::parseTo(const string& to, Reaction* reaction){
	string tmp = to;
	size_t plusPos = tmp.find('+');
	while(plusPos != string::npos){
		reaction->addTo(compartments_[tmp.substr(0,plusPos)]);
		tmp = tmp.substr(plusPos+1);
		plusPos = tmp.find('+');
	}
	reaction->addTo(compartments_[tmp]);
}

void Phyloepid::parseFrom(const string& from, Reaction* reaction){
	string tmp = from;
	size_t plusPos = tmp.find('+');
	while(plusPos != string::npos){
		reaction->addFrom(compartments_[tmp.substr(0,plusPos)]);
		tmp = tmp.substr(plusPos+1);
		plusPos = tmp.find('=');
	}
	reaction->addFrom(compartments_[tmp]);
}

void Phyloepid::initCompartments(){
	Compartment* tmpComp;
	vector<string> colNames = traj_.names();
	for(unsigned i = 3 ; i < colNames.size() ; i++){
		tmpComp = new Compartment();
		tmpComp->setName(colNames[i]);
		compartments_[colNames[i]] = tmpComp ;
		compartmentNames_.push_back(colNames[i]);
		vector<long> compTraj = traj_[colNames[i]];
		compTrajectories_[colNames[i]] = compTraj;

	}
}

void Phyloepid::updateCompartments(){
	for (map<string,Compartment*>::iterator it = compartments_.begin(); it!=compartments_.end() ; ++it){
		it->second->update();
	}
}

bool Phyloepid::run(){
	roots_.clear();
	bool ok = false;
	bool continue_ = true;
	double tmpTime;
	double tmpOldTime = -1.0;
	Reaction* tmpReaction;
	string strReaction;
	long nrep = 0;
	leafcount_ = 1;
	unsigned unrootedCount = 0;
	unsigned k = 0;

	//bool reacFind = false;

	vector<string> colNames = traj_.names();
	vector<double> times = traj_["Time"];
	vector<string> reactions = traj_["Reaction"];
	vector<long> nreps = traj_["Nrep"];

	unsigned tmpOldIndxTraj = 0;

	for(unsigned i = 0 ; i < nreps.size() && !ok && continue_ && leafcount_ >= 0 ; i++){
		tmpTime = times[i];
		strReaction = reactions[i];


		if (find(strReactions_.begin(), strReactions_.end(), strReaction) == strReactions_.end()){
		    warning(" Error : Reactions given in input must be similar than those in the trajectory. \nReaction ", strReaction, " at time ", tmpTime, " is not included the input reactions.", "\nPlease enter correct reactions. ");
			return false;
		}
		else{ //reaction found
			tmpReaction = reactions_[strReaction];
			nrep = nreps[i];

			if(tmpTime != tmpOldTime){
				k = 0;
				tmpOldIndxTraj = i;
				if(tmpReaction->isSampling()){
					k++;
					while(reactions_[reactions[i+k]]->isSampling()){
						k++;
					}
					tmpOldIndxTraj = i+k;
				}
				continue_ = updateDemeCompartments(tmpOldIndxTraj);
			}

			if(continue_){
				tmpOldTime = tmpTime;
				leafcount_ = tmpReaction->perform(nrep,strReaction,tmpOldTime,compTrajectories_,tmpOldIndxTraj, leafcount_, &roots_, reSampling_, fullTraj_);
				if (leafcount_ == -2){
					ok = false;
					leafcount_= 0;
				} else if( (leafcount_ != -1) & (leafcount_ != -2) ){
					unrootedCount = sumUnrootedNodes(); // compte le nombre de sous arbres non enracinés
					ok = (unrootedCount == 0); // determine si la simulation d'arbre est terminee en fonction de s'il reste encore des sous arbres non-enracinés
					if (verbose_ & !ok) {
						warning("%i unrooted nodes left.", unrootedCount);
					}
				} else{
					return false;
				}
			}
		}
	}
	return ok;
}

bool Phyloepid::updateDemeCompartments(unsigned indxTraj){
	bool ok = true ;
	bool update = true;
	long tmpSize;
	for(unsigned i = 0 ; i < compartmentNames_.size() ; i++){
		tmpSize =  compTrajectories_[compartmentNames_[i]][indxTraj] ;
		compartments_[compartmentNames_[i]]->setSize(tmpSize) ;
		update = compartments_[compartmentNames_[i]]->updateNodes() ;
		ok = (ok && update) ;
	}
	return ok;
}

long Phyloepid::sumUnrootedNodes(){
		long sum = 0;
		for (map<string,Compartment*>::iterator it = compartments_.begin() ; it != compartments_.end() ; ++it){
			sum += it->second->getNodeSize();
		}
		return sum;
}

Phyloepid::~Phyloepid() {}

vector<double> Phyloepid::getEdgeLengths() const{
	return roots_[0]->getBranchLengths() ;
}

vector<int> Phyloepid::getNbNodes() const{
	// first element id nmber of tips, second element is number of inner nodes
	pair<int,int> tmp;
	tmp = roots_[0]->getNbNodes();
	vector<int> nbnodes(2);
	nbnodes[0] = tmp.first;
	nbnodes[1] = tmp.second;
	return nbnodes;
}

vector<string> Phyloepid::getTipLabels() const{
	// vector<string> tiplabels ;
	// tiplabels = roots_[0]->getTipLabels();
	return roots_[0]->getTipLabels();
}

List Phyloepid::createTreeObject() const{
	vector<string> tiplabels = roots_[0]->getTipLabels();
	vector<double> tipHeigths = roots_[0]->getTipHeights();
	vector<string> nodelabels = roots_[0]->getNodeLabels();

	// int nbtips = tiplabels.size() ;
	int nbtips = roots_[0]->setLeavesID(0) ;
	int nbInodes = roots_[0]->setInnerNodesID(nbtips,0);
	// int nbEdges = nbtips + nbInodes -1 ;
	vector<double> edgeLengths = roots_[0]->getBranchLengths();
	map<string,vector<int>> edges = roots_[0]->getEdges();

	// IntegerMatrix mat(2,nbEdges);
	// mat.row(0) = edges["from"];
	// mat.row(1) = edges["to"];
	// mat(_,0) = edges["from"];
	// mat(_,1) = edges["to"];

	List phylogeny = List::create(Named("edge.length") = edgeLengths , _["tip.label"] = tiplabels, _["from"] = edges["from"], _["to"] = edges["to"], _["Nnode"] = nbInodes, _["node.label"] = nodelabels, _["tip.height"] = tipHeigths, _["seed"] = seed_);

	return phylogeny;
}

// static void phyloepid_finalizer(Phyloepid* ptr){
//     if (ptr){
//         delete ptr;
//     }
// }
RCPP_EXPOSED_CLASS(phyloepid)
RCPP_MODULE(phyloepid){
    Rcpp::class_<Phyloepid>( "Phyloepid" )
        .constructor<List,List,bool,bool,unsigned int,bool,NumericVector>("documentation for constructor")
	.method( "readReactions", &Phyloepid::readReactions, "reading model reactions")
	.method( "simulationTree", &Phyloepid::simulationTree, "simulation of the tree")
	.method( "getNexusTree", &Phyloepid::getNexusTree, "get simulated tree in Nexus format")
	.method( "getNewickTree", &Phyloepid::getNewickTree, "get simulated tree in Newick format" )
	.method( "getEdgeLengths", &Phyloepid::getEdgeLengths, "get branch lengths" )
	.method( "getNbNodes", &Phyloepid::getNbNodes, "get number of tips and number of inner nodes" )
	.method( "getTipLabels", &Phyloepid::getTipLabels, "get tip labels" )
	.method( "createTreeObject", &Phyloepid::createTreeObject, "create R tree object" )
	// .finalizer(&phyloepid_finalizer)
    ;
}
