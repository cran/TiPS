#include "Reaction.h"
using namespace std;

Reaction::Reaction() : type_(none), from_(), to_(), count_(0)
{}

Reaction::~Reaction()
{}

void Reaction::addFrom(Compartment* from)
{
	from_.push_back(from);
}

void Reaction::addTo(Compartment* to)
{
	to_.push_back(to);
}

int Reaction::perform(
	const long& nTimes,
	const string& strReaction,
	const double& time,
	map<string,vector<long>>& compTrajectories,
	unsigned indxTraj,
	unsigned leafcount,
	vector<Node*>* roots,
	bool isresampling,
	bool fullTraj)
{
	//try{
		int count = leafcount;
		// cout << "Leafcount Reaction.cpp 1 : " << strReaction << "  =  " << leafcount << endl;
		bool found = false;
		unsigned indexFrom;

		switch(type_){ // type de reaction à réaliser
			case birth: // naissance (ou coalescence)
				for(indexFrom = 0 ; indexFrom < from_.size() && !found ; indexFrom++){ //chrche le compartiment a l'origine de la transmission cad le donneur
					for(unsigned i = 0 ; i < to_.size() && !found ; i++){
						found = (from_[indexFrom]->getName() == to_[i]->getName());
					}
				}
				indexFrom--;
				if(!found){ //si aucun donneur -> enracinement
					if(from_[0]->getOldNodes() == 0){
						count = -1;
						warning("Error in rooting. No node left for rooting.");
					}
					else{
						performRooting(strReaction,time,roots);
					}
				}
				else{
					count = evalCoalescence(nTimes, indexFrom, strReaction, time, count, roots, compTrajectories, indxTraj, fullTraj);
				}
				break;
			case migration:
				if(fullTraj & (from_[0]->getOldNodes()==0) ){
					// si on veut simuler la phylogénie complète et qu'il y a 0 node à pour la mimgration
					count = -2; // type d'erreur
				} else{
					count = evalMigration(nTimes,strReaction,time, count, compTrajectories,indxTraj);
				}
				break;
			case death:
				// cout << "count Reaction.cpp death reac : " << count << endl;
				break;
			case sampledDeath:
				count = evalSampling(nTimes,strReaction,time,count,compTrajectories,indxTraj,false);
				break;
			case sampling:
				if( (to_[0]->getSize() - to_[0]->getOldNodes() - to_[0]->getOldUnsampledNodes()) >= nTimes ){
					// evalue si il reste des individus non echantillonnés
					count = evalSampling(nTimes,strReaction,time,count,compTrajectories,indxTraj,isresampling);
				}
				else{
					//warning("Error : Cannot sample compartment ", to_[0]->getName(), ", the number of individuals is not sufficient.");
					warning("Error : Cannot sample a compartment, the number of individuals is not sufficient.");
					count = -1;
				}
				break;
			case none:
				break;
			default:
				break;
		}
		return count;
}



int Reaction::performRooting(
	const string& strReaction,
	const double& time,
	vector<Node*>* roots)
{
		bool ok = true;
		int res = 0;
		Node* root = new Node("",time); //cree un nouveau noeud racine
		unsigned r1 = R::runif(0, from_[0]->getOldNodes()-1);
		// unsigned r1 = rand() % from_[0]->getOldNodes(); // pioche aléatoirement un noeud à enraciner
		root->addSon(from_[0]->popNode(r1)); //relie ce noeud au noeud racine
		ok = ok & from_[0]->decrementOldNodes();
		ok = ok & from_[0]->decrementSize();
		roots->push_back(root);
		if(!ok){
			res = -1;
		}
		return res;
}

unsigned Reaction::rhyper(const unsigned& nTimes, const unsigned& k, const unsigned& N){
	unsigned int count = 0;
	count = Rcpp::rhyper(1, k, (N-k), nTimes)[0];
	return count;
}

int Reaction::evalSampling(
		const long& nTimes,
		const string& strReaction,
		const double& time,
		const unsigned& leafcount,
		map<string,vector<long>>& compTrajectories,
		unsigned indxTraj,
		bool isresampling)
{
		int count = 0;
		bool ok = true;
		unsigned nbSampling = 0;	// nombre d'echantillonnages simples a realiser
		unsigned nbReSampling = 0;	// nombre de "re-"echantillonnage a realiser
		if( isresampling || (!isSampledDeath() && isresampling) ){
			nbReSampling = rhyper(nTimes, (double) to_[0]->getOldUnsampledNodes(), (to_[0]->getSize() - ((double) to_[0]->getOldNodes() - (double) to_[0]->getOldUnsampledNodes()) ) );
		}
		nbSampling = nTimes - nbReSampling;
		for(unsigned i = 0 ; i < nbSampling ; i++){
			ok = ok & performSampling(leafcount+i, strReaction, time); // realiser un simple echantillonnage
		}
		for (unsigned i = 0 ; i < nbReSampling && ok; i++){
			ok = ok & performReSampling(leafcount+i, strReaction, time); // realiser un "re-echantillonnage"
		}
		if(ok){
			count = leafcount+nTimes;
		} else{
			count = -1;
		}

	return count;
}

int Reaction::evalCoalescence(
		const long& nTimes,
		const unsigned& indexFrom,
		const string& strReaction,
		const double& time,
		const unsigned& leafcount,
		vector<Node*>* roots,
		map<string,vector<long>>& compTrajectories,
		unsigned indxTraj,
		bool fullTraj)
{
		int res = leafcount;
		bool ok = true ;

		unsigned nbCoal = 0;
		unsigned nbInviCoal = 0;

		unsigned nbDonnor = 0;
		unsigned nbRecipient = rhyper( nTimes, (double) from_[1-indexFrom]->getOldNodes(), from_[1-indexFrom]->getSize());

		if(from_[0] != from_[1]){ // si les deux individus obtenus par transmission appartiennent a 2 compartiments différents (ex : E ET I)
			nbDonnor = rhyper( nTimes, (double) from_[indexFrom]->getOldNodes(), from_[indexFrom]->getSize());
		}
		else{ // si le donneur genere un individu du même compartiment que lui
			nbDonnor = rhyper( nTimes, (double) from_[indexFrom]->getOldNodes() - nbRecipient, ( from_[1-indexFrom]->getSize() - nTimes));
		}

		nbCoal = rhyper((double) nbRecipient, (double) nbDonnor, (double) nTimes);
		nbInviCoal = nbRecipient - nbCoal;

		for(unsigned i = 0 ; i < nbCoal ; i++){
			ok = ok & performCoalescence(indexFrom, strReaction, time); // realiser une coalescence
		}
		for(unsigned i = 0 ; i < nbInviCoal ; i++){
			ok = ok & performInvisibleCoalescence(indexFrom, strReaction, time); // realiser une coalescence invisible
		}
		if(!ok){
			res = -1;
		}

		if(fullTraj & ((nbInviCoal + nbCoal) == 0) ){
			// dans les cas où il n'y a pas eu de feuilles crées : si coalescence -> on passe à la prochaine reaction
			res = -2; // type de warning -> continue dans la simulation
		}

		return res;
}

bool Reaction::performReSampling(
		const unsigned& leafcount,
		const string& strReaction,
		const double& time)
{
	bool ok = true;
	if(to_[0]->getOldUnsampledNodes() == 0){
		ok = false;
	    warning("Error in re-sampling. No nodes left.");
	} else{
		// cout << "resampling reaction" << endl;
		stringstream name, infos;
		name << strReaction << "_" << leafcount;
		// cout << "name of resampled node : " << name.str() << endl;
		Node* leaf = new Node(name.str(), time); // cree une nouvelle feuille
		infos << leaf->getInfos(); // ajoute des infos sur la feuille
		if (infos.str().size()>0)
			infos << ",";
		infos << "reaction_string=\"" << strReaction << "\"";
		infos << ",reaction_type=\"sampling\"";
		leaf->setInfos(infos.str());
		leaf->setIsSampled(true);
		Node* innerNode = new Node("", time); // cree un noeud interne pour connecter un noeud de l'arbre à la nouvelle feuille
		infos << ",reaction_specification=\"re-sampling\"";
		innerNode->setInfos(infos.str());
		innerNode->setIsSampled(true);
		unsigned r1 = R::runif(0,to_[0]->getOldUnsampledNodes()-1);
		// unsigned r1 = rand() % to_[0]->getOldUnsampledNodes(); // pioche aleatoirement un noeud de l'arbre correspondant a un individu non echantillonné
		innerNode->addSon(to_[0]->popNonSampledNode(r1)); // relie ce noeud au nouveau noeud interne
		ok = ok & to_[0]->decrementOldUnsampledNodes();
		ok = ok & to_[0]->decrementOldNodes();
		innerNode->addSon(leaf); // relie le nouveau noeud interne a la feuille
		to_[0]->addNode(innerNode);
		ok = ok & to_[0]->incrementNewNodes();
	}

	return ok;

}

bool Reaction::performSampling(
		const unsigned& leafcount,
		const string& strReaction,
		const double& time)
{
		// cout << "sampling reaction" << endl;
		bool ok = true ;
		stringstream name, infos;

		if(isSampledDeath()){
			name << to_[0]->getName() << "_" << count_ ;//leafcount;
		} else{
			name << strReaction << "_" << count_; //leafcount;
		}
		// cout << "name of sampled node : " << name.str() << endl;
		Node* leaf = new Node(name.str(), time); // cree une nouvelle feuille
		infos << leaf->getInfos(); // ajoute des infos
		if (infos.str().size()>0)
			infos << ",";
		infos << "reaction_string=\"" << strReaction << "\"";
		infos << ",reaction_type=\"sampling\"";
		leaf->setInfos(infos.str());
		leaf->setIsSampled(true);
		to_[0]->addNode(leaf);	// repertorie la feuille
		ok = ok & to_[0]->incrementNewNodes();

		count_++;

		return ok;
}

bool Reaction::performCoalescence(
		const unsigned& indexFrom,
		const string& strReaction,
		const double& time)
{
		bool ok = true;
		Node* innerNode = new Node("", time); // créé un nouveau noeud interne

		stringstream infos;
		infos << innerNode->getInfos(); // ajoute des infos sur le noeud au cas où on voudrait les imprimer dans l'arbre
		if (infos.str().size()>0)
			infos << ",";
		infos << "reaction_string=\"" << strReaction << "\"";
		infos << ",reaction_type=\"birth\"";
		infos << ",reaction_specification=\"coalescence\"";
		innerNode->setInfos(infos.str());

		unsigned rRecipient = R::runif(0,from_[1-indexFrom]->getOldNodes()-1);
		// unsigned rRecipient = rand() % from_[1-indexFrom]->getOldNodes(); // pioche le nœud receveur
		innerNode->addSon(from_[1-indexFrom]->popNode(rRecipient));
		ok = ok & from_[1-indexFrom]->decrementOldNodes();
		ok = ok & from_[1-indexFrom]->decrementSize();

		unsigned rDonnor = R::runif(0,from_[indexFrom]->getOldNodes()-1);
		// unsigned rDonnor = rand() % from_[indexFrom]->getOldNodes(); //pioche le nœud donneur
		innerNode->addSon(from_[indexFrom]->popNode(rDonnor));
		ok = ok & from_[indexFrom]->decrementOldNodes();
		from_[indexFrom]->addNode(innerNode); // coalesce les deux noeuds à leur noeud parent
		ok = ok & from_[indexFrom]->incrementNewNodes();

		return ok;
}

bool Reaction::performInvisibleCoalescence(
		const unsigned& indexFrom,
		const string& strReaction,
		const double& time)
{
		bool ok = true;
		// create new inner node
		Node* innerNode = new Node("", time); // créé un nouveau noeud interne

		// generate node infos
		stringstream infos; // ajoute des infos sur le noeud au cas où on voudrait les imprimer dans l'arbre
		infos << innerNode->getInfos();
		if (infos.str().size()>0)
			infos << ",";
		infos << "reaction_string=\"" << strReaction << "\"";
		infos << ",reaction_type=\"birth\"";
		infos << ",reaction_specification=\"invisible coalescence\"";
		innerNode->setInfos(infos.str());

		// select a node to relocate
		unsigned rRecipient = R::runif(0,from_[1-indexFrom]->getOldNodes()-1);
		// unsigned rRecipient = rand() % from_[1-indexFrom]->getOldNodes();
		// link it to the new inner node
		innerNode->addSon(from_[1-indexFrom]->popNode(rRecipient));
		ok = ok & from_[1-indexFrom]->decrementOldNodes();
		ok = ok & from_[1-indexFrom]->decrementSize(); // diminue le nombre de noeud disponibles
		from_[indexFrom]->addNode(innerNode);
		ok = ok & from_[indexFrom]->incrementNewNodes();

		return ok;
}

int Reaction::evalMigration(
		const long& nTimes,
		const std::string& strReaction,
		const double& time,
		const unsigned& leafcount,
		map<string,vector<long>>& compTrajectories,
		unsigned indxTraj)
{
		int res = leafcount;
		bool ok = true ;
		unsigned nbMigration = 0;
		nbMigration = rhyper(nTimes, (double)from_[0]->getOldNodes(), from_[0]->getSize());
		for(unsigned i = 0 ; i < nbMigration ; i++){
			ok = ok & performMigration(strReaction, time);
		}
		if(!ok){
			res = -1;
		}
		//from_[0]->setSize(from_[0]->getSize()-(nTimes-nbMigration));
		return res;
}

bool Reaction::performMigration(
		const string& strReaction,
		const double& time)
{
		bool ok = true;
		Node* innerNode = new Node("", time); // cree un nouveau noeud interne
		stringstream infos;
		infos << innerNode->getInfos();
		if (infos.str().size()>0)
			infos << ",";
		infos << "reaction_string=\"" << strReaction << "\"";
		infos << ",reaction_type=\"migration\"";
		innerNode->setInfos(infos.str());
		unsigned r1 = R::runif(0,from_[0]->getOldNodes()-1);
		// unsigned r1 = rand() % from_[0]->getOldNodes(); // pioche aleatoirement un noued a migrer
		innerNode->addSon(from_[0]->popNode(r1)); // relie ce noeud au nouveau noeud interne
		ok = ok & from_[0]->decrementOldNodes();
		ok = ok & from_[0]->decrementSize();
		to_[0]->addNode(innerNode); // ajoute le nouveau noeud interne au bon compartiment
		ok = ok & to_[0]->incrementNewNodes();

		return ok;
}
