#ifndef PHYLOEPID_H_
#define PHYLOEPID_H_

#include "Compartment.h"
#include "Reaction.h"
#include <chrono>
#include <sys/time.h>
#include <random>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

class Phyloepid {
	private:
		// Settings* sets_;
		std::map<std::string,Compartment*> compartments_;
		std::map<std::string,Reaction*> reactions_;
		std::vector<Node*> roots_;
		std::map<std::string, std::vector<long>> compTrajectories_;
		int leafcount_;
		bool fullTraj_;
		List traj_;
		Rcpp::String outTree_;
		std::mt19937 randomGenerator_;
		int nTrials_;
		std::vector<std::string> compartmentNames_;

		std::vector<std::string> strReactions_;
		bool reSampling_; // option permettant transmission malgré avoir été échantillonné
		unsigned nbdates_;
		bool verbose_;
		long seed_;

	public:
		// Phyloepid(Settings* sets);
		Phyloepid(List reactions, List traj, bool fulTree, bool isresampling, unsigned int nbdates, bool verbose, NumericVector options);
		bool simulationTree();
		std::string getNexusTree(const bool& withInfos);
		std::string getNewickTree(const bool& withInfos);
		void initRandomSeed();
		virtual ~Phyloepid();
		void getInfo();

		// std::map<std::string,Compartment*> getCompartments() const {return compartments_;}
		// std::map<std::string,Reaction*> getReactions() const {return reactions_;}
		// std::vector<Node*> getRoots() const {return roots_;}

		// void setCompartments(const std::map<std::string,Compartment*>& compartments) {compartments_=compartments;}
		// void setReactions(const std::map<std::string,Reaction*>& reactions) {reactions_=reactions;}
		// void setRoots(const std::vector<Node*>& roots) {roots_=roots;}

		bool getIsReSampling(){return reSampling_;}

		void initCompartments();
		void readInitialStates(List initials);
		void readReactions(List reactions);
		void parseFrom(const std::string& from, Reaction* reaction);
		void parseTo(const std::string& to, Reaction* reaction);

		// std::string compartmentsToString(const bool& initial);
		// std::string compartmentsToStringAll(const bool& initial);
		// std::string reactionsToString();
		long sumUnrootedNodes();
		// void initCompartments();
		void updateCompartments();
		bool updateDemeCompartments(unsigned indxTraj);

		bool run();

		// Créer objet R arbre
		std::vector<double> getEdgeLengths() const; // retourne la liste de toutes les longueurs de branches du sous-arbre constitué par le noeud courant
		int getNbTips() const; // retourne le nombre de feuilles
		std::vector<std::string> getTipLabels() const;
		std::vector<int> getNbNodes() const;

		List createTreeObject() const; //créer la liste pour la phylogénie
};

#endif /* PHYLOEPID_H_ */
