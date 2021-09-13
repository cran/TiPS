#ifndef REACTION_H_
#define REACTION_H_

#include <map>
#include <numeric>		// accumulate
#include <vector>
#include <iomanip>		// setprecision

#include <stdlib.h>		//srand, rand

#include<Rcpp.h>
using namespace Rcpp;


// #include <boost/math/distributions/hypergeometric.hpp> // hypergeometric distribution

#include "Compartment.h"

class Reaction {
	private:
		enum ReactionType{birth,migration,death,sampledDeath,sampling,none};
		ReactionType type_;
		std::vector<Compartment*> from_;
		std::vector<Compartment*> to_;
		int count_ ;

	public:
		Reaction();
		virtual ~Reaction();

		// Getters
		bool isBirth() const {return type_ == birth;}
		bool isMigration() const {return type_ == migration;}
		bool isDeath() const {return type_ == death;}
		bool isSampledDeath() const {return type_ == sampledDeath;}
		bool isSampling() const {return type_ == sampling;}
		std::vector<Compartment*> getFromCompartements() const {return from_;}
		std::vector<Compartment*> getToCompartements() const {return to_;}

		void getInfos(std::string& strReaction, const double& tmpTime, unsigned indx);

		int getCount() const {return count_;}

		// Setters
		void setIsBirth() {type_=birth;}
		void setIsMigration() {type_=migration;}
		void setIsDeath(const bool& fullTree) {if (fullTree) {type_=sampledDeath;} else {type_=death;}}
		void setIsSampling() {type_=sampling;}
		void setFromCompartements(const std::vector<Compartment*>& from) {from_=from;}
		void setToCompartements(const std::vector<Compartment*>& to) {to_=to;}

		void setCount(const int& count) {count_ = count;}

		// Adders
		void addFrom(Compartment* from);
		void addTo(Compartment* to);

		int perform(const long& nTimes,const std::string& strReaction,const double& time,std::map<std::string,std::vector<long>>& compTrajectories,unsigned indxTraj,unsigned count, std::vector<Node*>* roots, bool isresampling, bool fullTraj);
		// void evalRooting(const long& nTimes, const std::string& strReaction, const double& time, std::vector<Node*>* roots, std::ofstream& outLog);
		int evalCoalescence(const long& nTimes, const unsigned& indexFrom, const std::string& strReaction, const double& time, const unsigned& leafcount, std::vector<Node*>* roots, std::map<std::string,std::vector<long>>& compTrajectories,unsigned indxTraj, bool fullTraj);
		int evalMigration(const long& nTimes, const std::string& strReaction, const double& time, const unsigned& leafcount, std::map<std::string,std::vector<long>>& compTrajectories,unsigned indxTraj);
		int evalSampling(const long& nTimes, const std::string& strReaction, const double& time, const unsigned& leafcount, std::map<std::string,std::vector<long>>& compTrajectories,unsigned indxTraj, bool isresampling);

		bool performCoalescence(const unsigned& indexFrom, const std::string& strReaction, const double& time);
		bool performInvisibleCoalescence(const unsigned& indexFrom, const std::string& strReaction, const double& time);
		bool performReSampling(const unsigned& leafcount, const std::string& strReaction, const double& time);
		bool performSampling(const unsigned& leafcount, const std::string& strReaction, const double& time);
		int performRooting(const std::string& strReaction, const double& time, std::vector<Node*>* roots);
		bool performMigration(const std::string& strReaction, const double& time);

		unsigned rhyper(const unsigned& nTimes, const unsigned& k, const unsigned& N);
};

#endif /* REACTION_H_ */
