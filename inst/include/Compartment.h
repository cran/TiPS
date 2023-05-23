#ifndef COMPARTMENT_H_
#define COMPARTMENT_H_

#include <cstdlib>		// EXIT_FAILURE
#include <iostream>		// cout, endl, ...
#include <sstream>		// stringstream, ...
#include <string>		// string

#include "Node.h"

#include<Rcpp.h>
using namespace Rcpp;


class Compartment {
	private:
		std::string name_;
		long oldNodes_;
		long newNodes_;
		long oldUnsampledNodes_;
		std::vector<Node*> nodes_;
		long size_;
		bool isDeme_;

	public:
		Compartment();
		Compartment(const std::string& name);
		Compartment(const std::string& name, const long& size);
		virtual ~Compartment();

		// Getters
		std::string getName() const {return name_;}
		std::vector<Node*> getNodes() const {return nodes_;}
		Node* getNode(const unsigned& index) {return nodes_[index];}
		long getNodeSize() const {return nodes_.size();}
		long getOldNodes() const {return oldNodes_;}
		long getNewNodes() const {return newNodes_;}
		long getOldUnsampledNodes() const {return oldUnsampledNodes_;}
		long getSize() const { return size_; };
		bool getIsDeme() const { return isDeme_; };

		// Setters
		void setName(const std::string& name) {name_=name;}
		void setNodes(const std::vector<Node*>& nodes) {nodes_=nodes;}
		void setSize(const long& size) { size_=size; };
		void setIsDeme(const bool& isDeme) { isDeme_ = isDeme; }

		// Modifiers
		void init();
		void update();
		bool updateNodes();
		void addNode(Node* node);
		void insertNode(const unsigned& index, Node* node);
		Node* popNode(const long& index);
		Node* popNonSampledNode(const long& index);
		bool replaceNode(const unsigned& index, Node* node);
		bool decrementOldNodes();
		void incrementOldNodes();
		bool incrementNewNodes();
		bool decrementOldUnsampledNodes();
		bool decrementSize();

		// // Printers
		// std::string toString() const;
};

#endif /* COMPARTMENT_H_ */
