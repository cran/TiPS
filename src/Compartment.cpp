#include "Compartment.h"

using namespace std;

Compartment::Compartment() :
		name_(""),
		oldNodes_(0),
		newNodes_(0),
		oldUnsampledNodes_(0),
		nodes_(),
		size_(0)
{}

Compartment::Compartment(const std::string& name) :
		name_(name),
		oldNodes_(0),
		newNodes_(0),
		oldUnsampledNodes_(0),
		nodes_(),
		size_(0)
{}

Compartment::Compartment(const std::string& name, const long& i_size) :
		name_(name),
		oldNodes_(0),
		newNodes_(0),
		oldUnsampledNodes_(0),
		nodes_(),
		size_(i_size)
{}

Compartment::~Compartment(){}

void Compartment::addNode(Node* node){
	nodes_.push_back(node);
}

void Compartment::insertNode(const unsigned& index, Node* node){
	nodes_.insert(nodes_.begin()+index, node);
}

Node* Compartment::popNode(const long& index){
	Node* tmpNode = nodes_[index];
	nodes_.erase(nodes_.begin()+index);
	return tmpNode;
}

bool Compartment::replaceNode(const unsigned& index, Node* node){
	bool asw = (index >= 0 && index < nodes_.size());
	if (asw){
		nodes_[index] = node;
	}
	return asw;
}

bool Compartment::decrementOldNodes(){
	bool ok = true;
	if(oldNodes_ < 1){
		warning("Error: Compartment, variable oldNodes_ cannot have a negative size.");
	    ok = false;
	} else{
		oldNodes_--;
	}
	return ok ;
}

bool Compartment::incrementNewNodes(){
	bool ok = true;
	if (newNodes_ >= std::numeric_limits<long>::max()){
		warning("Error: Compartment, variable newNodes_ has reached the maximal LONG value.");
		ok = false;
	} else{
		newNodes_++;
	}
	return ok ;
}

Node* Compartment::popNonSampledNode(const long& index){
	unsigned i;
	unsigned count = 0;
	for (i=0 ; i<nodes_.size() && count<=index ; i++){
		if(!nodes_[i]->isSampled())
			count ++;
	}
	i--;
	Node* tmpNode = nodes_[i];
	nodes_.erase(nodes_.begin()+i);
	return tmpNode;
}

bool Compartment::updateNodes(){
	bool ok = true ;
	oldNodes_ = oldNodes_ + newNodes_;
	newNodes_ = 0;
	oldUnsampledNodes_ = 0;
	for (unsigned i=0 ; i<oldNodes_ ; i++){
		oldUnsampledNodes_ += (!nodes_[i]->isSampled());
	}
	if (oldNodes_ > size_){
		warning("Error: In compartment, variable activeNodes_ greater than compartment size.");
		ok = false ;
	}
	if (oldUnsampledNodes_ > oldNodes_){
		warning("Error: In compartment, variable activeUnsampledNodes_ greater than activeNodes_ .");
		ok = false ;
	}
	return ok ;
}

void Compartment::update(){
	oldNodes_ = oldNodes_ + newNodes_;
	newNodes_ = 0;
	oldUnsampledNodes_ = 0;
	for (unsigned i=0 ; i<oldNodes_ ; i++){
		oldUnsampledNodes_ += (!nodes_[i]->isSampled());
	}
}

bool Compartment::decrementOldUnsampledNodes(){
	bool ok = true;
	if (oldUnsampledNodes_ < 1){
		warning("Error: Compartment, variable oldUnsampledNodes_ cannot have a negative size.");
		ok = false;
	} else{
		oldUnsampledNodes_--;
	}
	return ok ;
}

bool Compartment::decrementSize(){
	bool ok = true;
	if (size_ < 1){
		warning("Error: In compartment, variable size_ cannot have a negative size.");
		ok = false;
	} else{
		size_--;
	}
	return ok ;
}

void Compartment::incrementOldNodes(){
	oldNodes_++;
}
