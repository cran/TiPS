#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <tuple>

#include "Node.h"

using namespace std;

Node::Node() :
    name_(""),
    father_(0),
    sons_(),
    distanceToFather_(0),
    height_(0),
    bootstrap_(0),
    isSampled_(false),
	infos_("")
{}
Node::Node(const string& name) :
    name_(name),
    father_(0),
    sons_(),
    distanceToFather_(0),
    height_(0),
    bootstrap_(0),
    isSampled_(false),
	infos_("")
{}

Node::Node(const std::string& name, const double& height):
	name_(name),
	father_(0),
	sons_(),
	distanceToFather_(0),
	height_(height),
	bootstrap_(0),
	isSampled_(false),
	infos_("")
{
	stringstream ss;
	ss << "height=" << setprecision(10) << height;
	infos_ = ss.str();
}

Node::Node(const double& distanceToFather, const double& height) :
    name_(""),
    father_(0),
    sons_(),
    distanceToFather_(distanceToFather),
    height_(height),
    bootstrap_(0),
    isSampled_(false),
	infos_("")
{
	stringstream ss;
	ss << "height=" << setprecision(10) << height;
	infos_ = ss.str();
}
Node::Node(const string& name, const double& distanceToFather, const double& height) :
    name_(name),
    father_(0),
    sons_(),
    distanceToFather_(distanceToFather),
    height_(height),
    bootstrap_(0),
    isSampled_(false),
	infos_("")
{
	stringstream ss;
	ss << "height=" << setprecision(10) << height;
	infos_ = ss.str();
}

Node::~Node(){}

void Node::setFather(Node* node)
{
    father_ = node;
    if (find(node->sons_.begin(), node->sons_.end(), this) == node->sons_.end())
    {
        node->sons_.push_back(this);
    }
}

unsigned Node::getNbSons() const
{
    return sons_.size();
}
bool Node::isRoot() const
{
    return father_ == 0;
}
bool Node::isLeaf() const
{
    return sons_.size() == 0;
}

bool Node::getIsSampled(){
    return isSampled_;
}

void Node::addSon(Node* node)
{
    if (find(sons_.begin(), sons_.end(), node) == sons_.end())
    {
        sons_.push_back(node);
    }
    node->father_ = this;
}
string Node::newick(const bool& withInfos) const
{
    string ss = "";
    if (!isLeaf())
    {
        ss += "(";
        ss += (sons_[0])->newick(withInfos);
        for (unsigned i=1 ; i<sons_.size() ; i++)
            ss += ","+(sons_[i])->newick(withInfos);
        ss += ")";
    }
    if (name_ != "" && !isRoot())
    {
        //cout << "Nom : " << name_ << endl;
        ss += "\""+getName()+"\"";
    }
    if (! isRoot())
    {
        stringstream ss2;
        ss2 << "" << setprecision(10) << getDistanceToFather();
        if (withInfos)
        	ss +="[&"+getInfos()+"]";
        ss += ":"+ss2.str();
    }
    return ss;
}

vector<double> Node::getBranchLengths() const {

    if( isLeaf() )
    {
        vector<double> tmp;
        tmp.push_back(distanceToFather_);
        return  tmp ;
    }
    if (isRoot())
    {
        vector<double> tmp;
        for(unsigned i=0 ; i<sons_.size() ; i++)
        {
            vector<double> tmp2 = (sons_[i])->getBranchLengths();
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp ;
    }
    else
    {
        vector<double> tmp;
        tmp.push_back(distanceToFather_);
        for(unsigned i=0 ; i<sons_.size() ; i++)
        {
            vector<double> tmp2 = (sons_[i])->getBranchLengths();
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return  tmp ;
    }
}


void Node::initializeHeights()
{
    if (!isRoot())
    {
        height_ = father_->getHeight() + getDistanceToFather();
    }
    for (unsigned i=0 ; i<sons_.size() ; i++)
    {
        (sons_[i])->initializeHeights();
    }
}

void Node::initializeDistances()
{
    if (!isRoot())
    {
        distanceToFather_ = height_ - father_->getHeight();
    }
    for (unsigned i=0 ; i<sons_.size() ; i++)
    {
        (sons_[i])->initializeDistances();
    }
}

void Node::clean(){
    // retirer les fils temporaires/invisibles type coalesecnce invisible devient juste une branche
	//cout << sons_.size() << endl;
	if (isRoot())
	{
		for (unsigned i=0 ; i<sons_.size() ; i++)
		{
			(sons_[i])->clean();
		}
	}
	else if (!isLeaf())
    {
    	if (sons_.size()==1)
    	{
    		//cout << "deleting singleton" << endl;
    	    Node* tmpSon = sons_[0];
    	    Node* tmpFather  = father_;
    	    tmpFather->removeSon(this);
    		removeSon(tmpSon);
    	    tmpSon->setFather(tmpFather);
    		tmpFather->clean();
    	}
    	else
    	{
    	    for (unsigned i=0 ; i<sons_.size() ; i++)
    	    {
    	        (sons_[i])->clean();
    	    }
    	}
    }
}

string Node::getNewick(const bool& withInfos) const
{
    stringstream ss;
    ss << newick(withInfos) << ";";
    return ss.str();
}
void Node::removeSon(Node* node)
{
    for (unsigned i=0 ; i<sons_.size() ; i++)
    {
        if (sons_[i] == node)
        {
            sons_.erase(sons_.begin() + i);
            node->father_ = 0;
            break;
        }
    }
}
vector<Node*> Node::getLeaves()
{
    if (isLeaf())
    {
        vector<Node*> tmp;
        tmp.push_back(this);
        return tmp;
    }
    else
    {
        vector<Node*> tmp;
        for (unsigned i=0 ; i< sons_.size() ; i++)
        {
            vector<Node*> tmp2 = (sons_[i])->getLeaves();
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp;
    }
}

vector<int> Node::getDepths(int d)
{
    if (isLeaf())
    {
        vector<int> tmp;
        tmp.push_back(d);
        return tmp;
    }
    else
    {
        vector<int> tmp;
        tmp.push_back(d);
        for (unsigned i=0 ; i< sons_.size() ; i++)
        {
            vector<int> tmp2 = (sons_[i])->getDepths(d+1);
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp;
    }
}

int Node::getNbLeaves()
{
    int tmp = 0;
    if (isLeaf())
    {
        tmp++;
    }
    for (unsigned i=0 ; i< sons_.size() ; i++)
    {
        tmp += (sons_[i])->getNbLeaves();
    }
    id_ = tmp;
    // cout << "id tip : " << tmp << endl;
    return tmp;
}

int Node::setInnerNodesID(int nbtips, int maxid){
    if (!isLeaf()){
        maxid++;
        id_ = maxid + nbtips ;
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        maxid = (sons_[i])->setInnerNodesID(nbtips, maxid);
    }
    return maxid;
}

int Node::setLeavesID(int maxid){
    if (isLeaf()){
        maxid++;
        id_ = maxid;
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        maxid = (sons_[i])->setLeavesID(maxid);
    }
    return maxid;
}

pair<int,int> Node::getNbNodes() {
    pair<int,int> nodes ; //first is tips, second is inner nodes
    pair<int,int> tmpNodes ;

    if (isLeaf()){
        nodes.first++;
    }
    else {
        nodes.second++;
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        tmpNodes = sons_[i]->getNbNodes();
        nodes.first += tmpNodes.first;
        nodes.second += tmpNodes.second;
    }
    return nodes ;
}

map<string,vector<int>> Node::getEdges(){
    map<string,vector<int>> edges;
    map<string,vector<int>> tmp;
    if (isLeaf()){
        edges["to"].push_back(id_);
    }
    else if(!isLeaf() && !isRoot()){
        edges["from"].push_back(id_);
        edges["to"].push_back(id_);
    }
    else if (isRoot()){
        edges["from"].push_back(id_);
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        tmp = sons_[i]->getEdges();
        edges["from"].insert(edges["from"].end(),tmp["from"].begin(),tmp["from"].end());
        edges["to"].insert(edges["to"].end(),tmp["to"].begin(),tmp["to"].end());
        if (i == (sons_.size()-2) && !isLeaf()){
            edges["from"].push_back(id_);
        }
    }
    return edges;
}

vector<string> Node::getNodeLabels(){
    vector<string> nodes;

    if (!isLeaf()){
        nodes.push_back(infos_);
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        vector<string> tmp = sons_[i]->getNodeLabels();
        //cout << "tmp size : " << tmp.size() << endl;
        nodes.insert(nodes.end(),tmp.begin(),tmp.end());
    }

    return nodes;
}

vector<string> Node::getTipLabels(){
    vector<string> tips;

    if (isLeaf()){
        tips.push_back(name_);
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        vector<string> tmp = sons_[i]->getTipLabels();
        //cout << "tmp size : " << tmp.size() << endl;
        tips.insert(tips.end(),tmp.begin(),tmp.end());
    }
    //essaie annotation id pour les tips
    // id_ = tips.size() ;
    //cout << "tips size : " << tips.size() << endl;
    return tips;
}

vector<double> Node::getTipHeights(){
    vector<double> tips;

    if (isLeaf()){
        tips.push_back(height_);
    }
    for (unsigned i=0 ; i< sons_.size() ; i++){
        vector<double> tmp = sons_[i]->getTipHeights();
        //cout << "tmp size : " << tmp.size() << endl;
        tips.insert(tips.end(),tmp.begin(),tmp.end());
    }
    //essaie annotation id pour les tips
    // id_ = tips.size() ;
    //cout << "tips size : " << tips.size() << endl;
    return tips;
}

vector<Node*> Node::getSampledLeaves()
{
    if (isLeaf())
    {
        if (isSampled())
        {
            vector<Node*> tmp;
            tmp.push_back(this);
            return tmp;
        }
        else
        {
            vector<Node*> tmp;
            return tmp;
        }
    }
    else
    {
        vector<Node*> tmp;
        for (unsigned i=0 ; i< sons_.size() ; i++)
        {
            vector<Node*> tmp2 = (sons_[i])->getSampledLeaves();
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp;
    }
}
vector<Node*> Node::getInnerNodes()
{
    if (isLeaf())
    {
        vector<Node*> tmp;
        return tmp;
    }
    else
    {
        vector<Node*> tmp;
        tmp.push_back(this);
        for (unsigned i=0 ; i< sons_.size() ; i++)
        {
            vector<Node*> tmp2 = (sons_[i])->getInnerNodes();
            tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
        }
        return tmp;
    }
}

