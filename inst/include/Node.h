#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>		// setprecision
#include <map>
#include <stdlib.h>

#include<Rcpp.h>
using namespace Rcpp;


class Node
{
    private:
        // Attributs
        std::string name_; // nom du noeud
        Node* father_; // pointeur vers le noeud père
        std::vector<Node*> sons_; // liste de pointeurs vers les noeuds fils
        double distanceToFather_; // distance au père
        double height_; // profondeur (distance à la racine ou somme des longueur de branches entre la racine et ce noeud)
        double bootstrap_;
        bool isSampled_; // est-ce que le noeud est échantillonné (est-ce qu'il est présent dans la phylogénie)
        std::string infos_;
        int id_;

    public:
        // Constructeurs
        Node();
        Node(const std::string& name);
        Node(const std::string& name, const double& height);
        Node(const double& distanceToFather, const double& height);
        Node(const std::string& name, const double& distanceToFather, const double& height);

        // Destructeur
        ~Node();

        // getters
        std::string getName() const { return name_; };
        Node* getFather() { return father_; };
        std::vector<Node*>& getSons() { return sons_; };
        double getDistanceToFather() const { return distanceToFather_; };
        double getHeight() const { return height_; };
        double getBootstrap() const { return bootstrap_; };
        bool isSampled() const { return isSampled_; };
        std::string getInfos() const {return infos_;}
        int getID() const {return id_;}

        // setters
        void setName(const std::string& name) { name_ = name; }
        void setDistanceToFather(const double& distance) { distanceToFather_ = distance; }
        void setHeight(const double& height) { height_ = height; }
        void setBootstrap(const double& bootstrap) { bootstrap_ = bootstrap; }
        void setIsSampled(const bool& iSampled) { isSampled_ = iSampled; }
        void setFather(Node* node);
        void addSon(Node* node);
        void setInfos(const std::string& infos) {infos_=infos;}
        void setID(const int& id) {id_=id;}

        // Autres
        unsigned getNbSons() const; // retourne le nombre de fils du noeud courant
        bool isRoot() const; // vrai si le noeud courant n'a pas de noeud père
        bool isLeaf() const; // vrai si la liste de fils est vide
        void removeSon(Node* node); // supprime la connection entre le noeud courant et un noeud fils
        void clean();
        bool getIsSampled();

        // Méthodes liées aux entrées
        void initializeHeights(); // initialise les attributs height_ de tous les noeuds présents sous le noeud courant (méthode a utiliser après readAndBuild)
        void initializeDistances(); // initialise les attributs distanceToFather_ de tous les noeuds présents sous le noeud courant (méthode a utiliser après readAndBuild)

        // Méthodes liées aux sorties
        std::string newick(const bool& withInfos) const; // construit la chaine newick par récursivité
        std::string getNewick(const bool& withInfos) const; // retourne la chaine newick terminée

        // Méthodes pour créer un arbre de classe phylo de R (ape)
        std::vector<double> getBranchLengths() const; // retourne la liste de toutes les longueurs de branches du sous-arbre constitué par le noeud courant
        std::vector<Node*> getLeaves(); // liste des feuilles du sous-arbre constitué par le noeud courant
        std::vector<int> getDepths(int d); // liste des profondeurs des noeuds du sous-arbre constitué par le noeud courant
        int getNbLeaves(); // nombre de feuilles du sous-arbre constitué par le noeud courant
        //std::vector<std::tuple<double,int>> createEvents(); // méthode générant une liste de tuples <date,événement> en associant dates (profondeurs) et événements ( (1) = embranchement / (-1) = feuille)
        //static bool compareTuples(const std::tuple<double,int>& t1, const std::tuple<double,int>& t2); // méthode permettant de comparer deux tuples <date,événement> afin de les trier
        //std::tuple<std::vector<double>,std::vector<int>> getEvents(); // méthode retournant un vecteur trié de tuples <date,événement>
        std::vector<Node*> getSampledLeaves(); // retourne une liste de pointeurs vers les feuilles échantillonnées
        std::vector<Node*> getInnerNodes(); // retourne une liste de pointeurs vers les noeud internes


        std::pair<int,int> getNbNodes();
        std::vector<std::string> getTipLabels();
        std::vector<double> getTipHeights();
        std::vector<std::string> getNodeLabels();
        std::map<std::string,std::vector<int>> getEdges();
        int setInnerNodesID(int nbtips, int maxid); // annotate inner nodes with id
        int setLeavesID(int maxid); // annotate tips with id


};
#endif // NODE_H
