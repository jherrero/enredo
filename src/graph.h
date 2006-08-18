#ifndef GRAPH_H
#define GRAPH_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/
#include <cstdlib>
#include <map>

typedef class Anchor Anchor;

class Graph{
public:
    Graph();

    ~Graph();
    void add_Anchor(Anchor *this_anchor);
    Anchor* get_Anchor(std::string id);
    bool populate_from_file(char *filename);
    void minimize();
    void print();
    void print_links(int min_anchors, int min_regions, int min_length);
    void merge_alternative_paths(int max_anchors);

protected:
    std::map<std::string, Anchor*> anchors;
    std::map<std::string, std::string*> species;
    std::map<std::string, std::string*> chrs;
};

#endif
