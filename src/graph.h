#ifndef GRAPH_H
#define GRAPH_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <fstream>

typedef class Anchor Anchor;

class Graph{
public:
    Graph();

    ~Graph();
    void add_Anchor(Anchor *this_anchor);
    Anchor* get_Anchor(std::string id);
    bool populate_from_file(char *filename, float min_score, int max_gap_length);
    void minimize();
    void print_anchors_histogram(std::ostream &out = std::cout);
    void print_stats(int histogram_size);
    unsigned long int print_links(std::ostream &out = std::cout, uint min_anchors = 1, uint min_regions = 1, uint min_length = 0);
    void merge_alternative_paths(uint max_anchors);
    void Graph::study_anchors(void);
    void Graph::simplify(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0);
    void split_unbalanced_links(float max_ratio);

protected:
    std::map<std::string, Anchor*> anchors;
    std::map<std::string, std::string*> species;
    std::map<std::string, std::string*> chrs;
};

#endif
