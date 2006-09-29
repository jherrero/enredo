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
    bool populate_from_file(char *filename, float min_score, int max_gap_length);
    void minimize();
    void print_anchors_histogram(void);
    void print_stats(int histogram_size);
    void print_links(uint min_anchors, uint min_regions, uint min_length);
    void merge_alternative_paths(uint max_anchors);

protected:
    std::map<std::string, Anchor*> anchors;
    std::map<std::string, std::string*> species;
    std::map<std::string, std::string*> chrs;
};

#endif
