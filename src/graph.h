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

//! A Graph is made of Anchor objects linked by Links. Each Anchor is a vertex and each Link is an edge

class Graph{
public:
    Graph();

    ~Graph();
    //! Adds an anchor in the graph
    void add_Anchor(Anchor *this_anchor);
    Anchor* get_Anchor(std::string id);
    bool populate_from_file(char *filename, float min_score, int max_gap_length, bool anchors_as_links);
    void minimize(std::string debug = "");
    void print_anchors_histogram(std::ostream &out = std::cout);
    void print_stats(int histogram_size);
    unsigned long int print_links(std::ostream &out = std::cout, uint min_anchors = 1, uint min_regions = 1, uint min_length = 0);
    int merge_alternative_paths(uint max_anchors, uint max_length = 10000, std::string debug = "");
    void study_anchors(void);
    int simplify(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0, std::string debug = "");
    int simplify_aggressive(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0, std::string debug = "");
    int split_unselected_links(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0, std::string debug = "");
    void split_unbalanced_links(float max_ratio, std::string debug = "");
    //! Looks for small palindromes and destroy them: A=B=C will become A-B-C-B-A
    uint resolve_small_palindromes(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0, std::string debug = "");
    //! Looks for small insertions and assimilates them in the paths
    uint assimilate_small_insertions(uint min_anchors = 1, uint min_regions = 1, uint min_length = 0,
                                     uint max_insertion_length = 10000, std::string debug = "");

protected:
    std::map<std::string, Anchor*> anchors;
    std::map<std::string, std::string*> species;
    std::map<std::string, std::string*> chrs;
};

#endif
