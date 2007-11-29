#ifndef LINK_H
#define LINK_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/

#include <string>
#include <list>
#include <vector>

using namespace std;

typedef class Anchor Anchor;

//! A tag in a Link represents a genomic region which goes through the path of the Link

struct tag {
  string *species;
  string *chr;
  uint start;
  uint end;
  short strand; // 1 when tag start => end corresponds to anchor_list.front() => anchor_list.back()
                // and -1 when start => end corresponds to anchor_list.back() => anchor_list.front()
};
    void print_tag(tag this_tag, ostream &out = cout);

//! A Link defines an edge in the Enredo graph

class Link{
public:
  Link(Anchor *anchor1, Anchor *anchor2);
  Link(Link *my_link);

    ~Link();
    void add_tag(string *species, string *chr, int start, int end, short strand);
    Link* merge(Link* other_link);
    bool try_to_concatenate_with(Link *other_link, short strand1 = 0, short strand2 = 0);
    void reverse();
    void print(ostream &out = cout);
    uint get_shortest_region_length();
    uint get_longest_region_length();
    bool is_an_alternative_path_of(Link* other_link);
    uint get_num_of_mismatches(Link* other_link);
    std::vector< std::list<tag>::iterator > get_matching_tags(Link *other_link, short strand1 = 0, short strand2 = 0,
                                                             bool allow_partial_match = false);
    Link* split(vector<bool> tags_to_split);
    Link* split(std::vector< std::list<tag>::iterator > tags_to_split);
    short get_strand_for_matching_tags(Anchor* anchor);
    bool is_valid(uint min_anchors, uint min_regions, uint min_length);
    bool is_bridge(uint min_anchors, uint min_regions, uint min_length, bool trim_link = true);

    list<Anchor*> anchor_list;

    list<tag> tags; //!< list of \link tag tags \endlink
};

#endif
