#ifndef LINK_H
#define LINK_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/

#include <string>
#include <list>

using namespace std;

typedef class Anchor Anchor;

struct tag {
  string *species;
  string *chr;
  uint start;
  uint end;
  short strand; // 1 when tag start => end corresponds to anchor_list.front() => anchor_list.back()
                // and -1 when start => end corresponds to anchor_list.back() => anchor_list.front()
};

class Link{
public:
  Link(Anchor *anchor1, Anchor *anchor2);

    ~Link();
    void add_tag(string *species, string *chr, int start, int end, short strand);
    Link* concatenate(Link* other_link);
    Link* merge(Link* other_link);
    bool can_be_concatenated_with(Link *other_link);
    void reverse();
    void print(ostream &out = cout);
    uint get_shortest_length();
    bool is_an_alternative_path_of(Link* other_link);
    uint get_num_of_mismatches(Link* other_link);

    list<Anchor*> anchor_list;

    list<tag> tags;
};

#endif
