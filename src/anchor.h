#ifndef ANCHOR_H
#define ANCHOR_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/

#include <iostream>
#include <string>
#include <list>
#include <set>
#include "link.h"
using namespace std;

//! Defines each vertex in the Graph

class Anchor{
public:
    Anchor(string id);

    ~Anchor();
    //! Returns the Link between this and another Anchor. Creates it if required
    Link* get_direct_Link(Anchor *other_anchor);
    //! Add a new Link at the end of the links list
    void add_Link(Link *link);
    //! Print the content of this Anchor
    void print(ostream &out = cout);
    //! Tries to transform (A=B) and (B=C) to (A=B=C) where B is this Anchor and A=B and B=C are compatible
    uint minimize(bool debug = false);

    string id; //!< the name of the Anchor as defined in the input file
    uint num; //!< the number of times this Anchor has been found in the input file
    std::list<Link*> links; //!< list of Link objects starting or ending in this Anchor
    std::set<std::string*> species; //!< sorted set of unique species in which this Anchor has been found

  protected:
};

#endif
