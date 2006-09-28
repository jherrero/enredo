#ifndef ANCHOR_H
#define ANCHOR_H

/**
	@author Javier Herrero <jherrero@ebi.ac.uk>
*/

#include <iostream>
#include <string>
#include <list>
#include "link.h"
using namespace std;


class Anchor{
public:
    Anchor(string id);

    ~Anchor();
    Link* get_direct_Link(Anchor *other_anchor);
    void add_Link(Link *link);
    void print();

    string id;
    uint num;
    std::list<Link*> links;

  protected:
};

#endif
