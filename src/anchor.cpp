#include "anchor.h"

Anchor::Anchor(string this_id)
{
  id = this_id;
  num = 1;
}


Anchor::~Anchor()
{
}


/*!
    \fn Anchor::get_direct_Link(Anchor *other_anchor)
 */
Link* Anchor::get_direct_Link(Anchor *other_anchor)
{
  for (std::list<Link*>::iterator it = links.begin(); it != links.end(); it++) {
    if ((*it)->anchor_list.size() != 2) {
      continue;
    }
    std::list<Anchor*>::iterator anchor_it = (*it)->anchor_list.begin(); 
    Anchor *anchor1 = (*anchor_it);
    anchor_it++;
    Anchor *anchor2 = (*anchor_it);
    if (anchor1->id == this->id) {
      if (anchor2 == other_anchor) {
        return *it;
      }
    } else if (anchor1 == other_anchor) {
      return *it;
    }
  }

  // other_anchor was not found => create a new Link
  Link* new_link = new Link(other_anchor, this);
  if (new_link) {
    this->add_Link(new_link);
    if (other_anchor != this) {
      other_anchor->add_Link(new_link);
    }
    return new_link;
  }
  cerr << "Cannot add a new Anchor" << endl;
  return NULL;
}


/*!
    \fn Anchor::add_Link(Link *link)
 */
void Anchor::add_Link(Link *link)
{
  this->links.push_back(link);
}


/*!
    \fn Anchor::print(ostream &out)
 */
void Anchor::print(ostream &out)
{
  out << "Anchor " << this->id << endl;
  for (std::list<Link*>::iterator p_link_it = this->links.begin(); p_link_it != this->links.end(); p_link_it++) {
    (*p_link_it)->print(out);
  }
}
