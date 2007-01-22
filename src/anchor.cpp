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
    \fn Graph::minimize_anchor(Anchor* this_anchor)
 */
uint Anchor::minimize(void)
{
  uint count = 0;
  bool merge_event;
  do {
    merge_event = false;
    for (std::list<Link*>::iterator p_link1 = this->links.begin(); p_link1 != this->links.end() and !merge_event; p_link1++) {
      short strand1;
      if ((*p_link1)->anchor_list.front() == (*p_link1)->anchor_list.back()) {
        strand1 = 0;
      } else if ((*p_link1)->anchor_list.back() == this) {
        strand1 = 1;
      } else if ((*p_link1)->anchor_list.front() == this) {
        strand1 = -1;
      }
      for (std::list<Link*>::iterator p_link2 = this->links.begin(); p_link2 != p_link1 and !merge_event; p_link2++) {
        if ((*p_link2)->tags.size() != (*p_link1)->tags.size()) {
          continue;
        }
        short strand2;
        if ((*p_link2)->anchor_list.front() == (*p_link2)->anchor_list.back()) {
          strand2 = 0;
        } else if ((*p_link2)->anchor_list.back() == this) {
          strand2 = -1;
        } else if ((*p_link2)->anchor_list.front() == this) {
          strand2 = 1;
        }
        if ((*p_link1)->try_to_concatenate_with(*p_link2, strand1, strand2)) {
          merge_event = true;
          count++;
        }
      }
    }
  } while (merge_event);

  return count;
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
