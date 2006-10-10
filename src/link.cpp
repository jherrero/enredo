#include <iostream>
#include "link.h"
#include "anchor.h"
#include <vector>

Link::Link(Anchor* anchor1, Anchor* anchor2)
{
  this->anchor_list.push_back(anchor1);
  this->anchor_list.push_back(anchor2);
//   cerr << "New link" << anchor1->id << ":" << anchor2->id << endl;
}


Link::~Link()
{
}


/*!
    \fn Link::add_tag(string species, string chr, int start, int end)
 */
void Link::add_tag(string *species, string *chr, int start, int end, short strand)
{
  tag this_tag;
  this_tag.species = species;
  this_tag.chr = chr;
  this_tag.start = start;
  this_tag.end = end;
  this_tag.strand = strand;

  tags.push_back(this_tag);
}


/*!
    \fn Link::concatenate(Link* other_link)
 */
Link* Link::concatenate(Link* other_link)
{
  // Get middle anchor
  Anchor *middle_anchor = NULL;
  if (this->anchor_list.back() == other_link->anchor_list.front() or this->anchor_list.back() == other_link->anchor_list.back()) {
    middle_anchor = this->anchor_list.back();
  } else if (this->anchor_list.front() == other_link->anchor_list.front() or this->anchor_list.front() == other_link->anchor_list.back()) {
    middle_anchor = this->anchor_list.front();
    this->reverse();
  } else {
    cerr << "Error while concatenating links!" << endl;
    return NULL;
  }
  if (middle_anchor != other_link->anchor_list.front()) {
    other_link->reverse();
  }
//   cout << "Middle anchor is " << middle_anchor->id << endl;

  // Concatenate tags
  bool * other_tag_has_been_concatenated;
  other_tag_has_been_concatenated = new bool [other_link->tags.size()];
  for (uint a = 0; a < other_link->tags.size(); a++) {
    other_tag_has_been_concatenated[a] = false;
  }
  for (list<tag>::iterator p_tag1 = this->tags.begin(); p_tag1 != this->tags.end(); p_tag1++) {
    bool this_tag_has_been_concatenated = false;
    int a = 0;
    for (list<tag>::iterator p_tag2 = other_link->tags.begin(); p_tag2 != other_link->tags.end(); p_tag2++) {
      if (p_tag1->species == p_tag2->species and p_tag1->chr == p_tag2->chr and p_tag1->start < p_tag2->end and p_tag2->start < p_tag1->end) {
        if (other_tag_has_been_concatenated[a] or p_tag1->strand * p_tag2->strand == -1) {
          continue;
        }
        if (p_tag1->strand == 0) {
          p_tag1->strand = p_tag2->strand;
        }
        if (p_tag1->start < p_tag2->start and p_tag1->end < p_tag2->end) {
          p_tag1->end = p_tag2->end;
        } else if (p_tag1->start > p_tag2->start and p_tag1->end > p_tag2->end) {
          p_tag1->start = p_tag2->start;
        } else {
          cerr << "Error while concatenating links! (1)" << endl;
          exit(2);
        }
        this_tag_has_been_concatenated = true;
        other_tag_has_been_concatenated[a] = true;
        break;
      }
      a++;
    }
    if (!this_tag_has_been_concatenated) {
      cerr << "Error while concatenating links! (2)" << endl;
      exit(2);
    }
  }
  for (uint a = 0; a < other_link->tags.size(); a++) {
    if (!other_tag_has_been_concatenated[a]) {
      cerr << "Error while concatenating links! (3)" << endl;
      exit(2);
    }
  }

  // Append anchors from other_link to this link
  for (list<Anchor*>::iterator anchor_it = ++other_link->anchor_list.begin(); anchor_it != other_link->anchor_list.end(); anchor_it++) {
//     cout << " Appending anchor " << (*anchor_it)->id << endl;
    this->anchor_list.push_back(*anchor_it);
  }

  if (this->anchor_list.front() != middle_anchor) {
    // Remove this and other_link from middle anchor
    for (list<Link*>::iterator p_link_it = middle_anchor->links.begin(); p_link_it != middle_anchor->links.end(); p_link_it++) {
      while (*p_link_it == this or *p_link_it == other_link) {
        p_link_it = middle_anchor->links.erase(p_link_it);
      }
    }
  } else {
    // Remove other_link only from middle anchor
    for (list<Link*>::iterator p_link_it = middle_anchor->links.begin(); p_link_it != middle_anchor->links.end(); p_link_it++) {
      if (*p_link_it == other_link) {
        p_link_it = middle_anchor->links.erase(p_link_it);
        break;
      }
    }
  }

  // Remove other_link from last anchor
  if (this->anchor_list.front() != this->anchor_list.back()) {
    for (list<Link*>::iterator p_link_it = other_link->anchor_list.back()->links.begin(); p_link_it != other_link->anchor_list.back()->links.end(); p_link_it++) {
      if (*p_link_it == other_link) {
        *p_link_it = this;
        break;
      }
    }
    delete(other_link);
  } else {
    for (list<Link*>::iterator p_link_it = other_link->anchor_list.back()->links.begin(); p_link_it != other_link->anchor_list.back()->links.end(); p_link_it++) {
      if (*p_link_it == other_link) {
        other_link->anchor_list.back()->links.erase(p_link_it);
        break;
      }
    }
  }

  return this;
}


/*!
    \fn Link::can_be_concatenated_with(Link *other_link)
 */
bool Link::can_be_concatenated_with(Link *other_link)
{
  std::vector<bool> other_tag_can_be_concatenated(this->tags.size(), false);
  for (list<tag>::iterator p_tag1 = this->tags.begin(); p_tag1 != this->tags.end(); p_tag1++) {
    bool this_tag_can_be_concatenated = false;
    int other_count = 0;
    for (list<tag>::iterator p_tag2 = other_link->tags.begin(); p_tag2 != other_link->tags.end(); p_tag2++) {
      if (p_tag1->species == p_tag2->species and p_tag1->chr == p_tag2->chr and p_tag1->start < p_tag2->end and p_tag2->start < p_tag1->end) {
        this_tag_can_be_concatenated = true;
        other_tag_can_be_concatenated[other_count] = true;
        continue;
      }
      other_count++;
    }
    if (!this_tag_can_be_concatenated) {
      return false;
    }
  }
  for (uint a = 0; a < other_link->tags.size(); a++) {
    if (!other_tag_can_be_concatenated[a]) {
      return false;
    }
  }
  return true;
}


/*!
    \fn Link::reverse()
 */
void Link::reverse()
{
  this->anchor_list.reverse();
  for (list<tag>::iterator p_tag = this->tags.begin(); p_tag != this->tags.end(); p_tag++) {
    // from 1 to -1; from -1 to 1 and from 0 to 0
    p_tag->strand *= -1;
  }
}


/*!
    \fn Link::print(ostream &out)
 */
void Link::print(ostream &out)
{
//   if (this->tags.size() == 1 or this->anchor_list.size() < 3) {
//     return;
//   }
  for (std::list<Anchor*>::iterator p_anchor_it = this->anchor_list.begin(); p_anchor_it != this->anchor_list.end(); p_anchor_it++) {
    out << " - " << (*p_anchor_it)->id;
  }
  out << "  (has " << this->tags.size() << " tags)" << endl;
  for (list<tag>::iterator p_tag = this->tags.begin(); p_tag != this->tags.end(); p_tag++) {
    out << "       " << *p_tag->species << ":" << *p_tag->chr << ":"
        << p_tag->start << ":" << p_tag->end
        << " [" << p_tag->strand << "] l=" << (p_tag->end - p_tag->start + 1) << endl;
  }
}


/*!
    \fn Link::get_shortest_length()
 */
uint Link::get_shortest_length()
{
  list<tag>::iterator p_tag = this->tags.begin();
  uint shortest_length = p_tag->end - p_tag->start + 1;
  for (p_tag++; p_tag != this->tags.end(); p_tag++) {
    uint length = p_tag->end - p_tag->start + 1;
    if (length < shortest_length) {
      shortest_length = length;
    }
  }
  return shortest_length;
}


/*!
    \fn Link::merge(Link* other_link)
 */
Link* Link::merge(Link* other_link)
{
  if (this->anchor_list.front() == other_link->anchor_list.back() and
      this->anchor_list.back() == other_link->anchor_list.front()) {
    other_link->reverse();
  }
  list<Anchor*>::iterator p_anchor_1 = this->anchor_list.begin();
  list<Anchor*>::iterator p_anchor_2 = other_link->anchor_list.begin();
//   cout << "Merging:" <<endl;
//   this->print();
//   other_link->print();
//   string kk;
//   cin >> kk;

  while (p_anchor_1 != this->anchor_list.end() and p_anchor_2 != other_link->anchor_list.end()) {
    if (*p_anchor_1 == *p_anchor_2) {
      p_anchor_1++;
      p_anchor_2++;
    } else {
      int dist = 1;
      bool match = false;
      list<Anchor*>::iterator p_this_anchor = p_anchor_1;
      p_this_anchor++;
      while (p_this_anchor != this->anchor_list.end()) {
        if (*p_this_anchor == *p_anchor_2) {
          match = true;
          break;
        }
        dist++;
        p_this_anchor++;
      }
      if (match) {
        p_anchor_1 = p_this_anchor;
        continue;
      }
      dist = 1;
      p_this_anchor = p_anchor_2;
      p_this_anchor++;
      while (p_this_anchor != other_link->anchor_list.end()) {
        if (*p_anchor_1 == *p_this_anchor) {
          match = true;
          break;
        }
        dist++;
        p_this_anchor++;
      }
      if (match) {
        for (int a=0; a<dist; a++) {
          this->anchor_list.insert(p_anchor_1, *p_anchor_2);
          p_anchor_2++;
        }
        continue;
      }
      p_anchor_1++;
      this->anchor_list.insert(p_anchor_1, *p_anchor_2);
      p_anchor_2++;
    }
  }
  this->tags.splice(this->tags.end(), other_link->tags);
//   cout << " --> this:" <<endl;
//   this->print();
//   other_link->print();
  for (list<Link*>::iterator p_link_it = other_link->anchor_list.front()->links.begin(); p_link_it != other_link->anchor_list.front()->links.end(); p_link_it++) {
    if (*p_link_it == other_link) {
      other_link->anchor_list.front()->links.erase(p_link_it);
      break;
    }
  }
  for (list<Link*>::iterator p_link_it = other_link->anchor_list.back()->links.begin(); p_link_it != other_link->anchor_list.back()->links.end(); p_link_it++) {
    if (*p_link_it == other_link) {
      other_link->anchor_list.back()->links.erase(p_link_it);
      break;
    }
  }
  delete(other_link);
  return this;
}


/*!
    \fn Link::is_an_alternative_path_of(Link* other_link)
 */
bool Link::is_an_alternative_path_of(Link* other_link)
{
  if (this->anchor_list.front() == other_link->anchor_list.front() and
      this->anchor_list.back() == other_link->anchor_list.back()) {
    return true;
  }
  if (this->anchor_list.front() == other_link->anchor_list.back() and
      this->anchor_list.back() == other_link->anchor_list.front()) {
    return true;
  }
  return false;
}


/*!
    \fn Link::get_num_of_mismatches(Link* other_link)
 */
uint Link::get_num_of_mismatches(Link* other_link)
{
  int distance = 0;
  if (this->anchor_list.front() == other_link->anchor_list.back() and
      this->anchor_list.back() == other_link->anchor_list.front()) {
    other_link->reverse();
  }
  list<Anchor*>::iterator p_anchor_1 = this->anchor_list.begin();
  list<Anchor*>::iterator p_anchor_2 = other_link->anchor_list.begin();

  while (p_anchor_1 != this->anchor_list.end() and p_anchor_2 != other_link->anchor_list.end()) {
    if (*p_anchor_1 == *p_anchor_2) {
      p_anchor_1++;
      p_anchor_2++;
    } else {
      int dist = 1;
      bool match = false;
      list<Anchor*>::iterator p_this_anchor = p_anchor_1;
      p_this_anchor++;
      while (p_this_anchor != this->anchor_list.end()) {
        if (*p_this_anchor == *p_anchor_2) {
          match = true;
          break;
        }
        dist++;
        p_this_anchor++;
      }
      if (match) {
        distance += dist;
        p_anchor_1 = p_this_anchor;
        continue;
      }
      dist = 1;
      p_this_anchor = p_anchor_2;
      p_this_anchor++;
      while (p_this_anchor != other_link->anchor_list.end()) {
        if (*p_anchor_1 == *p_this_anchor) {
          match = true;
          break;
        }
        dist++;
        p_this_anchor++;
      }
      if (match) {
        distance += dist;
        p_anchor_2 = p_this_anchor;
        continue;
      }
      p_anchor_1++;
      p_anchor_2++;
      distance+=2;
    }
  }
//   cout << distance << endl;
  return distance;
}
