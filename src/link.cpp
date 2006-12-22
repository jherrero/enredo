#include <iostream>
#include "link.h"
#include "anchor.h"
#include <vector>
#include <map>
#include <cstdlib>

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
    \fn Link::try_to_concatenate_with(Link *other_link, short strand1, short strand2)
 */
bool Link::try_to_concatenate_with(Link *other_link, short strand1, short strand2)
{
  if (strand1 == 0) {
    if (this->try_to_concatenate_with(other_link, 1, strand2)) {
      return true;
    } else {
      return this->try_to_concatenate_with(other_link, -1, strand2);
    }
  }
  if (strand2 == 0) {
    if (this->try_to_concatenate_with(other_link, strand1, 1)) {
      return true;
    } else {
      return this->try_to_concatenate_with(other_link, strand1, -1);
    }
  }

  std::vector< std::list<tag>::iterator > this_tag_links_to(this->tags.size(), other_link->tags.end());
  std::vector< std::list<tag>::iterator > other_tag_links_to(this->tags.size(), this->tags.end());
  uint this_tag_counter = 0;
  for (list<tag>::iterator p_tag1 = this->tags.begin(); p_tag1 != this->tags.end(); p_tag1++) {
    this_tag_counter++;
    short str1 = strand1 * p_tag1->strand;
    uint other_tag_counter = 0;
    for (list<tag>::iterator p_tag2 = other_link->tags.begin(); p_tag2 != other_link->tags.end(); p_tag2++) {
      other_tag_counter++;
      if (p_tag1->species == p_tag2->species and p_tag1->chr == p_tag2->chr and p_tag1->start < p_tag2->end and p_tag2->start < p_tag1->end) {
        short str2 = strand2 * p_tag2->strand;
        if (str1 == 1 and str2 == 1) {
          if (!(p_tag1->start < p_tag2->start and p_tag1->end < p_tag2->end)) {
            // link goes and come back: they will be concatenated when studying the other anchor
            continue;
          }
        } else if (str1 == -1 and str2 == -1) {
          if (!(p_tag2->start < p_tag1->start and p_tag2->end < p_tag1->end)) {
            // link goes and come back: they will be concatenated when studying the other anchor
            continue;
          }
        } else if (str1 != 0 and str2 != 0) {
          cerr << "Problem here (3) " << strand1 << " : " << strand2 << endl;
          this->print();
          other_link->print();
          continue;
        }
        if (other_tag_links_to[other_tag_counter - 1] == 0) {
          continue;
        }
        other_tag_links_to[other_tag_counter - 1] = p_tag1;
        this_tag_links_to[this_tag_counter - 1] = p_tag2;
        continue;
      }
    }
    if (this_tag_links_to[this_tag_counter - 1] == other_link->tags.end()) {
      return false;
    }
  }
  for (uint a = 0; a < other_tag_links_to.size(); a++) {
    if (other_tag_links_to[a] == this->tags.end()) {
      return false;
    }
  }

  if (strand1 == -1) {
    this->reverse();
  }
  if (strand2 == -1) {
    other_link->reverse();
  }

  // Concatenate the links
  this_tag_counter = 0;
  for (list<tag>::iterator p_tag1 = this->tags.begin(); p_tag1 != this->tags.end(); p_tag1++) {
    std::list<tag>::iterator p_tag2 = this_tag_links_to[this_tag_counter];
//     cout << "concatenating (" << *p_tag1->species << ":" << *p_tag1->chr << ":" << p_tag1->start << ":"
//         << p_tag1->end << ":" << p_tag1->strand << ")"
//         << " -- (" << *p_tag2->species << ":" << *p_tag2->chr << ":" << p_tag2->start << ":"
//         << p_tag2->end << ":" << p_tag2->strand << ")";
    if (p_tag1->strand != 0 and p_tag1->strand == -p_tag2->strand) {
      cerr << "Problem here (5)" << endl;
    } else if (p_tag1->strand == 1 or p_tag2->strand == 1) {
      p_tag1->end = p_tag2->end;
    } else if (p_tag1->strand == -1 or p_tag2->strand == -1) {
      p_tag1->start = p_tag2->start;
    }
    if (p_tag1->strand == 0 ) {
      p_tag1->strand = p_tag2->strand;
    }
//     cout << " => (" << *p_tag1->species << ":" << *p_tag1->chr << ":" << p_tag1->start << ":"
//         << p_tag1->end << ":" << p_tag1->strand << ")" << endl;
    this_tag_counter++;
  }

  Anchor *front_anchor = this->anchor_list.front();
  Anchor *middle_anchor = this->anchor_list.back();
  Anchor *back_anchor = other_link->anchor_list.back();
  if (middle_anchor != other_link->anchor_list.front()) {
    cerr << "Problem here (6)" << endl;
    exit(6);
  }

  // Delete other link from middle and back anchors
  for (list<Link*>::iterator p_link_it = middle_anchor->links.begin(); p_link_it != middle_anchor->links.end(); p_link_it++) {
    if (*p_link_it == other_link) {
      p_link_it = middle_anchor->links.erase(p_link_it);
    }
  }
  for (list<Link*>::iterator p_link_it = back_anchor->links.begin(); p_link_it != back_anchor->links.end(); p_link_it++) {
    if (*p_link_it == other_link) {
      p_link_it = back_anchor->links.erase(p_link_it);
    }
  }

  // Delete this link from middle anchor if needed
  if (front_anchor != middle_anchor) {
    // Remove this link from middle anchor
    for (list<Link*>::iterator p_link_it = middle_anchor->links.begin(); p_link_it != middle_anchor->links.end(); p_link_it++) {
      while (*p_link_it == this or *p_link_it == other_link) {
        p_link_it = middle_anchor->links.erase(p_link_it);
      }
    }
  }
  if (front_anchor != back_anchor) {
    // Add this link to back anchor
    back_anchor->add_Link(this);
  }

  // Append anchors from other_link to this link
  for (list<Anchor*>::iterator anchor_it = ++other_link->anchor_list.begin(); anchor_it != other_link->anchor_list.end(); anchor_it++) {
    this->anchor_list.push_back(*anchor_it);
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
    out << "       ";
    print_tag(*p_tag, out);
    out << endl;
  }
}


/*!
    \fn Link::print_tag(tag this_tag, ostream &out)
 */
void print_tag(tag this_tag, ostream &out)
{
  out << *this_tag.species << ":" << *this_tag.chr << ":"
      << this_tag.start << ":" << this_tag.end
      << " [" << this_tag.strand << "] l=" << (this_tag.end - this_tag.start + 1);
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
