#include <iostream>
#include "link.h"
#include "anchor.h"
#include <vector>
#include <map>
#include <cstdlib>
#include <iomanip>

Link::Link(Anchor* anchor1, Anchor* anchor2)
{
  this->anchor_list.push_back(anchor1);
  this->anchor_list.push_back(anchor2);
//   cerr << "New link" << anchor1->id << ":" << anchor2->id << endl;
}

/*!
    \fn Link::Link(Link my_link)
 */
Link::Link(Link *my_link)
{
  
  for (list<Anchor*>::iterator p_anchor_it = my_link->anchor_list.begin(); p_anchor_it != my_link->anchor_list.end(); p_anchor_it++) {
    this->anchor_list.push_back(*p_anchor_it);
  }
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
    \fn Link::get_matching_tags(Link *other_link, short strand1, short strand2, bool allow_partial_match)
 */
std::vector< std::list<tag>::iterator > Link::get_matching_tags(Link *other_link, short strand1, short strand2, bool allow_partial_match)
{
  std::vector< std::list<tag>::iterator > this_tag_links_to(this->tags.size(), other_link->tags.end());
  std::vector< std::list<tag>::iterator > other_tag_links_to(other_link->tags.size(), this->tags.end());
  if (strand1 == 0) {
    this_tag_links_to = this->get_matching_tags(other_link, 1, strand2);
    if (!this_tag_links_to.empty()) {
      return this_tag_links_to;
    } else {
      return this->get_matching_tags(other_link, -1, strand2);
    }
  }
  if (strand2 == 0) {
    this_tag_links_to = this->get_matching_tags(other_link, strand1, 1);
    if (!this_tag_links_to.empty()) {
      return this_tag_links_to;
    } else {
      return this->get_matching_tags(other_link, strand1, -1);
    }
  }

  uint this_tag_counter = 0;
  // p_tag1: iterator for tags in this link
  for (list<tag>::iterator p_tag1 = this->tags.begin(); p_tag1 != this->tags.end(); p_tag1++) {
    this_tag_counter++;
    short str1 = strand1 * p_tag1->strand;
    uint other_tag_counter = 0;
    // p_tag2: iterator for tags in the other link
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
          // testing wrong strand when 2 strands can be tested. return void vector
          this_tag_links_to.clear();
          return this_tag_links_to;
        }
        if (other_tag_links_to[other_tag_counter - 1] != this->tags.end()) {
          // This other tag already matches a tag. Skip this.
          continue;
        }
        other_tag_links_to[other_tag_counter - 1] = p_tag1;
        this_tag_links_to[this_tag_counter - 1] = p_tag2;
      }
    }
  }

  for (uint a = 0; a < other_tag_links_to.size(); a++) {
    // Check if any of the other tags has not been linked
    if (other_tag_links_to[a] == this->tags.end()) {
      if (allow_partial_match) continue;
      this_tag_links_to.clear();
      return this_tag_links_to;
    }
    // Check if any two of the other tags link to the same tag of this link
    for (uint b = a + 1; b < other_tag_links_to.size(); b++) {
      if (other_tag_links_to[a] == other_tag_links_to[b]) {
        this_tag_links_to.clear();
        return this_tag_links_to;
      }
    }
  }
  return this_tag_links_to;
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

  std::vector< std::list<tag>::iterator > this_tag_links_to = this->get_matching_tags(other_link, strand1, strand2);
  if (this_tag_links_to.empty()) {
    // vector will be empty if any of the tags in the other_link has no match in this link
    return false;
  } else {
    // check that all tags in this link have a match in the other_link
    for (uint a = 0; a < this_tag_links_to.size(); a++) {
      if (this_tag_links_to[a] == other_link->tags.end()) {
        cerr << "ERROR. One of the tags links to nothing (should have been detected before)";
        return false;
      }
    }
  } 

  if (strand1 == -1) {
    this->reverse();
  }
  if (strand2 == -1) {
    other_link->reverse();
  }

  // Concatenate the links
  uint this_tag_counter = 0;
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
  out << "block";
  for (std::list<Anchor*>::iterator p_anchor_it = this->anchor_list.begin(); p_anchor_it != this->anchor_list.end(); p_anchor_it++) {
    out << " - " << (*p_anchor_it)->id;
  }
  out << "  (made of " << this->tags.size() << " genomic regions)" << endl;
  for (list<tag>::iterator p_tag = this->tags.begin(); p_tag != this->tags.end(); p_tag++) {
    print_tag(*p_tag, out);
    out << endl;
  }
  out << endl;
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


/*!
    \fn Link::split(vector<bool> tags_to_split)
 */
Link* Link::split(vector<bool> tags_to_split)
{
  // Check that bool vector size matches the num of tags in this link
  if (tags_to_split.size() != this->tags.size()) {
    cerr << "Trying to split a link with " << this->tags.size() << " tags using a bool vector"
        << " containing only " << tags_to_split.size() << " values." << endl;
    return NULL;
  }

  Link* new_link = new Link(this);

  std::list<tag>::iterator p_tag_it = this->tags.begin();
  list<tag> tmp_tags;
  for (uint i=0; i < tags_to_split.size(); i++) {
    if (tags_to_split[i]) {
      new_link->tags.push_back(*p_tag_it);
    } else {
      tmp_tags.push_back(*p_tag_it);
    }
    p_tag_it++;
  }

  if (tmp_tags.size() == 0) {
    cerr << "Leaving empty link" << endl;
    exit(1);
  }
  this->tags = tmp_tags;
  new_link->anchor_list.front()->add_Link(new_link);
  if (new_link->anchor_list.front() != new_link->anchor_list.back()) {
    new_link->anchor_list.back()->add_Link(new_link);
  }

  return new_link;
}


/*!
    \fn Link::split(std::vector< std::list<tag>::iterator > tags_to_split)
 */
Link* Link::split(std::vector< std::list<tag>::iterator > tags_to_split)
{
  Link* new_link = new Link(this);

  std::list<tag>::iterator p_tag_it = this->tags.begin();
  list<tag> tmp_tags;
  for (uint i=0; i < this->tags.size(); i++) {
    bool tag_found = false;
    for (uint j=0; j < tags_to_split.size(); j++) {
      if (p_tag_it == tags_to_split[j]) {
        new_link->tags.push_back(*p_tag_it);
        tag_found = true;
        break;
      }
    }
    if (!tag_found) {
      tmp_tags.push_back(*p_tag_it);
    }
    p_tag_it++;
  }

  if (new_link->tags.size() != tags_to_split.size()) {
    cerr << "Couldn't assign all the tags" << endl;
    exit(1);
  }
  if (tmp_tags.size() == 0) {
    cerr << "Leaving empty link" << endl;
    exit(1);
  }
  this->tags = tmp_tags;
  new_link->anchor_list.front()->add_Link(new_link);
  if (new_link->anchor_list.front() != new_link->anchor_list.back()) {
    new_link->anchor_list.back()->add_Link(new_link);
  }

  return new_link;
}


/*!
    \fn Link::get_strand_for_matching_tags(Anchor* anchor)
 */
short Link::get_strand_for_matching_tags(Anchor* anchor)
{
  short strand;
  if (this->anchor_list.front() == this->anchor_list.back()) {
    // Cannot determine the strand. Use 0 which will try both strands
    strand = 0;
  } else if (this->anchor_list.back() == anchor) {
    strand = -1;
  } else if (this->anchor_list.front() == anchor) {
    strand = 1;
  }
  return strand;
}
