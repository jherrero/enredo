#include "graph.h"
#include "anchor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <set>

using namespace std;

Graph::Graph()
{
  anchors.clear();
}

Graph::~Graph()
{
}


/*!
    \fn Graph::add_Anchor(Anchor *this_anchor)
 */
void Graph::add_Anchor(Anchor *this_anchor)
{
  anchors[this_anchor->id] = this_anchor;
}


/*!
    \fn Graph::get_Anchor(string id)
 */
Anchor* Graph::get_Anchor(string id)
{
  if (anchors[id]) {
    anchors[id]->num ++;
    return anchors[id];
  }

  // id was not found => create a new Anchor
  Anchor* new_anchor = new Anchor(id);
  if (new_anchor) {
    add_Anchor(new_anchor);
    return new_anchor;
  }
  return NULL;
}


/*!
    \fn Graph::populate_from_file(string filename)
 */
bool Graph::populate_from_file(char *filename)
{
  ifstream inputfile (filename);
  string line;
  if (!inputfile.is_open()) {
    cerr << "Cannot open file " << filename << endl;
    return false;
  }

  int line_counter = 0;
  Anchor *last_anchor = NULL;
  string last_species;
  string last_chr;
  int last_start;
  int last_end;
  string last_strand;
  while (!inputfile.eof()) {
    getline(inputfile, line);
    if (line[0] == '#' or inputfile.eof()) {
      continue;
    }
    bool error = false;
    string anchor1_id;
    string species1;
    string chr1;
    int start1;
    int end1;
    string strand1;
    string anchor2_id;
    string species2;
    string chr2;
    int start2;
    int end2;
    string strand2;

    stringstream streamline;
    streamline << line;
    streamline
        >> anchor1_id >> species1 >> chr1 >> start1 >> end1 >> strand1;
    if (streamline.fail()) {
      cerr << "Error reading line (" << line_counter << ")<" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    if (start1 > end1) {
      cerr << "start cannot be larger than end in <" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    if (strand1 != "1" and strand1 != "0" and strand1 != "-1" and strand1 != "+" and strand1 != "-") {
      cerr << "strand must be 1 or -1 in <" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    Anchor *anchor1 = this->get_Anchor(anchor1_id);
    if (!anchor1) {
      cerr << "Out of memory" << endl;
      return false;
    }
    if (last_species == species1 and
        last_chr == chr1 and
        last_end < start1) {
      Link *this_link = anchor1->get_direct_Link(last_anchor);
      short this_strand;
      if (last_anchor == anchor1) {
        this_strand = 0;
      } else if (last_anchor == *this_link->anchor_list.begin()) {
        this_strand = 1;
      } else if (anchor1 == *this_link->anchor_list.begin()) {
        this_strand = -1;
      } else {
        cerr << "Error";
        exit(1);
      }
      string *this_species = species[last_species];
      if (!this_species) {
        cout << "New species " << last_species << endl;
        this_species = new string(last_species);
        species[last_species] = this_species;
      }
      string *this_chr = chrs[last_chr];
      if (!this_chr) {
        this_chr = new string(last_chr);
        chrs[last_chr] = this_chr;
      }
      this_link->add_tag(this_species, this_chr, last_start, end1, this_strand);
    }
    last_anchor = anchor1;
    last_species = species1;
    last_chr = chr1;
    last_start = start1;
    last_end = end1;
    last_strand = strand1;

    streamline
        >> anchor2_id >> species2 >> chr2 >> start2 >> end2 >> strand2;
    if (!streamline.fail()) {
      if (start2 > end2) {
        cerr << "start cannot be larger than end in <" << line << ">" << endl;
        inputfile.close();
        return false;
      }
      if (strand2 != "1" and strand2 != "0" and strand2 != "-1" and strand2 != "+" and strand2 != "-") {
        cerr << "strand must be 1 or -1 in <" << line << ">" << endl;
        inputfile.close();
        return false;
      }
      if (species1 != species2) {
        cerr << "species1 and species2 must match in <" << line << ">" << endl;
        inputfile.close();
        return false;
      }
      Anchor *anchor2 = this->get_Anchor(anchor2_id);
      if (anchor2 == NULL) {
        cerr << "Out of memory" << endl;
        return false;
      }
      if (last_species == species2 and
          last_chr == chr2 and
          last_end < start2) {
        Link *this_link = anchor2->get_direct_Link(last_anchor);
        short this_strand;
        if (last_anchor == anchor2) {
          this_strand = 0;
        } else if (last_anchor == *this_link->anchor_list.begin()) {
          this_strand = 1;
        } else if (anchor2 == *this_link->anchor_list.begin()) {
          this_strand = -1;
        } else {
          cerr << "Error";
          exit(1);
        }
        string *this_species = species[last_species];
        if (!this_species) {
          cout << "New species " << last_species << endl;
          this_species = new string(last_species);
          species[last_species] = this_species;
        }
        string *this_chr = chrs[last_chr];
        if (!this_chr) {
          this_chr = new string(last_chr);
          chrs[last_chr] = this_chr;
        }
        this_link->add_tag(this_species, this_chr, last_start, end2, this_strand);
      }
      last_anchor = anchor2;
      last_species = species2;
      last_chr = chr2;
      last_start = start2;
      last_end = end2;
      last_strand = strand2;
    }
    line_counter++;
    if (!(line_counter % 10000)) {
//       cerr << "Line " << line_counter << endl;
    }
  }
  inputfile.close();
  cout << "There are "<< anchors.size() << " Anchors in this graph." << endl;
  int hist[5] = {0,0,0,0,0};
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
//     cout << "Anchor " << it->second->id << endl;
    for (std::list<Link*>::iterator p_link = it->second->links.begin(); p_link != it->second->links.end(); p_link++) {
//       cout << "  " << (*p_link)->anchor1->id << " - " << (*p_link)->anchor2->id << endl;
//       for (list<tag>::iterator p_tag = (*p_link)->tags.begin(); p_tag != (*p_link)->tags.end(); p_tag++) {
//         cout << "       " << (*p_tag).species << ":" << (*p_tag).chr << ":" << (*p_tag).start << ":" << (*p_tag).end
//             << " [" << (*p_tag).strand << "]" << endl;
//       }
      int weight = 1;
      if ((*p_link)->anchor_list.front() == (*p_link)->anchor_list.back()) {
        weight = 2;
      }
      int size = (*p_link)->tags.size();
      if (size == 1) {
        hist[0]+=weight;
      } else if (size == 2) {
        hist[1]+=weight;
      } else if (size == 3) {
        hist[2]+=weight;
      } else if (size == 4) {
        hist[3]+=weight;
      } else {
        hist[4]+=weight;
      }
    }
  }
  cout << "Histogram of num. of regions per anchor" << endl;
  for (int a=0; a < 4; a++) {
    cout << a + 1 << ": " << hist[a]/2 << endl;
  }
  cout << "+: " << hist[4]/2 << endl;
  return true;
}




/*!
    \fn Graph::minimize()
*/
void Graph::minimize()
{
  int count = 0;
  cout << "Minimizing graph..." << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor *this_anchor = it->second;
//     if (it->second->id == "99959" or it->second->id == "289835" or it->second->id == "154794") {
//       cout << "=================== ANCHOR " << it->second->id << " ===========================" << endl;
//       it->second->print();
//     } else {
//       continue;
//     }
    bool merge_event;
    do {
      merge_event = false;
      for (std::list<Link*>::iterator p_link1 = this_anchor->links.begin(); p_link1 != this_anchor->links.end() and !merge_event; p_link1++) {
        for (std::list<Link*>::iterator p_link2 = this_anchor->links.begin(); p_link2 != p_link1 and !merge_event; p_link2++) {
          if ((*p_link2)->tags.size() != (*p_link1)->tags.size()) {
            continue;
          }
//           cout << " test ";
//           (*p_link1)->print();
//           cout << "  ?  ";
//           (*p_link2)->print();
          if ((*p_link1)->can_be_concatenated_with(*p_link2)) {
//             cout << " test ";
//             (*p_link1)->print();
//             cout << "  ?  ";
//             (*p_link2)->print();
//             cout << " YES " << endl;
            (*p_link1)->concatenate(*p_link2);
            merge_event = true;
            count++;
          } else {
//             cout << " NO" << endl;
          }
        }
      }
    } while (merge_event);
//     cout << endl;
  }
  cout << count << " merges." << endl;
  
}


/*!
    \fn Graph::print()
 */
void Graph::print()
{
  int count = 0;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
//     this_anchor->print();
    if (this_anchor->links.size()) {
      count++;
    }
  }
  cout << "Graph has " << count << " non-void anchors (" << anchors.size() << " in total)" << endl;
}


/*!
    \fn Graph::print_links(int min_anchors, int min_regions, int min_length)
 */
void Graph::print_links(int min_anchors, int min_regions, int min_length)
{
  int num_blocks = 0;
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->anchor_list.size() >= min_anchors and this_link->tags.size() >= min_regions and this_link->get_shortest_length() >= min_length) {
        if (all_links.insert(this_link).second) {
          num_blocks++;
        }
      }
    }
  }
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    (*p_link_it)->print();
  }
  cout << " Got " << all_links.size() << " blocks." << endl;
}


/*!
    \fn Graph::merge_alternative_paths(int max_anchors)
 */
void Graph::merge_alternative_paths(int max_anchors)
{
  int count = 0;
  cout << "Merging alternative paths..." << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    bool merge_event;
    Anchor *this_anchor = it->second;
    do {
      merge_event = false;
      for (std::list<Link*>::iterator p_link1 = this_anchor->links.begin(); p_link1 != this_anchor->links.end() and !merge_event; p_link1++) {
        for (std::list<Link*>::iterator p_link2 = this_anchor->links.begin(); p_link2 != p_link1 and !merge_event; p_link2++) {
          if ((*p_link1)->is_an_alternative_path_of(*p_link2)) {
            if ((*p_link1)->get_num_of_mismatches(*p_link2) <= max_anchors) {
              count++;
            (*p_link1)->merge(*p_link2);
            merge_event = true;
            }
//           } else {
//             cout << " NO" << endl;
          }
        }
      }
    } while (merge_event);
//     cout << endl;
  }
  cout << count << " merges." << endl;
}
