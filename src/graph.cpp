#include "graph.h"
#include "anchor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <set>
#include <vector>
#include <fstream>

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
bool Graph::populate_from_file(char *filename, float min_score, int max_gap_length)
{
  ifstream inputfile (filename);
  string line;
  if (!inputfile.is_open()) {
    cerr << "Cannot open file " << filename << endl;
    return false;
  }

  unsigned long long int line_counter = 0;
  uint long_gap_counter = 0;
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
    string this_anchor_id;
    string this_species;
    string this_chr;
    int this_start;
    int this_end;
    string this_strand;
    float this_score;

    stringstream streamline;
    streamline << line;
    streamline
        >> this_anchor_id >> this_species >> this_chr
        >> this_start >> this_end >> this_strand >> this_score;
    if (streamline.fail()) {
      cerr << "Error reading line (" << line_counter << ")<" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    if (this_start > this_end) {
      cerr << "start cannot be longer than end in <" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    if (this_strand != "+" and this_strand != "-") {
      cerr << "strand must be + or - in <" << line << ">" << endl;
      inputfile.close();
      return false;
    }
    if (this_score < min_score) {
      continue;
    }
    Anchor *anchor = this->get_Anchor(this_anchor_id);
    if (!anchor) {
      cerr << "Out of memory" << endl;
      return false;
    }
    if (!species[this_species]) {
      cout << "New species " << this_species << endl;
      species[this_species] = new string(this_species);
    }
    anchor->species.insert(species[this_species]);
    if (last_species == this_species and
        last_chr == this_chr and
        last_end < this_start) {
      if ((max_gap_length > 0) and (this_start - last_end - 1 > max_gap_length)) {
        long_gap_counter++;
      } else {
        Link *this_link = anchor->get_direct_Link(last_anchor);
        short this_link_strand;
        if (last_anchor == anchor) {
          this_link_strand = 0;
        } else if (last_anchor == *this_link->anchor_list.begin()) {
          this_link_strand = 1;
        } else if (anchor == *this_link->anchor_list.begin()) {
          this_link_strand = -1;
        } else {
          cerr << "Error";
          exit(1);
        }
        string *this_chr = chrs[last_chr];
        if (!this_chr) {
          this_chr = new string(last_chr);
          chrs[last_chr] = this_chr;
        }
        this_link->add_tag(species[this_species], this_chr, last_start, this_end, this_link_strand);
      }
    }
    last_anchor = anchor;
    last_species = this_species;
    last_chr = this_chr;
    last_start = this_start;
    last_end = this_end;
    last_strand = this_strand;

    line_counter++;
//     if (!(line_counter % 10000)) {
//       cerr << "Line " << line_counter << endl;
//     }
  }
  inputfile.close();
  cout << "Number of long gaps (larger than " << max_gap_length << "): " << long_gap_counter << endl;
  return true;
}




/*!
    \fn Graph::minimize()
*/
void Graph::minimize()
{
  int count = 0;
//   anchors["10_11557"]->print();
//   anchors["9_12874"]->print();
  cout << "Minimizing graph..." << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor *this_anchor = it->second;
//     if (it->second->id == "10_11557" or it->second->id == "9_12874") {
//       cout << "=================== ANCHOR " << it->second->id << " ===========================" << endl;
//       it->second->print();
//     } else {
//       continue;
//     }
    bool merge_event;
    do {
      merge_event = false;
      for (std::list<Link*>::iterator p_link1 = this_anchor->links.begin(); p_link1 != this_anchor->links.end() and !merge_event; p_link1++) {
        short strand1;
        if ((*p_link1)->anchor_list.front() == (*p_link1)->anchor_list.back()) {
          strand1 = 0;
        } else if ((*p_link1)->anchor_list.back() == this_anchor) {
          strand1 = 1;
        } else if ((*p_link1)->anchor_list.front() == this_anchor) {
          strand1 = -1;
        }
        for (std::list<Link*>::iterator p_link2 = this_anchor->links.begin(); p_link2 != p_link1 and !merge_event; p_link2++) {
          if ((*p_link2)->tags.size() != (*p_link1)->tags.size()) {
            continue;
          }
          short strand2;
          if ((*p_link2)->anchor_list.front() == (*p_link2)->anchor_list.back()) {
            strand2 = 0;
          } else if ((*p_link2)->anchor_list.back() == this_anchor) {
              strand2 = -1;
          } else if ((*p_link2)->anchor_list.front() == this_anchor) {
            strand2 = 1;
          }
          if ((*p_link1)->try_to_concatenate_with(*p_link2, strand1, strand2)) {
            merge_event = true;
            count++;
          }
        }
      }
    } while (merge_event);
  }
  cout << count << " merges." << endl;
  
}


/*!
    \fn Graph::print_anchors_histogram(std::ostream &out)
 */
void Graph::print_anchors_histogram(std::ostream &out)
{
  std::vector<unsigned long long int> hist_hits;
  std::map< std::string, unsigned long long int> hist_anchors_per_species;
  std::vector<unsigned long long int> hist_species(species.size(), 0);
  std::map< std::string, unsigned long long int> hist_patterns;

  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor *this_anchor = it->second;
    uint num = this_anchor->num;
    if (num > hist_hits.size()) {
      hist_hits.resize(num, 0);
    }
    hist_hits[num - 1]++;
    if (this_anchor->species.size() > 0 ) {
      hist_species[this_anchor->species.size() - 1]++;
      std::string pattern = "";
      for (std::map<std::string, std::string*>::iterator it2 = species.begin(); it2 != species.end(); it2++) {
        for (std::set<std::string*>::iterator it3 = this_anchor->species.begin(); it3 != this_anchor->species.end(); it3++) {
          if (it2->second == *it3) {
            hist_anchors_per_species[**it3]++;
            if (pattern != "") {
              pattern.append(" + ");
            }
            pattern.append(**it3);
          }
        }
      }
//       if (!hist_patterns[pattern]) {
//         hist_patterns[pattern] = 1;
//       } else {
        hist_patterns[pattern]++;
//       }
    } else {
      cerr << "Hmmm... There is an anchor (" << this_anchor->id << ") mapping on no species! :-S" << endl;
    }
  }
  out << endl << "Histogram of num. of Anchors per species" << endl;
  for (std::map< std::string, unsigned long long int>::iterator it = hist_anchors_per_species.begin(); it != hist_anchors_per_species.end(); it++) {
    cout << it->first << ": " << it->second << " (" << (100.0f * it->second) / anchors.size() << "%)" << endl;
  }
  out << endl << "Histogram of num. of species per Anchor (in how many species each Anchor is found)" << endl;
  uint sum = 0;
  ios::fmtflags current_flags = out.flags();
  out.setf(ios::fixed);
  out.precision(1);
  for (uint a = 0; a < hist_species.size(); a++) {
    sum += hist_species[a];
  }
  for (uint a = 0; a < hist_species.size(); a++) {
    out << a + 1 << ": " << hist_species[a] << " (" << (100.0f * hist_species[a]) / sum << "%)" << endl;
  }
  out << endl << "The same, by set of species (in no particular order)" << endl;
  for (std::map< std::string, unsigned long long int>::iterator it = hist_patterns.begin(); it != hist_patterns.end(); it++) {
    cout << it->first << ": " << it->second << " (" << (100.0f * it->second) / sum << "%)" << endl;
  }
  out << endl << "Histogram of num. of hits per Anchor (how many times each Anchor is found)" << endl;
  sum = 0;
  for (uint a = 0; a < hist_hits.size(); a++) {
    sum += hist_hits[a];
  }
  for (uint a = 0; a < hist_hits.size(); a++) {
    if (hist_hits[a]) {
      out << a + 1 << ": " << hist_hits[a] << " (" << (100.0f * hist_hits[a]) / sum << "%)" << endl;
    }
  }
  out.flags(current_flags);
}


/*!
    \fn Graph::print_stats(int histogram_size)
 */
void Graph::print_stats(int histogram_size)
{
  unsigned long long int non_void_anchors_counter = 0;
  unsigned long long int links_counter = 0;
  std::vector<unsigned long long int> hist(histogram_size, 0);
  std::map< std::string, std::vector<unsigned long long int> > hist_per_species;
  std::map< std::string, list<uint> > lengths;
  std::map< std::string, unsigned long long int > total_length;
  std::vector< std::map< std::string, list<uint> > > lengths_per_cardinality(histogram_size, lengths);
  std::map< std::string, std::vector< uint > > length_per_species_cardinality;

  for (std::map<std::string, std::string*>::iterator it = species.begin(); it != species.end(); it++) {
    hist_per_species[it->first].resize(histogram_size, 0);
  }
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    if (this_anchor->links.size()) {
      non_void_anchors_counter++;
    }
    for (std::list<Link*>::iterator p_link = this_anchor->links.begin(); p_link != this_anchor->links.end(); p_link++) {
      if ((*p_link)->anchor_list.front()->id != this_anchor->id) {
        /* Each link will be counted twice, once when counting links for the
        first anchor and once when counting links for the last anchor
        (except when the fist and the last anchors are the same one). */
        continue;
      }
      links_counter++;

      int size = (*p_link)->tags.size();
      if (size > histogram_size) {
        size = histogram_size;
      }
      hist[size - 1]++;

      if (size > 0) {
        std::map< std::string, unsigned int > cardinality_per_species;
        for (std::list<tag>::iterator p_tag = (*p_link)->tags.begin(); p_tag != (*p_link)->tags.end(); p_tag++) {
          string this_species = *(p_tag->species);
          uint this_length = p_tag->end - p_tag->start + 1;
          hist_per_species[this_species][size - 1]++;
          lengths[this_species].push_back(this_length);
          lengths_per_cardinality[size - 1][this_species].push_back(this_length);
          cardinality_per_species[this_species]++;
        }
        for (std::list<tag>::iterator p_tag = (*p_link)->tags.begin(); p_tag != (*p_link)->tags.end(); p_tag++) {
          string this_species = *(p_tag->species);
          uint this_length = p_tag->end - p_tag->start + 1;
          if (length_per_species_cardinality[this_species].size() < cardinality_per_species[this_species]) {
            length_per_species_cardinality[this_species].resize(cardinality_per_species[this_species], 0);
          }
          length_per_species_cardinality[this_species][cardinality_per_species[this_species]-1] += this_length;
        }
      }

    }
  }
  ios::fmtflags current_flags = cout.flags();
  cout.setf(ios::fixed);
  cout.precision(1);
  cout << "Graph has " << non_void_anchors_counter << " non-void anchors ("
      << anchors.size() << " in total) and " << links_counter << " links (edges)" << endl << endl;

  cout << "Duplications according to graph (length in bp)" << endl;
  cout << "|! species\t|! 1x\t|! 2x\t|! 3x\t|! 4x\t|! 5x\t|" << endl;
  for (std::map<std::string, std::string*>::iterator it = species.begin(); it != species.end(); it++) {
    cout << "|! " << it->first;
    for (uint a = 0; a < 5; a++) {
      if (a < length_per_species_cardinality[it->first].size()) {
        cout << "\t| " << length_per_species_cardinality[it->first][a];
      } else {
        cout << "\t| 0";
      }
    }
    cout << "\t|" << endl;
  }
  cout << endl;

  for (std::map<std::string, list<uint> >::iterator it = lengths.begin(); it != lengths.end(); it++) {
    cout << "N50 for " << it->first << ": ";
    list<uint> *this_list = &(it->second);
    this_list->sort();
    unsigned long long int sum = 0;
    for (std::list<uint>::iterator p_length = this_list->begin(); p_length != this_list->end(); p_length++) {
      sum += *p_length;
    }
    unsigned long long int acc = 0;
    for (std::list<uint>::iterator p_length = this_list->begin(); p_length != this_list->end(); p_length++) {
      acc += *p_length;
      if (acc * 2 > sum) {
        cout << *p_length;
        break;
      }
    }
    total_length[it->first] = sum;
    cout << " (total length = " << sum << ")" << endl;
  }
  cout << endl;

  cout << "Detailed N50 stats per species and cardinality" << endl;
  for (std::map<std::string, list<uint> >::iterator it = lengths.begin(); it != lengths.end(); it++) {
    cout << "|>|>|>|>|>|! " << it->first << " |" << endl;
    cout << "|! Link cardinality |! Total num. of Links |! num of links |! Total Length |! length% |! N50 |" << endl;
    for (int a=0; a < histogram_size; a++) {
      if (a == histogram_size - 1) {
        cout << "|! >" << a << " | " << hist[a];
      } else {
        cout << "|! " << a + 1 << " | " << hist[a];
      }

      cout << " | " << hist_per_species[it->first][a] << " | ";
      list<uint> *this_list = &lengths_per_cardinality[a][it->first];
      this_list->sort();
      unsigned long long int sum = 0;
      for (std::list<uint>::iterator p_length = this_list->begin(); p_length != this_list->end(); p_length++) {
        sum += *p_length;
      }
      cout << sum << " | " << 100.0f * sum / total_length[it->first] << "% | ";
      unsigned long long int acc = 0;
      for (std::list<uint>::iterator p_length = this_list->begin(); p_length != this_list->end(); p_length++) {
        acc += *p_length;
        if (acc * 2 > sum) {
          cout << *p_length << " |";
          break;
        }
      }
      cout << endl;
    }
    cout << endl;

  }
  cout.flags(current_flags);
}


/*!
    \fn Graph::print_links(ostream &out, int min_anchors, int min_regions, int min_length)
 */
unsigned long int Graph::print_links(ostream &out, uint min_anchors, uint min_regions, uint min_length)
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
    (*p_link_it)->print(out);
  }
  return num_blocks;
}


/*!
    \fn Graph::merge_alternative_paths(int max_anchors)
 */
void Graph::merge_alternative_paths(uint max_anchors)
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
//           (*p_link1)->print();
//           (*p_link2)->print();
          if ((*p_link1)->is_an_alternative_path_of(*p_link2)) {
            if ((*p_link1)->get_num_of_mismatches(*p_link2) <= max_anchors) {
              count++;
              (*p_link1)->merge(*p_link2);
              merge_event = true;
            }
          }
        }
      }
    } while (merge_event);
//     cout << endl;
  }
  cout << count << " merges." << endl;
}


/*!
    \fn Graph::simplify()
 */
void Graph::simplify(uint min_anchors, uint min_regions, uint min_length)
{
  cout << "Simplifying graph..." << endl;
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->anchor_list.size() >= min_anchors and this_link->tags.size() >= min_regions and this_link->get_shortest_length() >= min_length) {
        all_links.insert(this_link).second;
      }
    }
  }
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    //TODO: implement me! Look for links that won't be selected as syntenic regions and try to break them
    // in order to potentially enlarge neighbouring links
  }
}


/*!
    \fn Graph::study_anchors()
 */
void Graph::study_anchors(void)
{
  int count = 0;
  uint inversion_count = 0;
  uint end_path_count = 0;
  uint tandem_count = 0;
  uint pure_branching_count = 0;
  uint singular_branching_count = 0;
  uint pure_other_count = 0;
  uint singular_other_count = 0;
  bool debug = false;

  cout << "Studying anchors..." << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor *this_anchor = it->second;
//     if (this_anchor->id != "11_23281" and this_anchor->id != "1_9852") {
//       continue;
//     }
    
    if (this_anchor->links.size() == 0) {
      continue;
    }
    std::map< uint, uint > path;
    vector <tag> all_tags;
    vector <uint> link_nums;
    vector <uint> tag_nums;
    vector <short> strands;
    uint link_num = 1;
    for (std::list<Link*>::iterator p_link1 = this_anchor->links.begin(); p_link1 != this_anchor->links.end(); p_link1++) {
      uint tag_num = 1;
      short link_strand;
      if ((*p_link1)->anchor_list.front() == this_anchor) {
        if ((*p_link1)->anchor_list.back() == this_anchor) {
          link_strand = 0;
        } else {
          link_strand = 1;
        }
      } else if ((*p_link1)->anchor_list.back() == this_anchor) {
        link_strand = -1;
      } else {
        cerr << endl << "ERROR #12 (anchor=" << this_anchor->id << ")" << endl;
        (*p_link1)->print(cerr);
        exit(12);
      }
      for (std::list<tag>::iterator p_tag1 = (*p_link1)->tags.begin(); p_tag1 != (*p_link1)->tags.end(); p_tag1++) {
        all_tags.push_back(*p_tag1);
        link_nums.push_back(link_num);
        tag_nums.push_back(tag_num);
        strands.push_back(link_strand);
        tag_num++;
      }
      link_num++;
    }
    if (debug) {
      this_anchor->print();
      cout << endl;
      for (uint c = 0; c < all_tags.size(); c++) {
        cout << "  " << strands[c] << " (" << link_nums[c] << "." << tag_nums[c] << ") ";
        print_tag(all_tags[c]);
        cout << endl;
      }
    }
    vector <uint> match (all_tags.size(), 0);
    vector <uint> inv_match (all_tags.size(), 0);
    for (uint a = 0; a < all_tags.size()-1; a++) {
      for (uint b = a + 1; b < all_tags.size(); b++) {
        tag *p_tag1 = &all_tags[a];
        tag *p_tag2 = &all_tags[b];
        if (p_tag1->species != p_tag2->species or p_tag1->chr != p_tag2->chr) {
          continue;
        }
        short strand1 = p_tag1->strand * strands[a];
        short strand2 = p_tag2->strand * strands[b];
        if (p_tag1->start < p_tag2->end and p_tag2->start < p_tag1->end) {
          if (link_nums[a] == link_nums[b]) {
            if (strand1 == 1 and strand2 == -1) {
              if (p_tag2->start < p_tag1->start and p_tag2->end < p_tag1->end) {
//                 cout << "INVERSION: ";
                inv_match[a]++;
                inv_match[b]++;
              }
            } else if (strand1 == -1 and strand2 == 1) {
              if (p_tag1->start < p_tag2->start and p_tag1->end < p_tag2->end) {
//                 cout << "INVERSION: ";
                inv_match[a]++;
                inv_match[b]++;
              }
            } else if (strand1 == 0 and strand2 == 0) {
//               cout << "Loop: ";
              match[a]++;
              match[b]++;
            } else {
              cerr << endl << "ERROR #13" << endl;
              exit(1);
            }
          } else {
            if (strand1 == 1 and strand2 == -1) {
              if (p_tag2->start < p_tag1->start and p_tag2->end < p_tag1->end) {
//                 cout << "Match: ";
                match[a]++;
                match[b]++;
                path[a] = b;
                path[b] = a;
              }
            } else if (strand1 == -1 and strand2 == 1) {
              if (p_tag1->start < p_tag2->start and p_tag1->end < p_tag2->end) {
//                 cout << "Match: ";
                match[a]++;
                match[b]++;
                path[a] = b;
                path[b] = a;
              }
            } else {
              // At least one of the tags is for a loop
              match[a]++;
              match[b]++;
              path[a] = b;
              path[b] = a;
            }
          }
          if (debug) {
                        cout << "(" << link_nums[a] << "." << tag_nums[a] << ") ";
                        print_tag(*p_tag1);
                        cout << " -- (" << link_nums[b] << "." << tag_nums[b] << ") ";
                        print_tag(*p_tag2);
                        cout << endl;
          }
        }
      }
    }
    uint this_inversion_count = 0;
    uint this_end_path_count = 0;
    uint this_tandem_count = 0;
    vector <bool> left(match.size(), false);
    for (uint a = 0; a < match.size(); a++) {
      if (match[a] == 0) {
        if (inv_match[a] == 1) {
//           cout << "[inversion] (" << link_nums[a] << "." << tag_nums[a] << ")";
          this_inversion_count++;
        } else if (inv_match[a] > 1) {
          this_anchor->print();
          cout << "[ARGHHH] (" << link_nums[a] << "." << tag_nums[a] << ")";
          string resp;
          cin >> resp;
        } else {
//           cout << "[end_path] (" << link_nums[a] << "." << tag_nums[a] << ")";
          this_end_path_count++;
        }
      } else if (match[a] > 1) {
//         cout << "[tandem] (" << link_nums[a] << "." << tag_nums[a] << ")";
        this_tandem_count++;
      } else if (strands[a] == 0) {
//         cout << "[end_path] (" << link_nums[a] << "." << tag_nums[a] << ")";
        this_end_path_count++;
      } else {
        left[a] = true;
//         other_count++;
      }
    }
    inversion_count += this_inversion_count;
    end_path_count += this_end_path_count;
    tandem_count += this_tandem_count;

    // look for branching point
    bool branching;
    for (uint a = 1; a < this_anchor->links.size() + 1; a++) {
      branching = true;
      for (uint b = 0; b < match.size(); b++) {
        if (!left[b]) continue;
        if (link_nums[b] == a) continue;
        if (!path.count(b) or link_nums[path[b]] != a) {
          branching = false;
          continue;
        }
      }
      if (branching) {
        for (uint b = 0; b < match.size(); b++) {
          if (!left[b]) continue;
          if (link_nums[b] != a) continue;
//           cout << "[bifurcation " << a << "] (" << link_nums[b] << "." << tag_nums[b] << ") => ("
//               << link_nums[path[b]] << "." << tag_nums[path[b]] << ")" << endl;
          left[b] = false;
          left[path[b]] = false;
        }
        continue;
      }
    }
    if (branching) {
      if (this_inversion_count or this_end_path_count or this_tandem_count) {
        singular_branching_count++;
      } else {
        pure_branching_count++;
      }
    } else {
      if (this_inversion_count or this_end_path_count or this_tandem_count) {
        singular_other_count++;
      } else {
        pure_other_count++;
      }
      for (uint b = 0; b < match.size(); b++) {
        if (!left[b]) continue;
//       cout << "[left] (" << link_nums[b] << "." << tag_nums[b] << ")";
//       if (path.count(b)) {
//         cout << " => (" << link_nums[path[b]] << "." << tag_nums[path[b]] << ")";
//         if (!left[path[b]]) {
//           if (match[path[b]] > 1) {
//             cout << " --tandem--";
//           } else {
//             cout << " -- other--";
//           }
//         }
//       }
//       cout << endl;
//       this_anchor->print();
//       cout << endl;
//       for (uint c = 0; c < all_tags.size(); c++) {
//         cout << "  " << strands[c] << " (" << link_nums[c] << "." << tag_nums[c] << ") ";
//         print_tag(all_tags[c]);
//         cout << endl;
//       }

//       cout << "[Enter] ";
//       string resp;
//       cin >> resp;
      }
    }

  }
  cout << end_path_count << " end of path." << endl;
  cout << (inversion_count/2) << " inversions." << endl;
  cout << tandem_count << " tandem anchors." << endl;
  cout << pure_branching_count << " bifurcating anchors." << endl;
  cout << singular_branching_count << " bifurcating anchors + singularity." << endl;
  cout << pure_other_count << " other situations." << endl;
  cout << singular_other_count << " other situations + singularity." << endl;
  cout << count << " merges." << endl;
}
