#include "graph.h"
#include "anchor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <set>
#include <vector>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

bool DEBUG = false;

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
bool Graph::populate_from_file(char *filename, float min_score, int max_gap_length, bool anchors_as_links)
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
    } else if (line == "--") {
      last_species = "";
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
        if (DEBUG) {
          cout << " ** LONG GAP **   " << this_species << ":" << this_chr << ":" << last_end << ".." <<
              this_start << "    " << last_anchor->id << " <--> " << anchor->id << endl;
        }
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
void Graph::minimize(std::string debug)
{
  uint count = 0;
//   anchors["10_11557"]->print();
//   anchors["9_12874"]->print();
  cout << "Minimizing graph..." << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor *this_anchor = it->second;
    uint num_merges;
    if (debug == "ALL" or (it->second->id == debug or it->second->id == debug)) {
      cout << "=================== ANCHOR " << it->second->id << " ===========================" << endl;
      it->second->print();
      num_merges = this_anchor->minimize(true);
    } else {
      num_merges = this_anchor->minimize(false);
    }
    if (num_merges > 0 and (debug == "ALL" or it->second->id == debug)) {
      cout << "------------ NEW ANCHOR AFTER MINIMIZATION " << it->second->id << " -----------" << endl;
      it->second->print();
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
    }
    count += num_merges;
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
    \fn Graph::print_links(ostream &out, int min_anchors, int min_regions, int min_length, bool allow_bridges)
 */
unsigned long int Graph::print_links(ostream &out, uint min_anchors, uint min_regions, uint min_length, bool allow_bridges)
{
  int num_blocks = 0;
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->is_valid(min_anchors, min_regions, min_length)) {
        if (all_links.insert(this_link).second) {
          num_blocks++;
        }
      } else if (allow_bridges and this_link->is_bridge(min_anchors, min_regions, min_length)) {
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
    \fn Graph::merge_alternative_paths(int max_anchors, std::string debug)
 */
int Graph::merge_alternative_paths(uint max_anchors, uint max_length, std::string debug)
{
  int count = 0;
  cout << "Merging alternative paths... (max anchors: " << max_anchors << ")" << endl;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    bool merge_event;
    Anchor *this_anchor = it->second;
    if (debug == "ALL" or this_anchor->id == debug) {
      cout << "Anchor before merging..." << endl;
      this_anchor->print();
    }
    do {
      merge_event = false;
      for (std::list<Link*>::iterator p_link1 = this_anchor->links.begin(); p_link1 != this_anchor->links.end() and !merge_event; p_link1++) {
        for (std::list<Link*>::iterator p_link2 = this_anchor->links.begin(); p_link2 != p_link1 and !merge_event; p_link2++) {
          if ((*p_link1)->is_an_alternative_path_of(*p_link2)) {
            if ((max_anchors == 0) or ((*p_link1)->get_num_of_mismatches(*p_link2) <= max_anchors)) {
              count++;
              if (debug == "ALL" or this_anchor->id == debug) {
                cout << "Merging these two paths:" << endl;
                (*p_link1)->print();
                (*p_link2)->print();
                cout << "^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
              }
              (*p_link1)->merge(*p_link2);
              merge_event = true;
            }
          }
        }
      }
    } while (merge_event);
    if (debug == "ALL" or this_anchor->id == debug) {
      cout << "Anchor after merging..." << endl;
      this_anchor->print();
    }
//     cout << endl;
  }
  cout << count << " merges." << endl;

  return count;
}


/*!
    \fn Graph::simplify(uint min_anchors, uint min_regions, uint min_length, std::string debug)
    Tries to get rid of small or spurious path matches. The idea is that two paths may join in a very small
    region and break a longer co-linear region. Let A-B-C-D and X-B-C-Y be two paths in the Graph. The
    minimize step would join both of them in B-C. If (B-C) does not fullfill all the requirements to be
    considered as a co-linear region, this method will try to split it in order to allow the minimization of
    (A-B-C-D) and (X-B-C-Y).
    -# Get all the unselected links with enough regions to be split and still be eligible as selected link
    -# See if any of these blocks can be split in order to enlarge adjacent blocks
       This happens when for any of the "front" Link, all the Link.tags are compatible with
       the ones in this Link and with all the ones of the "back" Link
 */
int Graph::simplify(uint min_anchors, uint min_regions, uint min_length, std::string debug)
{
  cout << "Simplifying graph..." << endl;
  // Get set of links that won't be selected as syntenic regions but contain enough regions to be split
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->tags.size() > min_regions and !this_link->is_valid(min_anchors, min_regions, min_length)) {
        all_links.insert(this_link);
      }
    }
  }

  uint split_count = 0;

  // See if any of these blocks can be split in order to enlarge adjacent blocks
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    Link* this_link = *p_link_it;
    Anchor *front_anchor = this_link->anchor_list.front();
    Anchor *back_anchor = this_link->anchor_list.back();

    // ==========================================================
    // Get the list of links for the front and the back anchors
    // ==========================================================
    set<Link*> front_links;
    set<Link*> back_links;
    for (list<Link*>::iterator p_front_link_it = front_anchor->links.begin(); p_front_link_it != front_anchor->links.end(); p_front_link_it++) {
      Link *this_front_link = *p_front_link_it;
      if (this_front_link != this_link
          and this_front_link->tags.size() < this_link->tags.size()
          and this_front_link->tags.size() >= min_regions) {
        front_links.insert(this_front_link);
      }
    }
    for (list<Link*>::iterator p_back_link_it = back_anchor->links.begin(); p_back_link_it != back_anchor->links.end(); p_back_link_it++) {
      Link *this_back_link = *p_back_link_it;
      if (this_back_link != this_link
          and this_back_link->tags.size() < this_link->tags.size()
          and this_back_link->tags.size() >= min_regions) {
        back_links.insert(this_back_link);
      }
    }
    // ==========================================================

    // ==========================================================
    // Print debug info if requested
    // ==========================================================
    bool print_debug_info = false;
    if (debug == "ALL" or (front_anchor->id == debug or back_anchor->id == debug)) {
      print_debug_info = true;
      cout << "------------------------------------" << endl;
      cout << " Test Simplifying ";
      this_link->print();
      if (front_anchor == back_anchor) {
        cout << "Front and back links (loop edge):" << endl;
        for (std::set<Link*>::iterator p_front_it = front_links.begin(); p_front_it != front_links.end(); p_front_it++) {
          if ((*p_front_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_front_it)->print();
        }
      } else {
        cout << "Front links:" << endl;
        for (std::set<Link*>::iterator p_front_it = front_links.begin(); p_front_it != front_links.end(); p_front_it++) {
          if ((*p_front_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_front_it)->print();
        }
        cout << "Back links:" << endl;
        for (std::set<Link*>::iterator p_back_it = back_links.begin(); p_back_it != back_links.end(); p_back_it++) {
          if ((*p_back_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_back_it)->print();
        }
      }
      cout << "------------------------------------" << endl;
    }
    // ==========================================================

    // test empty front or back links (or just 1?)
    if (front_links.size() == 0 or back_links.size() == 0) {
      continue;
    }
    // ==========================================================
    // spot matches front -- partial tested -- back
    //  this uses two interlaced loops
    // ==========================================================
    bool split = false;
    do {
      split = false;
      // FIRST LOOP: back anchor
      for (std::set<Link*>::iterator p_back_it = back_links.begin(); 
           !split and (p_back_it != back_links.end()); p_back_it++) {
        Link* back_link = *p_back_it;
        // the get_matching_tags method assumes this link comes before the other one.
        // Use strand 1 to test connectivity between the END of this link and the front_link:
        //  FrontAnchor---(this_link)---BackAnchor  <=?=> BackAnchor---(other_link)---AnyAnchor
        short this_strand = 1;
        short back_strand = back_link->get_strand_for_matching_tags(back_anchor);
        std::vector< std::list<tag>::iterator > this_tag_links_to_back =
            this_link->get_matching_tags(back_link, this_strand, back_strand, false);
        if (this_tag_links_to_back.empty()) {
          if (print_debug_info) cout << "empty back matching tag" << back_link->tags.size() << endl;
          continue;
        }

        // SECOND LOOP: front anchor
        for (std::set<Link*>::iterator p_front_it = front_links.begin();
             !split and (p_front_it != front_links.end()); p_front_it++) {
          Link* front_link = *p_front_it;
          if (front_link == back_link) {
//             cerr << "Self-link" << endl;
            continue;
            exit(11);
          }
          // the get_matching_tags method assumes this link comes before the other one.
          // Use strand -1 to test connectivity between the BEGINNING of this link and the front_link:
          //  BackAnchor---(this_link)---FrontAnchor  <=?=> FrontAnchor---(other_link)---AnyAnchor
          this_strand = -1;
          short front_strand = front_link->get_strand_for_matching_tags(front_anchor);
          std::vector< std::list<tag>::iterator > this_tag_links_to_front =
              this_link->get_matching_tags(front_link, this_strand, front_strand, false);
          if (this_tag_links_to_front.empty()) {
            if (print_debug_info) cout << "empty front matching tag" << front_link->tags.size() << endl;
            continue;
          }
          if (this_tag_links_to_front.size() != this_tag_links_to_back.size()) {
            cerr << "THIS SHOULD NEVER HAPPEN: lists for front and back do not have the same length: " <<
                this_tag_links_to_front.size() << " -- " << this_tag_links_to_back.size() << endl;
            exit(1);
          }
          bool should_be_split = true;
          uint front_matches = 0;
          uint back_matches = 0;
          uint number_of_matches = 0;
          vector<bool> tags_to_split(this_tag_links_to_front.size(), false);
          for (uint i=0; i< this_tag_links_to_front.size(); i++) {
            // this_tag_links_to_front[i] = front_link->tags.end() when no match has been found
            if (this_tag_links_to_front[i] != front_link->tags.end()) {
              front_matches++;
            }
            // this_tag_links_to_back[i] = back_link->tags.end() when no match has been found
            if (this_tag_links_to_back[i] != back_link->tags.end()) {
              back_matches++;
            }
            if (this_tag_links_to_front[i] != front_link->tags.end() and this_tag_links_to_back[i] != back_link->tags.end()) {
              number_of_matches++;
              tags_to_split[i] = true;
            } else if (this_tag_links_to_front[i] != front_link->tags.end() or this_tag_links_to_back[i] != back_link->tags.end()) {
              should_be_split = false;
  //             continue;
            }
          }
          if (print_debug_info) {
            std::list<tag>::iterator p_center_tag_it = this_link->tags.begin();
            for (uint i=0; i< this_tag_links_to_front.size(); i++) {
              cout << setw(2) << i+1 << " : ";
              if (this_tag_links_to_front[i] != front_link->tags.end()) {
                print_tag(*this_tag_links_to_front[i]);
              } else {
                cout << " ----------------------------------------- ";
              }
              cout << " :: ";
              print_tag(*(p_center_tag_it++));
              cout << " :: ";
              if (this_tag_links_to_back[i] != back_link->tags.end()) {
                print_tag(*this_tag_links_to_back[i]);
              } else {
                cout << " ----------------------------------------- ";
              }
              cout << endl;
            }
            cout << " ==> " << number_of_matches << " full match(es); details = "
                << front_matches << "/" << front_link->tags.size() << "; " << this_link->tags.size() << "; "
                << back_matches << "/" << back_link->tags.size() << "; should be split: "
                << (should_be_split and number_of_matches and number_of_matches < this_tag_links_to_front.size())
                << endl;
          }
          if (should_be_split and number_of_matches and number_of_matches < this_tag_links_to_front.size()) {
            // Split this link. A posterior loop of minimization will take care of joining the edges
            Link* new_link = this_link->split(tags_to_split);
            if (new_link) {
              split = true;
              split_count++;
              if (print_debug_info) {
                cout << "This link has been split in two parts:" << endl;
                new_link->print();
                this_link->print();
              }
            }
          }
        }
      }
    } while (split);
  }

  cout << split_count << " splits." << endl;

  return split_count;
}


/*!
    \fn Graph::simplify_aggressive(uint min_anchors, uint min_regions, uint min_length, std::string debug)
    Tries to get rid of small or spurious path matches. The idea is that two paths may join in a very small
    region and break a longer co-linear region. Let A-B-C-D, A-B-X and Y-C-D be three paths in the Graph. The
    minimize step would join the first two in A-B and the first and the third in C-D. If (A-B) and (C-D) do not
    fullfill all the requirements to be considered as a co-linear region, this method will try to split them
    in order to allow the minimization of (A-B-C-D).
    -# Get all the \link Link links \endlink with enough regions to be eligible as selected link
    -# See if two adjacent unselected \link Link links \endlink can be split in order to lengthen this one
 */
int Graph::simplify_aggressive(uint min_anchors, uint min_regions, uint min_length, std::string debug)
{
  cout << "Simplifying graph (aggressive method)..." << endl;
  // Get set of links that contain enough regions to be split
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      // AGGRESSIVE MODE: Choose all co-linear regions, even good ones.
      if (this_link->tags.size() >= min_regions) {
        all_links.insert(this_link);
      }
    }
  }

  uint split_count = 0;

  // See if any of these blocks can be split in order to enlarge adjacent blocks
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    Link* this_link = *p_link_it;
    Anchor *front_anchor = this_link->anchor_list.front();
    Anchor *back_anchor = this_link->anchor_list.back();

    // ==========================================================
    // Get the list of links for the front and the back anchors
    // ==========================================================
    set<Link*> front_links;
    set<Link*> back_links;
    for (list<Link*>::iterator p_front_link_it = front_anchor->links.begin();
         p_front_link_it != front_anchor->links.end(); p_front_link_it++) {
      Link *this_front_link = *p_front_link_it;
//       this_link->anchor_list.size() < min_anchors or this_link->get_shortest_length() < min_length
      if (this_front_link != this_link
          and this_front_link->tags.size() >= this_link->tags.size()
          and (this_front_link->anchor_list.size() < min_anchors or
          this_front_link->get_shortest_region_length() < min_length)) {
        front_links.insert(this_front_link);
      }
    }
    for (list<Link*>::iterator p_back_link_it = back_anchor->links.begin(); p_back_link_it != back_anchor->links.end(); p_back_link_it++) {
      Link *this_back_link = *p_back_link_it;
      if (this_back_link != this_link
          and this_back_link->tags.size() >= this_link->tags.size()
          and (this_back_link->anchor_list.size() < min_anchors or this_back_link->get_shortest_region_length() < min_length)) {
        back_links.insert(this_back_link);
      }
    }
    // ==========================================================

    // ==========================================================
    // Print debug info if requested
    // ==========================================================
    bool print_debug_info = false;
    if (debug == "ALL" or (front_anchor->id == debug or back_anchor->id == debug)) {
      print_debug_info = true;
      cout << "------------------------------------" << endl;
      cout << " Test Simplifying ";
      this_link->print();
      if (front_anchor == back_anchor) {
        cout << "Front and back links (loop edge):" << endl;
        for (std::set<Link*>::iterator p_front_it = front_links.begin(); p_front_it != front_links.end(); p_front_it++) {
          if ((*p_front_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_front_it)->print();
        }
      } else {
        cout << "Front links:" << endl;
        for (std::set<Link*>::iterator p_front_it = front_links.begin(); p_front_it != front_links.end(); p_front_it++) {
          if ((*p_front_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_front_it)->print();
        }
        cout << "Back links:" << endl;
        for (std::set<Link*>::iterator p_back_it = back_links.begin(); p_back_it != back_links.end(); p_back_it++) {
          if ((*p_back_it)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*p_back_it)->print();
        }
      }
      cout << "------------------------------------" << endl;
    }
    // ==========================================================

    // test empty front or back links (or just 1?)
    if (front_links.size() == 0 or back_links.size() == 0) {
      continue;
    }
    // ==========================================================
    // spot matches front -- partial tested -- back
    //  this uses two interlaced loops
    // ==========================================================
    bool split = false;
    do {
      split = false;
      // FIRST LOOP: back anchor
      for (std::set<Link*>::iterator p_back_it = back_links.begin(); 
           !split and (p_back_it != back_links.end()); p_back_it++) {
        Link* back_link = *p_back_it;
        // the get_matching_tags method assumes this link comes before the other one.
        // Use strand 1 to test connectivity between the END of this link and the front_link:
        //  FrontAnchor---(this_link)---BackAnchor  <=?=> BackAnchor---(other_link)---AnyAnchor
        short this_strand = 1;
        short back_strand = back_link->get_strand_for_matching_tags(back_anchor);
        std::vector< std::list<tag>::iterator > this_tag_links_to_back =
            this_link->get_matching_tags(back_link, this_strand, back_strand, true);
        if (this_tag_links_to_back.empty()) {
          if (print_debug_info) cout << "empty back matching tag" << back_link->tags.size() << endl;
          continue;
        }

        // SECOND LOOP: front anchor
        for (std::set<Link*>::iterator p_front_it = front_links.begin();
             !split and (p_front_it != front_links.end()); p_front_it++) {
          Link* front_link = *p_front_it;
          if (front_link == back_link) {
//             cerr << "Self-link" << endl;
            continue;
            exit(11);
          }
          // the get_matching_tags method assumes this link comes before the other one.
          // Use strand -1 to test connectivity between the BEGINNING of this link and the front_link:
          //  BackAnchor---(this_link)---FrontAnchor  <=?=> FrontAnchor---(other_link)---AnyAnchor
          this_strand = -1;
          short front_strand = front_link->get_strand_for_matching_tags(front_anchor);
          std::vector< std::list<tag>::iterator > this_tag_links_to_front =
              this_link->get_matching_tags(front_link, this_strand, front_strand, true);
          if (this_tag_links_to_front.empty()) {
            if (print_debug_info) cout << "empty front matching tag" << front_link->tags.size() << endl;
            continue;
          }
          if (this_tag_links_to_front.size() != this_tag_links_to_back.size()) {
            cerr << "THIS SHOULD NEVER HAPPEN: lists for front and back do not have the same length: " <<
                this_tag_links_to_front.size() << " -- " << this_tag_links_to_back.size() << endl;
            exit(1);
          }
          bool should_be_split = true;
          uint front_matches = 0;
          uint back_matches = 0;
          uint number_of_matches = 0;
          vector<bool> tags_to_split(this_tag_links_to_front.size(), false);
          for (uint i=0; i< this_tag_links_to_front.size(); i++) {
            // this_tag_links_to_front[i] = front_link->tags.end() when no match has been found
            if (this_tag_links_to_front[i] != front_link->tags.end()) {
              front_matches++;
            }
            // this_tag_links_to_back[i] = back_link->tags.end() when no match has been found
            if (this_tag_links_to_back[i] != back_link->tags.end()) {
              back_matches++;
            }
            if (this_tag_links_to_front[i] != front_link->tags.end() and this_tag_links_to_back[i] != back_link->tags.end()) {
              number_of_matches++;
              tags_to_split[i] = true;
            } else if (this_tag_links_to_front[i] != front_link->tags.end() or this_tag_links_to_back[i] != back_link->tags.end()) {
              should_be_split = false;
  //             continue;
            }
          }
          bool print_this_debug_info = false;
          should_be_split = false;
          if ((front_matches and front_matches == front_link->tags.size() and front_matches == number_of_matches) or
               (back_matches and back_matches == back_link->tags.size() and back_matches == number_of_matches)) {
//             print_this_debug_info = true;
//             cout << "SHOW THIS MATCH" << endl;
          }
          if (front_matches == number_of_matches and back_matches == number_of_matches and
             number_of_matches == this_link->tags.size()) {
//             print_this_debug_info = true;
            should_be_split = true;
//             cout << "SHOW THIS MATCH 2" << endl;
          }
          if (print_debug_info or print_this_debug_info) {
            std::list<tag>::iterator p_center_tag_it = this_link->tags.begin();
            for (uint i=0; i< this_tag_links_to_front.size(); i++) {
              cout << setw(2) << i+1 << " : ";
              if (this_tag_links_to_front[i] != front_link->tags.end()) {
                print_tag(*this_tag_links_to_front[i]);
              } else {
                cout << " ----------------------------------------- ";
              }
              cout << " :: ";
              print_tag(*(p_center_tag_it++));
              cout << " :: ";
              if (this_tag_links_to_back[i] != back_link->tags.end()) {
                print_tag(*this_tag_links_to_back[i]);
              } else {
                cout << " ----------------------------------------- ";
              }
              cout << endl;
            }
            cout << " ==> " << number_of_matches << " full match(es); details = "
                << front_matches << "/" << front_link->tags.size() << "; " << this_link->tags.size() << "; "
                << back_matches << "/" << back_link->tags.size() << "; should be split: "
                << (should_be_split and number_of_matches and number_of_matches < this_tag_links_to_front.size())
                << endl;
          }
          if (should_be_split and number_of_matches) {
            // Split this link. A posterior loop of minimization will take care of joining the edges
            if (front_matches < front_link->tags.size()) {
              Link* new_front_link = front_link->split(this_tag_links_to_front);
              if (new_front_link) {
                split = true;
                split_count++;
                if (print_debug_info) {
                  cout << "This link has been split in two parts (1/2):" << endl;
                  new_front_link->print();
                  front_link->print();
                }
              }
            }
            if (back_matches < back_link->tags.size()) {
              Link* new_back_link = back_link->split(this_tag_links_to_back);
              if (new_back_link) {
                split = true;
                split_count++;
                if (print_debug_info) {
                  cout << "This link has been split in two parts (2/2):" << endl;
                  new_back_link->print();
                  back_link->print();
                }
              }
            }
          }
        }
      }
    } while (split);
  }

  cout << split_count << " aggresive splits." << endl;

  return split_count;
}




/*!
    \fn Graph::split_unselected_links(uint min_anchors, uint min_regions, uint min_length, std::string debug)
    Tries to get rid of small or spurious path matches. Any unselected \link Link link \endlink will
    be split in order to be minimized. Hopefully this will allow to join
    -# Get all the unselected \link Link links \endlink
    -# Split them in unique regions
 */
int Graph::split_unselected_links(uint min_anchors, uint min_regions, uint min_length, std::string debug)
{
  cout << "Splitting unselected links..." << endl;
  // Get set of links that won't be selected as syntenic regions but contain enough regions to be splitted
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->tags.size() > min_regions and !this_link->is_valid(min_anchors, min_regions, min_length)) {
        all_links.insert(this_link);
      }
    }
  }

  uint split_count = 0;
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    Link* this_link = *p_link_it;
    while (this_link->tags.size() > 1) {
      vector<bool> tags_to_split(this_link->tags.size(), false);
      tags_to_split[0] = true;
      Link* new_link = this_link->split(tags_to_split);
      if (new_link) {
        split_count++;
      }
    }
  }

  cout << split_count << " splits (unselected links)" << endl;

  return split_count;
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



/*!
    \fn Graph::split_unbalanced_links(float max_ratio, std::string debug)
 */
void Graph::split_unbalanced_links(float max_ratio, std::string debug)
{
  if (max_ratio <= 1.0) {
    return;
  }

  cout << "Edit unbalanced links..." << endl;

  uint unbalanced_segments_counter = 0;
  uint unbalanced_links_counter = 0;

  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (this_link->tags.size() > 1) {
        std::map<std::string, uint> longest_segment;
        std::map<std::string, uint> shortest_segment;
        for (map<std::string, string*>::iterator p_species_it = this->species.begin(); p_species_it != this->species.end(); p_species_it++) {
          longest_segment[p_species_it->first] = 0;
          shortest_segment[p_species_it->first] = 0;
        }
//         uint longest_segment = 0;
//         uint shortest_segment = 0;
        for (list<tag>::iterator p_tag_it = this_link->tags.begin(); p_tag_it != this_link->tags.end(); p_tag_it++) {
          uint length = p_tag_it->end - p_tag_it->start + 1;
          if (length > longest_segment[*p_tag_it->species]) {
            longest_segment[*p_tag_it->species] = length;
          }
          if (shortest_segment[*p_tag_it->species] == 0) {
            shortest_segment[*p_tag_it->species] = length;
          } else if (length < shortest_segment[*p_tag_it->species]) {
            shortest_segment[*p_tag_it->species] = length;
          }
        }
        bool all_segments_are_balanced = true;
        for (map<std::string, string*>::iterator p_species_it = this->species.begin(); p_species_it != this->species.end(); p_species_it++) {
          if (shortest_segment[p_species_it->first] * max_ratio < longest_segment[p_species_it->first]) {
            all_segments_are_balanced = false;
          }
        }
        if (all_segments_are_balanced) {
          continue;
        }
        if (debug == "ALL" or this_anchor->id == debug) {
          cout << "Removing unbalanced links:" << endl;
          this_link->print();
        }
        list<tag> tmp_tags;
        for (list<tag>::iterator p_tag_it = this_link->tags.begin(); p_tag_it != this_link->tags.end(); p_tag_it++) {
          uint length = p_tag_it->end - p_tag_it->start + 1;
          if ((length * max_ratio) < longest_segment[*p_tag_it->species]) {
            unbalanced_segments_counter++;
            if (debug == "ALL" or this_anchor->id == debug) {
              cout << "Drop this: ";
              print_tag(*p_tag_it);
              cout << endl;
            }
          } else {
            tmp_tags.push_back(*p_tag_it);
          }
        }
        if (debug == "ALL" or this_anchor->id == debug) {
          cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
        }
        if (tmp_tags.size() == 0) {
          cout << "Leaving empty link" << endl;
          exit(1);
        }
        unbalanced_links_counter++;
        this_link->tags = tmp_tags;
      }
    }
  }
  cout << "removed " << unbalanced_segments_counter << " unbalanced segments in " << unbalanced_links_counter << " blocks" << endl;
}


/*!
    \fn Graph::resolve_small_palindromes(uint min_anchors, uint min_regions, uint min_length, std::string debug)
    Tries to get rid of small palindromes. Any unselected \link Link link \endlink
    will be tested to see if it corresponds to a palindrome. If it is, it will be resolved by removing
    the hairpin duplication and making a longer segment without the duplications
    -# Get all the unselected \link Link links \endlink
 */
uint Graph::resolve_small_palindromes(uint min_anchors, uint min_regions, uint min_length, std::string debug)
{
  cout << "Resolving small palindromes..." << endl;
  // Get set of circular links
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (!this_link->is_valid(min_anchors, min_regions, min_length) and fmod((float)this_link->tags.size(), 2.0f) == 0) {
        all_links.insert(this_link);
      }
    }
  }

  uint palindromes_count = 0;

  // See if any of these blocks is a small insertion breaking a large block
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    Link* this_link = *p_link_it;
    if (this_link->anchor_list.front()->id == debug and this_link->anchor_list.back()->id == debug) this_link->print();
    std::vector< std::list<tag>::iterator > this_tag_links_to_itself =
        this_link->get_matching_tags(this_link, 1, -1, false);
    vector<bool> tags_to_split(this_link->tags.size(), false);
    if (!this_tag_links_to_itself.empty()) {
//       this_link->print();
      std::list<tag>::iterator p_tag_it = this_link->tags.begin();
      for (uint i=0; i< this_link->tags.size(); i++) {
//         cout << i+1;
//         print_tag(*p_tag_it);
//         cout << " >>> ";
//         print_tag(*this_tag_links_to_itself[i]);
//         cout << endl;
        if (p_tag_it->start < this_tag_links_to_itself[i]->start) {
          tags_to_split[i] = true;
        }
        p_tag_it++;
      }
//       cout << endl;
    } else {
      this_tag_links_to_itself = this_link->get_matching_tags(this_link, -1, 1, false);
      if (!this_tag_links_to_itself.empty()) {
//         this_link->print();
        std::list<tag>::iterator p_tag_it = this_link->tags.begin();
        for (uint i=0; i< this_link->tags.size(); i++) {
//           cout << i+1;
//           print_tag(*p_tag_it);
//           cout << " >>> ";
//           print_tag(*this_tag_links_to_itself[i]);
//           cout << endl;
          if (p_tag_it->start < this_tag_links_to_itself[i]->start) {
            tags_to_split[i] = true;
          }
          p_tag_it++;
        }
//         cout << endl;
      }
    }
    if (!this_tag_links_to_itself.empty()) {
      Link* new_link = this_link->split(tags_to_split);
//       this_link->print();
//       new_link->print();
      if (this_link->try_to_concatenate_with(new_link)) {
        palindromes_count++;
      } else {
        cout << "FAIL!!!" << endl;
      }
      for (list<tag>::iterator p_tag1 = this_link->tags.begin(); p_tag1 != this_link->tags.end(); p_tag1++) {
        p_tag1->strand = 0;
      }
//       this_link->print();
    }
  }

  cout << "removed " << palindromes_count << " palindromic segments" << endl;

  return palindromes_count;
}


/*!
    \fn Graph::assimilate_small_insertions(uint min_anchors, uint min_regions, uint min_length, uint max_insertion_length, std::string debug)
    Tries to get rid of small or spurious insertions. Any unselected \link Link link \endlink
    will be tested to see if it corresponds to a small insertion. If such insertion is shorter
    than the max. insertion length, it will be assimilated into the path. A posterior loop of
    graph minimization will finish the work.
    -# Get all the unselected \link Link links \endlink
 */
uint Graph::assimilate_small_insertions(uint min_anchors, uint min_regions, uint min_length,
                                        uint max_insertion_length, std::string debug)
{
  cout << "Assimilating small insertions (max. insertion length: " << max_insertion_length << ")..." << endl;
  // Get set of circular links
  set<Link*> all_links;
  for (std::map<std::string, Anchor*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
    Anchor * this_anchor = it->second;
    for (list<Link*>::iterator p_link_it = this_anchor->links.begin(); p_link_it != this_anchor->links.end(); p_link_it++) {
      Link * this_link = *p_link_it;
      if (!this_link->is_valid(min_anchors, min_regions, min_length)) {
        all_links.insert(this_link);
      }
    }
  }

  uint assimilate_count = 0;

  // See if any of these blocks is a small insertion breaking a large block
  for (std::set<Link*>::iterator p_link_it = all_links.begin(); p_link_it != all_links.end(); p_link_it++) {
    Link* this_link = *p_link_it;
    Anchor *front_anchor = this_link->anchor_list.front();
    Anchor *back_anchor = this_link->anchor_list.back();

    // ==========================================================
    // Get the list of links for the front and the back anchors
    // ==========================================================
    set<Link*> front_links;
    set<Link*> back_links;
    for (list<Link*>::iterator p_front_link_it = front_anchor->links.begin(); p_front_link_it != front_anchor->links.end(); p_front_link_it++) {
      Link *this_front_link = *p_front_link_it;
      if (this_front_link != this_link and this_front_link->is_valid(min_anchors, min_regions, min_length)
         and this_front_link->tags.size() > this_link->tags.size()) {
        front_links.insert(this_front_link);
      }
    }
    for (list<Link*>::iterator p_back_link_it = back_anchor->links.begin(); p_back_link_it != back_anchor->links.end(); p_back_link_it++) {
      Link *this_back_link = *p_back_link_it;
      if (this_back_link != this_link and this_back_link->is_valid(min_anchors, min_regions, min_length)
          and this_back_link->tags.size() > this_link->tags.size()) {
        back_links.insert(this_back_link);
      }
    }
    // ==========================================================

    // ==========================================================
    // Print debug info if requested
    // ==========================================================
    bool print_debug_info = false;
    if (debug == "ALL" or (front_anchor->id == debug or back_anchor->id == debug)) {
      print_debug_info = true;
      cout << "------------------------------------" << endl;
      cout << " Test Assimilating ";
      this_link->print();
      if (front_anchor == back_anchor) {
        cout << "Front and back links (loop edge):" << endl;
        for (std::set<Link*>::iterator it1 = front_links.begin(); it1 != front_links.end(); it1++) {
          if ((*it1)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*it1)->print();
        }
      } else {
        cout << "Front links:" << endl;
        for (std::set<Link*>::iterator it1 = front_links.begin(); it1 != front_links.end(); it1++) {
          if ((*it1)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*it1)->print();
        }
        cout << "Back links:" << endl;
        for (std::set<Link*>::iterator it1 = back_links.begin(); it1 != back_links.end(); it1++) {
          if ((*it1)->is_valid(min_anchors, min_regions, min_length)) cout << "valid ";
          (*it1)->print();
        }
      }
      cout << "------------------------------------" << endl;
    }
    // ==========================================================

    // test empty front or back links (or just 1?)
    if (front_links.size() == 0 or back_links.size() == 0) {
      continue;
    }

    bool has_been_assimilated = false;

    // ==========================================================
    // spot small insertions
    // ==========================================================
    // FIRST LOOP: back anchor
    for (std::set<Link*>::iterator p_back_it = back_links.begin(); !has_been_assimilated and p_back_it != back_links.end(); p_back_it++) {
      Link* back_link = *p_back_it;

      // the get_matching_tags method assumes this link comes before the other one.
      // Use strand 1 to test connectivity between the END of this link and the front_link:
      //  FrontAnchor---(this_link)---BackAnchor  <=?=> BackAnchor---(other_link)---AnyAnchor
      short this_strand = 1;
      short back_strand = back_link->get_strand_for_matching_tags(back_anchor);
      std::vector< std::list<tag>::iterator > this_tag_links_to_back =
          this_link->get_matching_tags(back_link, this_strand, back_strand, true);
      if (this_tag_links_to_back.empty()) {
        if (print_debug_info) cout << "empty back matching tag " << back_link->tags.size() << " " << back_strand << endl;
        continue;
      }
      uint back_matches = 0;
      // I only want to go on if all the tags in this link match one from the back link
      for (uint i=0; i< this_link->tags.size(); i++) {
        // this_tag_links_to_front[i] = front_link->tags.end() when no match has been found
        if (this_tag_links_to_back[i] != back_link->tags.end()) back_matches++;
      }
      if (back_matches != this_link->tags.size()) continue;


      // SECOND LOOP: front anchor
      for (std::set<Link*>::iterator p_front_it = front_links.begin(); !has_been_assimilated and p_front_it != front_links.end(); p_front_it++) {
        Link* front_link = *p_front_it;
        if (front_link == back_link) continue;
        // I am looking for insertions between front and back link, therefore they must
        // have the same number of tags
        if (front_link->tags.size() != back_link->tags.size()) continue;

        // the get_matching_tags method assumes this link comes before the other one.
        // Use strand -1 to test connectivity between the BEGINNING of this link and the front_link:
        //  BackAnchor---(this_link)---FrontAnchor  <=?=> FrontAnchor---(other_link)---AnyAnchor
        this_strand = -1;
        short front_strand = front_link->get_strand_for_matching_tags(front_anchor);
        std::vector< std::list<tag>::iterator > this_tag_links_to_front =
            this_link->get_matching_tags(front_link, this_strand, front_strand, true);
        if (this_tag_links_to_front.empty()) {
          if (print_debug_info) cout << "empty front matching tag" << front_link->tags.size() << endl;
          continue;
        }
        if (this_tag_links_to_front.size() != this_tag_links_to_back.size()) {
          cerr << "THIS SHOULD NEVER HAPPEN: lists for front and back do not have the same length: " <<
              this_tag_links_to_front.size() << " -- " << this_tag_links_to_back.size() << endl;
          exit(1);
        }
        uint front_matches = 0;
        // I only want to go on if all the tags in this link match one from the back link
        for (uint i=0; i< this_tag_links_to_front.size(); i++) {
          // this_tag_links_to_front[i] = front_link->tags.end() when no match has been found
          if (this_tag_links_to_front[i] != front_link->tags.end()) front_matches++;
        }
        if (front_matches != this_link->tags.size()) continue;

        std::vector< std::list<tag>::iterator > front_tag_links_to_this =
            front_link->get_matching_tags(this_link, -front_strand, 1, true);
        std::vector< std::list<tag>::iterator > front_tag_links_to_back =
            front_link->get_matching_tags(back_link, -front_strand, back_strand, true);
        bool is_an_insertion = true;
        if (front_tag_links_to_back.empty() or front_tag_links_to_this.empty()) continue;
        for (uint i=0; i< front_link->tags.size(); i++) {
          if (front_tag_links_to_back[i] != back_link->tags.end() and front_tag_links_to_this[i] != this_link->tags.end()) {
            cout << "------------------------------------" << endl;
            cout << " WEIRD INSERTION ";
            front_link->print();
            this_link->print();
            back_link->print();
            cout << i << " matches this and back" << endl;
            is_an_insertion = false;
          } else if (front_tag_links_to_back[i] == back_link->tags.end() and
                     front_tag_links_to_this[i] == this_link->tags.end()) {
            is_an_insertion = false;
          }
        }
        if (!is_an_insertion) continue;
        uint longest_region_length = this_link->get_longest_region_length();
        if (longest_region_length > max_insertion_length) continue;
        if (print_debug_info) {
          cout << "------------------------------------" << endl;
          cout << " INSERTION ";
          front_link->print();
          this_link->print();
          back_link->print();
        }
        std::list<tag>::iterator p_front_tag_it = front_link->tags.begin();
        for (uint i=0; i< front_link->tags.size(); i++) {
          if (front_tag_links_to_this[i] != this_link->tags.end()) {
            if (p_front_tag_it->end < front_tag_links_to_this[i]->end) {
              p_front_tag_it->end = front_tag_links_to_this[i]->end;
            } else if (p_front_tag_it->start > front_tag_links_to_this[i]->start) {
              p_front_tag_it->start = front_tag_links_to_this[i]->start;
            }
          }
          p_front_tag_it++;
        }
        delete(this_link);
        has_been_assimilated = true;
        assimilate_count++;

      }
    }

  }

  cout << assimilate_count << " assimilated insertion." << endl;

  return assimilate_count;
}
