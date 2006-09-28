

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "graph.h"

using namespace std;

void print_help(void);

int main(int argc, char *argv[])
{
  Graph my_graph;
  char *filename = NULL;
  uint max_gap_length = 10000;
  float min_score = 0.0f;
  uint min_length = 100000;
  uint min_regions = 2;
  uint min_anchors = 3;
  uint path_dissimilarity = 0;
  int histogram_size = 10;
  bool print_all = false;
  bool help = false;
  bool ret;
  string this_arg;

  for (int a = 1; a < argc; a++) {
    this_arg = argv[a];
    if ((this_arg == "--max-gap-length") and (a < argc - 1)) {
      a++;
      max_gap_length = atoi(argv[a]);
    } else if ((this_arg == "--min-score") and (a < argc - 1)) {
      a++;
      min_score = atof(argv[a]);
    } else if ((this_arg == "--min-length") and (a < argc - 1)) {
      a++;
      min_length = atoi(argv[a]);
    } else if ((this_arg == "--min-regions") and (a < argc - 1)) {
      a++;
      min_regions = atoi(argv[a]);
    } else if ((this_arg == "--min-anchors") and (a < argc - 1)) {
      a++;
      min_anchors = atoi(argv[a]);
    } else if (this_arg == "--all") {
      print_all = true;
    } else if ((this_arg == "--histogram-size") and (a < argc - 1)) {
      a++;
      histogram_size = atoi(argv[a]);
    } else if ((this_arg == "--help") or (this_arg == "-h")) {
      help = true;
    } else if (!filename) {
      filename = argv[a];
    } else {
      cerr << "Unknown option: " << this_arg << endl;
    }
  }

  if (help or !filename) {
    print_help();
    exit(0);
  } else if (print_all) {
    min_length = 1;
    min_regions = 1;
    min_anchors = 1;
  }

  cout << "SyntenyGraph v" << VERSION << endl;
  cout << endl
      << " Parameters:" << endl
      << "====================================" << endl
      << "Input-file: " << filename << endl
      << "min-score: " << min_score << endl
      << "max-gap-length: " << max_gap_length << endl
      << "min-length: " << min_length << endl
      << "min-regions: " << min_regions << endl
      << "min-anchors: " << min_anchors << endl;

  cout << endl
      << " Reading input file:" << endl
      << "====================================" << endl;
  ret = my_graph.populate_from_file(filename, min_score, max_gap_length);
  if (!ret) {
    cerr << "EXIT (Error while reading file)" << endl;
    exit(1);
  }


  cout << endl
      << " Stats before minimizing the Graph:" << endl
      << "====================================" << endl;
  my_graph.print_stats(histogram_size);

  my_graph.minimize();
  for (uint a = 1; a < path_dissimilarity; a++) {
    my_graph.merge_alternative_paths(a);
    my_graph.minimize();
  }

  cout << endl
      << " Stats after minimizing the Graph:" << endl
      << "===================================" << endl;
  my_graph.print_stats(histogram_size);

  cout << endl
      << " Resulting blocks:" << endl
      << "===================================" << endl;
  my_graph.print_links(min_anchors, min_regions, min_length);

  return EXIT_SUCCESS;
}


void print_help(void)
{
  cout << "SyntenyGraph v" << VERSION << endl;
  cout << endl;
  cout << "Usage: sytnenygraph [options] anchors_file.txt" << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << " --max-gap-length: maximum allowed gap between two anchors" << endl;
  cout << endl;
  cout << " --min-length: minimum length of final syntenic block (def: 100000)" << endl;
  cout << " --min-regions: minimum number of region in the syntenic block (def: 2)" << endl;
  cout << " --min-anchors: minimum number of anchors in the syntenic block (def: 3)" << endl;
  cout << " --all: print all the syntenic block (overwrite previous values)" << endl;
  cout << endl;
  cout << " --histogram-size: size for histogram of num. of regions pero link (def: 10)" << endl;
  cout << endl;
  cout << " --help: prints this help" << endl;
  cout << endl;
  cout << "See README file for more details." << endl;
  cout << endl;
}
