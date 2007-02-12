

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
  char *input_filename = NULL;
  char *output_filename = NULL;
  uint max_gap_length = 100000;
  float min_score = 0.0f;
  uint min_length = 100000;
  uint min_regions = 2;
  uint min_anchors = 3;
  uint path_dissimilarity = 0;
  bool simplify_graph = false;
  int histogram_size = 10;
  bool print_all = false;
  bool print_stats = true;
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
    } else if ((this_arg == "--max-path-dissimilarity") and (a < argc - 1)) {
      a++;
      path_dissimilarity = atoi(argv[a]);
    } else if ((this_arg == "--simplify-graph")) {
      simplify_graph = true;
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
    } else if (this_arg == "--no-stats") {
      print_stats = false;
    } else if (((this_arg == "--output-file") or (this_arg == "--output") or (this_arg == "-o"))and (a < argc - 1)) {
      a++;
      output_filename  = argv[a];
    } else if ((this_arg == "--help") or (this_arg == "-h")) {
      help = true;
    } else if (!input_filename) {
      input_filename = argv[a];
    } else {
      cerr << "Unknown option: " << this_arg << endl;
      exit(1);
    }
  }

  if (help or !input_filename) {
    print_help();
    exit(0);
  }

  cout << "Enredo v" << VERSION << endl;
  cout << endl
      << " Parameters:" << endl
      << "====================================" << endl
      << "Input-file: " << input_filename << endl
      << "min-score: " << min_score << endl
      << "max-gap-length: " << max_gap_length << endl
      << "max-path-dissimilarity: " << path_dissimilarity << endl
      << "min-length: " << min_length << endl
      << "min-regions: " << min_regions << endl
      << "min-anchors: " << min_anchors << endl
      << "simplify-graph: " << (simplify_graph?"yes":"no") << endl
      << "print-all: " << (print_all?"yes":"no") << endl;

  cout << endl
      << " Reading input file:" << endl
      << "====================================" << endl;
  ret = my_graph.populate_from_file(input_filename, min_score, max_gap_length);
  if (!ret) {
    cerr << "EXIT (Error while reading file)" << endl;
    exit(1);
  }
  if (print_stats) {
    my_graph.print_anchors_histogram();
  }


  cout << endl;
  if (print_stats) {
    cout << " Stats before minimizing the Graph:" << endl
        << "====================================" << endl;
    my_graph.print_stats(histogram_size);
  }

  my_graph.minimize();
  for (uint a = 0; a < path_dissimilarity; a++) {
    my_graph.merge_alternative_paths(a + 1);
    my_graph.minimize();
  }

  if (print_stats) {
    cout << endl
        << " Stats after minimizing the Graph:" << endl
        << "===================================" << endl;
    my_graph.print_stats(histogram_size);
  }

  my_graph.study_anchors();
  if (simplify_graph) {
    my_graph.simplify(min_anchors, min_regions, min_length);
    my_graph.minimize();
    if (print_stats) {
      cout << endl
          << " Stats after simplifying the Graph:" << endl
          << "===================================" << endl;
//       my_graph.print_stats(histogram_size);
    }
    my_graph.study_anchors();
  }
  if (print_all) {
    min_length = 1;
    min_regions = 1;
    min_anchors = 1;
  }

  cout << endl
      << " Resulting blocks:" << endl
      << "===================================" << endl;
  unsigned long int num_of_blocks;
  if (output_filename) {
    ofstream output_stream(output_filename);
    if (!output_stream.is_open()) {
      cerr << "EXIT (Cannot open file <" << output_filename << "> for output)" << endl;
      exit(1);
    }
    cout << "Results in file <" << output_filename << ">" << endl;
    output_stream << "## Enredo v" << VERSION << endl
        << "#" << endl
        << "#  Parameters:" << endl
        << "# ====================================" << endl
        << "# Input-file: " << input_filename << endl
        << "# min-score: " << min_score << endl
        << "# max-gap-length: " << max_gap_length << endl
        << "# max-path-dissimilarity: " << path_dissimilarity << endl
        << "# min-length: " << min_length << endl
        << "# min-regions: " << min_regions << endl
        << "# min-anchors: " << min_anchors << endl
        << "# simplify-graph: " << (simplify_graph?"yes":"no") << endl
        << "# print-all: " << (print_all?"yes":"no") << endl
        << endl;
    num_of_blocks = my_graph.print_links(output_stream, min_anchors, min_regions, min_length);
    output_stream.close();
  } else {
    num_of_blocks = my_graph.print_links(cout, min_anchors, min_regions, min_length);
  }
  cout << " Got " << num_of_blocks << " blocks." << endl;

  return EXIT_SUCCESS;
}


void print_help(void)
{
  cout << "Enredo v" << VERSION << endl
      << endl
      << "Usage: enredo [options] anchors_file.txt" << endl
      << endl
      << "Options:" << endl
      << " --max-gap-length: maximum allowed gap between two anchors (def: 100000)"  << endl
      << " --min-score: minimum score required to accept a hit" << endl
      << endl
      << " --max-path-dissimilarity: merge alternative paths in the graph if their" << endl
      << "       dissimilarity is up to this threshold (def: 0)" << endl
      << " --simplify-graph: try to split small edges in order to lengthen" << endl
      << "       other syntenic blocks (def. false)" << endl
      << endl
      << " --min-length: minimum length of final syntenic block (def: 100000)" << endl
      << " --min-regions: minimum number of region in the syntenic block (def: 2)" << endl
      << " --min-anchors: minimum number of anchors in the syntenic block (def: 3)" << endl
      << " --all: print all the syntenic block (overwrite previous values)" << endl
      << endl
      << " --histogram-size: size for histogram of num. of regions pero link (def: 10)" << endl
      << endl
      << " --help: prints this help" << endl
      << endl
      << "See README file for more details." << endl
      << endl;
}
