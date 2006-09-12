

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
  int min_length = 100000;
  int min_regions = 2;
  int min_anchors = 3;
  bool help = false;
  bool ret;
  string this_arg;

  for (int a = 1; a < argc; a++) {
    this_arg = argv[a];
    if ((this_arg == "--min-length") and (a < argc - 1)) {
      a++;
      min_length = atoi(argv[a]);
    } else if ((this_arg == "--min-regions") and (a < argc - 1)) {
      a++;
      min_regions = atoi(argv[a]);
    } else if ((this_arg == "--min-anchors") and (a < argc - 1)) {
      a++;
      min_anchors = atoi(argv[a]);
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
  }

  ret = my_graph.populate_from_file(filename);
  if (!ret) {
    cerr << "EXIT (Error while reading file)" << endl;
    exit(1);
  }
//   my_graph.print();

  my_graph.minimize();
  my_graph.merge_alternative_paths(1);
//   string kk;
//   cin >> kk;
  my_graph.minimize();
  my_graph.merge_alternative_paths(2);
//   cin >> kk;
  my_graph.minimize();
//   my_graph.merge_alternative_paths(3);
//   cin >> kk;
//   my_graph.minimize();

  my_graph.print();
  my_graph.print_links(min_anchors, min_regions, min_length);

  return EXIT_SUCCESS;
}


void print_help(void)
{
  cout << "SyntenyGraph v0.1" << endl;
  cout << endl;
  cout << "Usage: sytnenygraph [options] anchors_file.txt" << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << " --min-length: minimum length of final syntenic block (def: 100000)" << endl;
  cout << " --min-regions: minimum number of region in the syntenic block (def: 2)" << endl;
  cout << " --min-anchors: minimum number of anchors in the syntenic block (def: 3)" << endl;
  cout << " --help: prints this help" << endl;
  cout << endl;
  cout << "See README file for more details." << endl;
  cout << endl;
}
