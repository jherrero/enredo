

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "graph.h"

using namespace std;

int main(int argc, char *argv[])
{
  Graph my_graph;
  char *filename = NULL;
  int min_length = 100000;
  int min_regions = 2;
  int min_anchors = 3;
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
    } else if (!filename) {
      filename = argv[a];
    } else {
      cerr << "Unknown option: " << this_arg << endl;
    }
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
