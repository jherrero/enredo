

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>

using namespace std;

bool read_file(char *filename, float min_score);
bool print_file(char *input_filename, char *output_filename, float min_score);
void print_help(void);

// Global variable
std::map<std::string, std::string* > anchors;

int main(int argc, char *argv[])
{
  char *filename = NULL;
  char *output_filename = NULL;
  float min_score = 0.0f;
  bool help = false;
  bool ret;
  string this_arg;

  for (int a = 1; a < argc; a++) {
    this_arg = argv[a];
    if (((this_arg == "--output") or (this_arg == "-o")) and (a < argc - 1)) {
      a++;
      output_filename = argv[a];
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

  ret = read_file(filename, min_score);
  ret = print_file(filename, output_filename, min_score);

  return EXIT_SUCCESS;
}

bool read_file(char *filename, float min_score)
{
  ifstream inputfile (filename);
  string line;
  if (!inputfile.is_open()) {
    cerr << "Cannot open file " << filename << endl;
    return false;
  }

  unsigned long long int line_counter = 0;
  string last_anchor_id = "";
  string last_species = "";
  string last_chr = "";
  int last_start = 0;
  int last_end = 0;
  string last_strand = "";

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
      cerr << "start cannot be larger than end in <" << line << ">" << endl;
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

    if (
        last_species == this_species and
        last_chr == this_chr and
        this_start <= last_end) {
      // These two anchors overlap
      cerr << this_anchor_id << " overlaps with " << last_anchor_id << ": ";
      if (!anchors[this_anchor_id] and !anchors[last_anchor_id]) {
        cerr << "new link " << this_anchor_id << " => " << last_anchor_id << endl;
        anchors[this_anchor_id] = new string(last_anchor_id);
      } else if (!anchors[this_anchor_id] and anchors[last_anchor_id]) {
        if (*anchors[last_anchor_id] != this_anchor_id) {
          cerr << "set " << this_anchor_id << " to " << *anchors[last_anchor_id] << endl;
          anchors[this_anchor_id] = anchors[last_anchor_id];
        } else {
          cerr << "(already existing)" << endl;
        }
      } else if (anchors[this_anchor_id] and !anchors[last_anchor_id]) {
        if (*anchors[this_anchor_id] != last_anchor_id) {
          cerr << "set " << last_anchor_id << " to " << *anchors[this_anchor_id] << endl;
          anchors[last_anchor_id] = anchors[this_anchor_id];
        } else {
          cerr << "(already existing)" << endl;
        }
      } else if (anchors[this_anchor_id] and anchors[last_anchor_id]) {
        if (*anchors[this_anchor_id] != *anchors[last_anchor_id]) {
          cerr << "both link to different anchors (" << *anchors[this_anchor_id] << " and " << *anchors[last_anchor_id] << "); set all to " << *anchors[this_anchor_id] << endl;
          string former_merged_id = *anchors[last_anchor_id];
          anchors[former_merged_id] = anchors[this_anchor_id];
          for (std::map<std::string, std::string*>::iterator it = anchors.begin(); it != anchors.end(); it++) {
            if (it->second and (*(it->second) == former_merged_id)) {
              it->second = anchors[this_anchor_id];
            }
          }
        } else {
          cerr << "both link to " << *anchors[this_anchor_id] << endl;
        }
      }

    }

    // Set last_* variables to current ones before next loop
    last_anchor_id = this_anchor_id;
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
  return true;
}

bool print_file(char *input_filename, char *output_filename, float min_score) {
  ifstream inputfile (input_filename);
  ofstream outputfile;
  ostream *out;
  if (output_filename) {
    outputfile.open(output_filename);
    if (!outputfile.is_open()) {
      cerr << "Cannot open file " << output_filename << endl;
      return false;
    }
    out = &outputfile;
  } else {
    out = &cout;
  }
  string line;
  if (!inputfile.is_open()) {
    cerr << "Cannot open file " << output_filename << endl;
    return false;
  }

  unsigned long long int line_counter = 0;
  string last_anchor_id = "";
  string last_species = "";
  string last_chr = "";
  int last_start = 0;
  int last_end = 0;
  string last_strand = "+";
  float last_score = 0.0f;
  while (!inputfile.eof()) {
    getline(inputfile, line);
    if (inputfile.eof()) {
      if (anchors[last_anchor_id]) {
        *out << *anchors[last_anchor_id] << "\t" << last_species << "\t"
            << last_chr << "\t" << last_start << "\t" << last_end << "\t"
            << last_strand << "\t" << last_score << endl;
      } else if (last_anchor_id != "") {
        *out << last_anchor_id << "\t" << last_species << "\t"
            << last_chr << "\t" << last_start << "\t" << last_end << "\t"
            << last_strand << "\t" << last_score << endl;
      }
      break;
    }
    if (line[0] == '#') {
      *out << line << endl;
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
    if (this_score < min_score) {
      *out << "#LOW_SCORE:" << line << endl;
      continue;
    }
    if (
        last_species == this_species and
        last_chr == this_chr and
        this_start <= last_end) {
      if (this_end > last_end) {
        last_end = this_end;
      }
      if (this_score > last_score) {
        last_score = this_score;
      }
    } else {
      if (anchors[last_anchor_id]) {
        *out << *anchors[last_anchor_id] << "\t" << last_species << "\t"
            << last_chr << "\t" << last_start << "\t" << last_end << "\t"
            << last_strand << "\t" << last_score << endl;
      } else if (last_anchor_id != "") {
        *out << last_anchor_id << "\t" << last_species << "\t"
            << last_chr << "\t" << last_start << "\t" << last_end << "\t"
            << last_strand << "\t" << last_score << endl;
      }

      // Set last_* variables to current ones before next loop
      last_anchor_id = this_anchor_id;
      last_species = this_species;
      last_chr = this_chr;
      last_start = this_start;
      last_end = this_end;
      last_strand = this_strand;
      last_score = this_score;
    }

    line_counter++;
  }
  inputfile.close();
  if (output_filename) {
    outputfile.close();
  }
  return true;
}

void print_help(void)
{
  cout << "MergeOverlap v" << VERSION << endl;
  cout << endl;
  cout << "Usage: mergeoverlap anchors_file.txt" << endl;
  cout << endl;
  cout << " --help: prints this help" << endl;
  cout << endl;
  cout << "See README file for more details." << endl;
  cout << endl;
}
