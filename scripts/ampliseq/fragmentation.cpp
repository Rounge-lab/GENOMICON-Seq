// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <random>
#include <fstream>

using namespace Rcpp;

struct Fragment {
  int start;
  int end;
};

// [[Rcpp::export]]
List process_sequence_cpp(int sequence_length, int fragment_length_min, int fragment_length_max, int seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(fragment_length_min, fragment_length_max);
  std::uniform_int_distribution<> dis2(1, sequence_length - fragment_length_min + 1);
  
  std::vector<Fragment> fragments;
  while (sequence_length >= fragment_length_min) {
    int frag_length = dis(gen);
    frag_length = std::min(frag_length, sequence_length);
    int start_pos = dis2(gen);
    int end_pos = start_pos + frag_length - 1;
    
    fragments.push_back(Fragment{start_pos, end_pos});
    
    sequence_length -= frag_length;
  }
  
  // Convert the vector of Fragments to a List of IntegerVectors and return it
  List result;
  for (const Fragment& frag : fragments) {
    result.push_back(IntegerVector::create(frag.start, frag.end));
  }
  return result;
}

// [[Rcpp::export]]
void expand_table_cpp(DataFrame table, std::string table_name, int genome_length, int fragment_length_min, int fragment_length_max, int seed, std::string output_folder) {
  // Replace _sample_table with _fragment_coordinates in the table name
  size_t pos = table_name.find("_sample_table");
  if (pos != std::string::npos) {
    table_name.replace(pos, 13, "_fragment_coordinates");
  }
  
  // Create the full file path
  std::string file_path = output_folder + "/" + table_name + ".csv";
  
  std::ofstream file(file_path);
  
  // Write the column names
  file << "genome_name_copy;coordinates\n";
  
  // Extract the columns from the DataFrame
  CharacterVector genome_names = table["genome_name"];
  IntegerVector copy_numbers = table["copy_number"];
  
  for (int i = 0; i < table.nrows(); ++i) {
    for (int j = 0; j < copy_numbers[i]; ++j) {
      List fragments = process_sequence_cpp(genome_length, fragment_length_min, fragment_length_max, seed + i + j);
      file << genome_names[i] << "_c" << j + 1 << ";";
      for (int k = 0; k < fragments.size(); ++k) {
        IntegerVector frag = fragments[k];
        if (k > 0) {
          file << ",";
        }
        file << frag[0] << ":" << frag[1];
      }
      file << "\n";
    }
  }
  file.close();
}



