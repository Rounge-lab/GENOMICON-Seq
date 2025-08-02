// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <random>
#include <fstream>

using namespace Rcpp;

struct Fragment { int start; int end; };

// -----------------------------------------------------------------------------
// helper: create one fragment (or two sub-fragments if it wraps)
// -----------------------------------------------------------------------------
std::vector<Fragment> make_fragment(int genome_len,
                                    int frag_min,
                                    int frag_max,
                                    std::mt19937 &gen,
                                    bool circular)
{
    std::uniform_int_distribution<> len_dis(frag_min, frag_max);
    int frag_len = len_dis(gen);

    std::uniform_int_distribution<> start_dis(1, genome_len);
    int start = start_dis(gen);

    // -------- linear genome: retry until fragment fits completely ------------
    if (!circular) {
        while (start + frag_len - 1 > genome_len)
            start = start_dis(gen);
        return { {start, start + frag_len - 1} };
    }

    // -------- circular genome -------------------------------------------------
    int end = start + frag_len - 1;

    // fits without wrapping
    if (end <= genome_len)
        return { {start, end} };

    // split into two segments
    int end_first  = genome_len;
    int len_second = end - genome_len;          // overflow length
    return { {start, end_first}, {1, len_second} };
}

// -----------------------------------------------------------------------------
// C++ worker exported to R
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void expand_table_cpp(Rcpp::DataFrame  table,
                      std::string      table_name,
                      int              genome_length,
                      int              frag_min,
                      int              frag_max,
                      int              seed,
                      bool             circular,        // <-- NEW
                      std::string      output_folder)
{
    // rename output file
    size_t pos = table_name.find("_sample_table");
    if (pos != std::string::npos)
        table_name.replace(pos, 13, "_fragment_coordinates");

    std::string file_path = output_folder + "/" + table_name + ".csv";
    std::ofstream file(file_path);
    file << "genome_name_copy;coordinates\n";

    Rcpp::CharacterVector genome_names = table["genome_name"];
    Rcpp::IntegerVector   copy_numbers = table["copy_number"];

    std::mt19937 gen(seed);

    // loop over rows and genome copies
    for (int i = 0; i < table.nrows(); ++i) {
        for (int j = 0; j < copy_numbers[i]; ++j) {

            std::vector<Fragment> frags =
                make_fragment(genome_length, frag_min, frag_max, gen, circular);

            file << genome_names[i] << "_c" << j + 1 << ";";

            for (size_t k = 0; k < frags.size(); ++k) {
                if (k) file << ",";
                file << frags[k].start << ":" << frags[k].end;
            }
            file << "\n";
        }
    }
    file.close();
}
