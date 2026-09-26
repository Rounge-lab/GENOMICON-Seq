// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace Rcpp;

struct Fragment {
    int start;
    int end;
};


// -----------------------------------------------------------------------------
// Fragment one genome copy.
//
// The genome is treated as a set of remaining unfragmented regions.
// From each region:
//   1. choose a random fragment length
//   2. choose a random start position
//   3. record the fragment
//   4. split the remaining sequence into left and right regions
//   5. continue until no remaining region can produce a fragment >= frag_min
//
// Circularity is handled outside this function through the "new_start"
// coordinates generated in fragmentation.R and subsequently converted to
// genomic coordinates by csv_to_bedgz.
// -----------------------------------------------------------------------------

std::vector<Fragment> make_fragments(int genome_len,
                                     int frag_min,
                                     int frag_max,
                                     std::mt19937 &gen)
{
    std::vector<Fragment> fragments;

    // Nothing can be fragmented if the genome is shorter than minimum size.
    if (genome_len < frag_min) {
        return fragments;
    }

    // Regions that are still available for fragmentation.
    std::vector<Fragment> remaining_regions;

    remaining_regions.push_back({1, genome_len});


    while (!remaining_regions.empty()) {

        // Take one remaining region.
        Fragment region = remaining_regions.back();
        remaining_regions.pop_back();

        int region_len = region.end - region.start + 1;

        // Region is too short to produce another fragment.
        if (region_len < frag_min) {
            continue;
        }


        // Fragment cannot be longer than the remaining region.
        int current_max = std::min(frag_max, region_len);

        std::uniform_int_distribution<int> length_distribution(
            frag_min,
            current_max
        );

        int fragment_length = length_distribution(gen);


        // Pick a valid random start within this remaining region.
        int latest_start =
            region.end - fragment_length + 1;

        std::uniform_int_distribution<int> start_distribution(
            region.start,
            latest_start
        );

        int fragment_start = start_distribution(gen);
        int fragment_end =
            fragment_start + fragment_length - 1;


        // Store generated fragment.
        fragments.push_back({
            fragment_start,
            fragment_end
        });


        // -------------------------------------------------------------
        // Sequence remaining BEFORE the selected fragment.
        // -------------------------------------------------------------

        int left_start = region.start;
        int left_end   = fragment_start - 1;

        int left_length =
            left_end - left_start + 1;

        if (left_length >= frag_min) {
            remaining_regions.push_back({
                left_start,
                left_end
            });
        }


        // -------------------------------------------------------------
        // Sequence remaining AFTER the selected fragment.
        // -------------------------------------------------------------

        int right_start = fragment_end + 1;
        int right_end   = region.end;

        int right_length =
            right_end - right_start + 1;

        if (right_length >= frag_min) {
            remaining_regions.push_back({
                right_start,
                right_end
            });
        }
    }


    // Sorting is not required for fragmentation itself,
    // but keeps fragment coordinates ordered and output predictable.
    std::sort(
        fragments.begin(),
        fragments.end(),
        [](const Fragment &a, const Fragment &b) {
            return a.start < b.start;
        }
    );


    return fragments;
}


// -----------------------------------------------------------------------------
// Main Rcpp worker
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
void expand_table_cpp(Rcpp::DataFrame table,
                      std::string table_name,
                      int genome_length,
                      int frag_min,
                      int frag_max,
                      int seed,
                      bool circular,
                      std::string output_folder)
{
    // Keep argument for compatibility with fragmentation.R.
    //
    // Circular genomes are handled by assigning each genome copy a random
    // new start position in fragmentation.R. The downstream BED conversion
    // performs coordinate rotation/wrapping.
    (void)circular;


    // -------------------------------------------------------------------------
    // Rename output file:
    // genome_sample_table -> genome_fragment_coordinates
    // -------------------------------------------------------------------------

    size_t pos = table_name.find("_sample_table");

    if (pos != std::string::npos) {
        table_name.replace(
            pos,
            13,
            "_fragment_coordinates"
        );
    }


    std::string file_path =
        output_folder + "/" + table_name + ".csv";

    std::ofstream file(file_path);


    if (!file.is_open()) {
        Rcpp::stop(
            "Unable to open fragmentation output file: " +
            file_path
        );
    }


    file << "genome_name_copy;coordinates\n";


    // -------------------------------------------------------------------------
    // Input sample table
    // -------------------------------------------------------------------------

    Rcpp::CharacterVector genome_names =
        table["genome_name"];

    Rcpp::IntegerVector copy_numbers =
        table["copy_number"];


    // One RNG for the whole fragmentation replicate.
    std::mt19937 gen(seed);


    // -------------------------------------------------------------------------
    // Process every genome type and every physical genome copy.
    // -------------------------------------------------------------------------

    for (int i = 0; i < table.nrows(); ++i) {

        for (int j = 0; j < copy_numbers[i]; ++j) {


            // Fragment the complete genome copy.
            std::vector<Fragment> fragments =
                make_fragments(
                    genome_length,
                    frag_min,
                    frag_max,
                    gen
                );


            // Genome-copy identifier.
            file
                << genome_names[i]
                << "_c"
                << j + 1
                << ";";


            // Store all coordinates belonging to this copy.
            //
            // Example:
            //
            // chr7_NMG_c1;1:178,179:421,422:615,...
            //
            // csv_to_bedgz will later turn every coordinate pair into
            // an individual fragment and assign _f1, _f2, _f3, ...
            //
            for (size_t k = 0; k < fragments.size(); ++k) {

                if (k > 0) {
                    file << ",";
                }

                file
                    << fragments[k].start
                    << ":"
                    << fragments[k].end;
            }


            file << "\n";
        }
    }


    file.close();
}
