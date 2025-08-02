#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <cstdio>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <limits>

std::vector<std::string> split(const std::string& str, char delim);
void concatenate_csv(const std::string& directory, const std::string& output_file);
double calculate_sum(const std::string& input_file);
void process_csv_to_txt(const std::string& input_file, const std::string& output_file, double sum);

// Split a string by a delimiter into a vector of strings
std::vector<std::string> split(const std::string& str, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

void concatenate_csv(const std::string& directory, const std::string& output_file) {
    std::ofstream out(output_file, std::ios::app); 

    if (!out.is_open()) {
        std::cerr << "Unable to open output file: " << output_file << std::endl;
        return;
    }

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(directory.c_str())) != nullptr) {
        while ((ent = readdir(dir)) != nullptr) {
            std::string filename = ent->d_name;

            if (filename.find("_temp_") != std::string::npos) {
                std::ifstream in(directory + "/" + filename);

                if (!in.is_open()) {
                    std::cerr << "Unable to open input file: " << directory + "/" + filename << std::endl;
                    continue;
                }

                std::string header;
                std::getline(in, header); 

                std::string line;
                while (std::getline(in, line)) {
                    out << line << "\n";
                }

                in.close();

                if (std::remove((directory + "/" + filename).c_str()) != 0) {
                    std::cerr << "Error deleting file: " << directory + "/" + filename << std::endl;
                }
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Cannot open directory: " << directory << std::endl;
    }

    out.close();
}

double calculate_sum(const std::string& input_file) {
    std::ifstream in(input_file);
    if (!in.is_open()) {
        std::cerr << "Unable to open input file for sum calculation: " << input_file << std::endl;
        return 0;
    }

    std::string line;
    double sum = 0.0;
    while (std::getline(in, line)) {
        auto tokens = split(line, ';');
        if (tokens.size() >= 6) {
            try {
                sum += std::stod(tokens[5]);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number encountered in column 6: " << tokens[5] << std::endl;
            }
        }
    }

    in.close();
    return sum;
}

void process_csv_to_txt(const std::string& input_file, const std::string& output_file, double sum) {
    std::ifstream in(input_file);
    std::ofstream out(output_file, std::ios::trunc); // Overwrite mode

    if (!in.is_open()) {
        std::cerr << "Unable to open input file: " << input_file << std::endl;
        return;
    }
    if (!out.is_open()) {
        std::cerr << "Unable to open output file: " << output_file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        auto tokens = split(line, ';');
        if (tokens.size() >= 6) {
            std::string combined = tokens[0] + "_" + tokens[4];
            std::replace(combined.begin(), combined.end(), ' ', '_');
            std::replace(combined.begin(), combined.end(), ',', '_');

            double normalizedValue = std::stod(tokens[5]) / sum;
            out << combined << "\t" << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << normalizedValue << "\n";
        }
    }

    in.close();
    out.close();
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <directory> <merged_output_file> <processed_output_txt>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string merged_output_file = argv[2];
    std::string processed_output_txt = argv[3];

    concatenate_csv(directory, merged_output_file);

    double sum = calculate_sum(merged_output_file);
    if (sum == 0) {
        std::cerr << "Sum of the sixth column is zero, cannot proceed with division." << std::endl;
        return 1;
    }

    process_csv_to_txt(merged_output_file, processed_output_txt, sum);

    return 0;
}
