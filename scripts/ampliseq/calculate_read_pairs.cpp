#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <boost/program_options.hpp>
#include <random>


namespace po = boost::program_options;

std::map<std::string, double> readAbundanceFile(const std::string& filename) {
    std::map<std::string, double> abundanceData;
    std::ifstream file(filename);
    std::string line, name;
    double value;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> name >> value)) { break; }
        abundanceData[name] = value;
    }
    return abundanceData;
}

std::map<std::string, int> readFastaLengthsFile(const std::string& filename) {
    std::map<std::string, int> fastaLengths;
    std::ifstream file(filename);
    std::string line, name;
    int length;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> name >> length)) { break; }
        fastaLengths[name] = length;
    }
    return fastaLengths;
}

double invertedSigmoid(int length, double k_value, double midpoint) {
    return 1.0 / (1.0 + exp(k_value * (static_cast<double>(length) - midpoint)));
}

double calculateMedian(std::map<std::string, int>& data) {
    std::vector<int> lengths;
    for (const auto& pair : data) {
        lengths.push_back(pair.second);
    }
    std::sort(lengths.begin(), lengths.end());
    int n = lengths.size();
    return n % 2 == 0 ? (lengths[n / 2 - 1] + lengths[n / 2]) / 2.0 : lengths[n / 2];
}

double calculateQuartile(std::map<std::string, int>& data, int quartilePercent) {
    std::vector<int> lengths;
    for (const auto& pair : data) {
        lengths.push_back(pair.second);
    }

    if (lengths.empty()) {
        return 0.0; 
    }

    std::sort(lengths.begin(), lengths.end());

    if (quartilePercent == 1) {
        int index = static_cast<int>((0.25 * lengths.size()) - 1);
        return lengths[index];
    } else if (quartilePercent == 2) {
        int index = static_cast<int>((0.50 * lengths.size()) - 1);
        return lengths[index];
    } else if (quartilePercent == 3) {
        int index = static_cast<int>((0.75 * lengths.size()) - 1);
        return lengths[index];
    } else if (quartilePercent == 4) {
        return static_cast<double>(*std::max_element(lengths.begin(), lengths.end()));
    }

    return 0.0; 
}


std::vector<int> multinomial(int n, const std::map<std::string, double>& abundanceData, unsigned int seed) {
    std::mt19937 gen(seed); 
    std::vector<double> pvals;

    for (const auto& pair : abundanceData) {
        pvals.push_back(pair.second);
    }

    std::discrete_distribution<> d(pvals.begin(), pvals.end());

    std::vector<int> counts(pvals.size(), 0);
    for (int i = 0; i < n; ++i) {
        int number = d(gen);
        counts[number]++;
    }
    return counts;
}


int main(int argc, char* argv[]) {
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("abundance_file", po::value<std::string>()->required(), "Path to abundance file")
            ("optimal_length_mode", po::value<std::string>()->required(), "Mode for optimal length")
            ("k_value", po::value<float>(), "k value for calculations")
            ("fasta_lengths", po::value<std::string>(), "Path to fasta lengths file")
            ("seed", po::value<int>()->required(), "Seed value")
            ("read_number", po::value<int>()->required(), "Read number value")
            ("output", po::value<std::string>()->required(), "Output file path");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if(vm.count("optimal_length_mode")) {
            std::string mode = vm["optimal_length_mode"].as<std::string>();
            if (mode != "NO") {
                desc.add_options()
                    ("k_value", po::value<float>()->required())
                    ("fasta_lengths", po::value<std::string>()->required());
            }
        }

        po::notify(vm);

        bool useFixedLength = false;
        int fixedLength = 0;
        bool useMedian = false;
        bool useQuartile = false;
        int quartileValue = 0;

        std::string optimalLengthMode = vm["optimal_length_mode"].as<std::string>();

        if (optimalLengthMode.find("fixed") != std::string::npos) {
            useFixedLength = true;
            std::istringstream iss(optimalLengthMode);
            std::string token;
            iss >> token; 
            if (iss >> fixedLength) {
                std::cout << "Using Fixed Length Mode with Length: " << fixedLength << "\n";
            } else {
                throw std::runtime_error("Fixed length not provided");
            }
        } else if (optimalLengthMode == "default") {
            useFixedLength = true;
            fixedLength = 350;
            std::cout << "Using Default Fixed Length Mode with Length: 350\n";
        } else if (optimalLengthMode == "median") {
            useMedian = true;
            std::cout << "Using Median Length Mode\n";
        } else if (optimalLengthMode.find("quartile") != std::string::npos) {
            useQuartile = true;
            std::istringstream iss(optimalLengthMode);
            std::string token;
            iss >> token; 
            if (iss >> quartileValue) {
                std::cout << "Using Quartile Length Mode with Quartile: " << quartileValue << "\n";
            } else {
                throw std::runtime_error("Quartile value not provided");
            }
        } else {
            std::cout << "Not using Optimal Length Mode\n";
        }

        auto abundanceData = readAbundanceFile(vm["abundance_file"].as<std::string>());
        std::map<std::string, int> fastaLengths;
        if (optimalLengthMode != "NO") {
            fastaLengths = readFastaLengthsFile(vm["fasta_lengths"].as<std::string>());
        }

        double midpoint = 0.0;
        if (useFixedLength) {
            midpoint = static_cast<double>(fixedLength);
        } else if (useMedian) {
            midpoint = calculateMedian(fastaLengths);
        } else if (useQuartile) {
            midpoint = calculateQuartile(fastaLengths, quartileValue);
        }

        std::map<std::string, double> p_value_lengths;
        if (!fastaLengths.empty() && optimalLengthMode != "NO") {
            for (const auto& pair : fastaLengths) {
                p_value_lengths[pair.first] = invertedSigmoid(pair.second, vm["k_value"].as<float>(), midpoint);
            }
            fastaLengths.clear(); 
        }

        if (optimalLengthMode != "NO") {
            for (auto& pair : abundanceData) {
                if (p_value_lengths.find(pair.first) != p_value_lengths.end()) {
                    pair.second *= p_value_lengths[pair.first];
                }
            }
        }

        double totalAbundance = std::accumulate(abundanceData.begin(), abundanceData.end(), 0.0, 
            [](double sum, const std::pair<std::string, double>& p) { return sum + p.second; });

        if (totalAbundance != 1.0) {
            for (auto& pair : abundanceData) {
                pair.second /= totalAbundance;
            }
        }

        int totalReadPairs = vm["read_number"].as<int>() / 2;

        auto seed = vm["seed"].as<int>();
        auto readPairCountsVector = multinomial(totalReadPairs, abundanceData, seed);

        std::map<std::string, int> readPairCounts;
        auto it = abundanceData.begin();
        for (int count : readPairCountsVector) {
            if (count > 0) { // Filtering out fragments with 0 read pairs
                readPairCounts[it->first] = count;
            }
            ++it;
        }

        std::string outputFile = vm["output"].as<std::string>();

        std::ofstream outStream(outputFile);

        if (!outStream.is_open()) {
            throw std::runtime_error("Cannot open output file: " + outputFile);
        }

        for (const auto& pair : readPairCounts) {
            outStream << pair.first << ";" << pair.second << "\n";
        }

        outStream.close();

    } catch(std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
