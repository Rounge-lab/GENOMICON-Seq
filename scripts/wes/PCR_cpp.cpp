#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <thread>
#include <random>
#include <functional>
#include <mutex>
#include <cmath> 

void applyBinomial(std::vector<int>& copy_number, double efficiency, std::mt19937& gen) {
    for (auto& cn : copy_number) {
        std::binomial_distribution<> d(cn, efficiency);
        cn += d(gen);
    }
}

unsigned int deriveSeedForFragment(const std::string& fragmentName, unsigned int masterSeed) {
    std::hash<std::string> hasher;
    size_t hash = hasher(fragmentName);
    return static_cast<unsigned int>(hash ^ masterSeed);
}

void processBatch(const std::vector<std::string>& batch, int threadId, const std::string& outputDir, int numberOfPcrCycles, int midpointCycles, double kParameter, long mainSeed) {
    std::ostringstream outputFilePath;
    outputFilePath << outputDir << "/worker_" << threadId << "_temp.csv";

    std::ofstream outputFile(outputFilePath.str());
    outputFile << "Names;Mutation_Positions;Copies;genome_coordinates\n";

    for (const auto& row : batch) {
        std::stringstream ss(row);
        std::string token;
        std::vector<std::string> rowData;
        while (getline(ss, token, ';')) {
            rowData.push_back(token);
        }

        std::string name = rowData[0];
        std::string coordinates = rowData[1];
        int length = std::stoi(rowData[2]);
        int copy_number = 1;
        std::string mutation_positions = "No_mutations";

        unsigned int fragmentSeed = deriveSeedForFragment(name, static_cast<unsigned int>(mainSeed));
        std::mt19937 gen(fragmentSeed);

        std::vector<int> copy_numbers(numberOfPcrCycles, copy_number);
        for (int cycle = 1; cycle <= numberOfPcrCycles; ++cycle) {
            double efficiency = 1 / (1 + exp(kParameter * (cycle - midpointCycles)));
            applyBinomial(copy_numbers, efficiency, gen);
        }

        outputFile << name << ";" << mutation_positions << ";" << copy_numbers.back() << ";" << coordinates << "\n";
    }

    outputFile.close();
}


void processBatch(const std::vector<std::string>& batch, int threadId, const std::string& outputDir, int numberOfPcrCycles, int midpointCycles, double kParameter, long mainSeed);

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <input_file_path> <number_of_cpus> <seed> <number_of_pcr_cycles> <midpoint_cycles> <k_parameter> <output_dir>\n";
        return 1;
    }

    std::string inputFilePath = argv[1];
    int numberOfCpus = std::stoi(argv[2]);
    long mainSeed = std::stol(argv[3]);
    int numberOfPcrCycles = std::stoi(argv[4]);
    int midpointCycles = std::stoi(argv[5]);
    double kParameter = std::stod(argv[6]);
    std::string outputDir = argv[7];

    std::ifstream inputFile(inputFilePath);
    std::string line;
    std::vector<std::vector<std::string>> batches(numberOfCpus);
    int batchIndex = 0;

    while (getline(inputFile, line)) {
        batches[batchIndex % numberOfCpus].push_back(line);
        batchIndex++;
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < numberOfCpus; ++i) {
        threads.emplace_back(processBatch, std::ref(batches[i]), i, outputDir, numberOfPcrCycles, midpointCycles, kParameter, mainSeed);
    }

    for (auto& t : threads) {
        t.join();
    }

    return 0;
}

