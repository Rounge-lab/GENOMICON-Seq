#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <random>

struct Fragment {
    int start;
    int end;
};

void readCSV(const std::string& filePath, std::vector<std::string>& genomeNames, std::vector<int>& copyNumbers) {
    std::ifstream file(filePath);
    std::string line, name, temp;
    int copyNumber;

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::getline(ss, name, ';');
        std::getline(ss, temp, ';');
        ss >> copyNumber;

        if (!name.empty() && name[0] == '"') {
            name = name.substr(1, name.size() - 2);
        }

        genomeNames.push_back(name);
        copyNumbers.push_back(copyNumber);
    }
}




void expandNames(const std::vector<std::string>& genomeNames, std::vector<int>& copyNumbers, std::vector<std::string>& expandedNames) {
    for (size_t i = 0; i < genomeNames.size(); ++i) {
        for (int j = 1; j <= copyNumbers[i]; ++j) {
            expandedNames.push_back(genomeNames[i] + "_c" + std::to_string(j));
        }
    }
}

void distributeWork(const std::vector<std::string>& expandedNames, int numCores, std::vector<std::vector<std::string>>& batches) {
    size_t totalSize = expandedNames.size();
    size_t batchSize = totalSize / numCores;
    size_t remainingElements = totalSize % numCores;

    auto it = expandedNames.begin();
    for (int i = 0; i < numCores; ++i) {
        size_t currentBatchSize = batchSize + (i < remainingElements ? 1 : 0);
        std::vector<std::string> batch(it, it + currentBatchSize);
        batches.push_back(batch);
        it += currentBatchSize;
    }
}


void generateSeeds(const std::vector<std::string>& names, int globalSeed, std::vector<int>& seeds) {
    for (const auto& name : names) {
        int seed = std::hash<std::string>{}(name) ^ globalSeed;
        seeds.push_back(seed);
    }
}

void generateFragments(const std::string& name, int minLen, int maxLen, int genomeLen, int seed, std::vector<Fragment>& fragments) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(minLen, maxLen);
    std::uniform_int_distribution<> dis2(1, genomeLen - minLen + 1);

    while (genomeLen >= minLen) {
        int frag_length = dis(gen);
        frag_length = std::min(frag_length, genomeLen);
        int start_pos = dis2(gen);
        int end_pos = start_pos + frag_length - 1;

        fragments.push_back(Fragment{start_pos, end_pos});
        genomeLen -= frag_length;
    }
}

void writeOutput(const std::string& batchName, const std::vector<std::string>& names, std::vector<std::vector<Fragment>>& allFragments) {
    std::ofstream file(batchName);
    for (size_t i = 0; i < names.size(); ++i) {
        file << names[i] << ";";
        for (const auto& frag : allFragments[i]) {
            file << frag.start << ":" << frag.end << ",";
        }
        file << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " [CSV path] [fasta length] [min fragment length] [max fragment length] [number of cores] [seed] [output path]" << std::endl;
        return 1;
    }

    std::string csvPath = argv[1];
    int fastaLength = std::stoi(argv[2]);
    int minFragmentLength = std::stoi(argv[3]);
    int maxFragmentLength = std::stoi(argv[4]);
    int numCores = std::stoi(argv[5]);
    int seed = std::stoi(argv[6]);
    std::string outputPath = argv[7];

    std::vector<std::string> genomeNames, expandedNames;
    std::vector<int> copyNumbers;
    readCSV(csvPath, genomeNames, copyNumbers);
    expandNames(genomeNames, copyNumbers, expandedNames);

    std::vector<std::vector<std::string>> batches;
    distributeWork(expandedNames, numCores, batches);

    std::vector<std::thread> threads;
    for (size_t i = 0; i < batches.size(); ++i) {
        threads.push_back(std::thread([&, i]() {
            std::vector<int> seeds;
            generateSeeds(batches[i], seed, seeds);

            std::vector<std::vector<Fragment>> allFragments;
            for (size_t j = 0; j < batches[i].size(); ++j) {
                const auto& name = batches[i][j];
                std::vector<Fragment> fragments;
                int nameSeed = seeds[j];
                generateFragments(name, minFragmentLength, maxFragmentLength, fastaLength, nameSeed, fragments);
                allFragments.push_back(fragments);
            }

            std::string batchName = outputPath + "_batch" + std::to_string(i) + ".csv";
            writeOutput(batchName, batches[i], allFragments);
        }));
    }


    for (auto& thread : threads) {
        thread.join();
    }

    return 0;
}
