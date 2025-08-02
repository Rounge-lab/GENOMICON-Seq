#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <zlib.h>
#include <string>
#include <iterator>

std::vector<std::string> readGzFile(const std::string& filename) {
    gzFile gz = gzopen(filename.c_str(), "rb");
    if (!gz) {
        std::cerr << "Cannot open gzipped file" << std::endl;
        return {};
    }

    std::vector<std::string> lines;
    char buffer[128];
    std::string currentLine;

    while (int bytesRead = gzread(gz, buffer, sizeof(buffer) - 1)) {
        buffer[bytesRead] = '\0';
        currentLine += buffer;

        size_t newlinePos;
        while ((newlinePos = currentLine.find('\n')) != std::string::npos) {
            lines.push_back(currentLine.substr(0, newlinePos));
            currentLine = currentLine.substr(newlinePos + 1);
        }
    }

    gzclose(gz);
    return lines;
}

void writeGzFile(const std::vector<std::string>& lines, const std::string& filename) {
    gzFile gz = gzopen(filename.c_str(), "wb");
    if (!gz) {
        std::cerr << "Cannot open gzipped file for writing" << std::endl;
        return;
    }

    for (const auto& line : lines) {
        gzputs(gz, line.c_str());
        gzputs(gz, "\n");
    }

    gzclose(gz);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <input.bed.gz> <output_base_name> <number_of_splits> <seed>" << std::endl;
        return 1;
    }

    std::string inputFilename = argv[1];
    std::string outputBaseName = argv[2];
    int numberOfSplits = std::stoi(argv[3]);
    int seed = std::stoi(argv[4]);

    std::vector<std::string> lines = readGzFile(inputFilename);

    std::unordered_map<std::string, std::vector<std::string>> fragmentLines;
    for (const auto& line : lines) {
        std::istringstream iss(line);
        std::vector<std::string> tokens(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
        if (tokens.size() >= 4) {
            fragmentLines[tokens[3]].push_back(line);
        }
    }

    std::vector<std::string> uniqueFragments;
    for (const auto& pair : fragmentLines) {
        uniqueFragments.push_back(pair.first);
    }

    std::shuffle(uniqueFragments.begin(), uniqueFragments.end(), std::default_random_engine(seed));

    int maxFragmentsPerSplit = std::ceil(static_cast<double>(uniqueFragments.size()) / numberOfSplits);

    std::vector<std::vector<std::string>> splits(numberOfSplits);

    int currentSplit = 0;
    int fragmentsInCurrentSplit = 0;
    for (const auto& fragment : uniqueFragments) {
        splits[currentSplit].insert(splits[currentSplit].end(), fragmentLines[fragment].begin(), fragmentLines[fragment].end());
        fragmentsInCurrentSplit++;
        if (fragmentsInCurrentSplit >= maxFragmentsPerSplit) {
            currentSplit++;
            fragmentsInCurrentSplit = 0;
        }
    }

    for (int i = 0; i < numberOfSplits; ++i) {
        std::string outputFilename = outputBaseName + "_" + std::to_string(i+1) + ".bed.gz";
        writeGzFile(splits[i], outputFilename);
    }

    return 0;
}