#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <string>
#include <sstream>
#include <set>
#include <unordered_set>



using Identifier = std::vector<int>;

struct VectorHash {
    size_t operator()(const Identifier& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

using FragmentMap = std::unordered_map<Identifier, int, VectorHash>;

void readInputTable(const std::string& filePath, std::unordered_map<std::string, int>& fragments) {
    std::ifstream file(filePath);
    std::string line;

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fragmentName, temp, lengths;
        std::getline(iss, fragmentName, ';');
        std::getline(iss, temp, ';'); 
        std::getline(iss, lengths, ';');
        fragments[fragmentName] = std::stoi(lengths); 
    }
}


void recalibrateIdentifiers(const FragmentMap& sourceMap, FragmentMap& tempMap, int length) {
    tempMap.clear(); 

    for (const auto& pair : sourceMap) {
        if (pair.first.size() == 1 && pair.first[0] == 0) {
            tempMap[pair.first] = pair.second;
            continue;
        }

        Identifier recalibratedIdentifier;

        for (int pos : pair.first) {
            recalibratedIdentifier.push_back(length + 1 - pos);
        }

        std::sort(recalibratedIdentifier.begin(), recalibratedIdentifier.end());

        tempMap[recalibratedIdentifier] = pair.second;
    }
}

double calculateEfficiency(int cycle, int midpointCycle, float kParameter) {
    return 1.0 / (1.0 + exp(kParameter * (cycle - midpointCycle)));
}



void applyBinomial(FragmentMap& map, double efficiency, std::minstd_rand& generator) {
    for (auto& pair : map) {
        std::binomial_distribution<int> distribution(pair.second, efficiency);
        int successfulCopies = distribution(generator);

        pair.second = successfulCopies;
    }

    for (auto it = map.begin(); it != map.end();) {
        if (it->second == 0) {
            it = map.erase(it);
        } else {
            ++it;
        }
    }
}

std::vector<int> poissonDistribution(int copyNumber, double lambda, std::minstd_rand& generator) {
    std::vector<int> mutations(copyNumber, 0);
    std::poisson_distribution<int> distribution(lambda);

    for (int& mutation : mutations) {
        mutation = distribution(generator);
    }

    return mutations;
}


std::vector<int> pickMutationPositions(int numMutations, int length, const Identifier& excludedPositions, std::minstd_rand& generator) {
    std::unordered_set<int> positionsSet;
    std::uniform_int_distribution<int> distribution(1, length);

    while (positionsSet.size() < static_cast<size_t>(numMutations)) {
        int pos = distribution(generator);
        if (std::find(excludedPositions.begin(), excludedPositions.end(), pos) == excludedPositions.end() &&
            positionsSet.find(pos) == positionsSet.end()) {
            positionsSet.insert(pos);
        }
    }

    std::vector<int> positions(positionsSet.begin(), positionsSet.end());
    std::sort(positions.begin(), positions.end());
    return positions;
}




void simulateMutations(FragmentMap& mapTemp, double polymeraseErrorRate, int fragmentLength, std::minstd_rand& generator) {
    FragmentMap newMutationsMap;
    double lambda = polymeraseErrorRate * fragmentLength;
    
    for (const auto& pair : mapTemp) {
        const auto& identifier = pair.first;
        int copyNumber = pair.second;
        std::vector<int> mutationResults = poissonDistribution(copyNumber, lambda, generator);

        int noMutationCount = std::count(mutationResults.begin(), mutationResults.end(), 0);
        mapTemp[identifier] = noMutationCount;

        // Process mutations
        for (int mutationCount : mutationResults) {
            if (mutationCount > 0) {
                auto newPositions = pickMutationPositions(mutationCount, fragmentLength, identifier, generator);
                
                Identifier newIdentifier;
                if (identifier == Identifier{0}) {
                    newIdentifier = newPositions; 
                } else {
                    newIdentifier = identifier;
                    newIdentifier.insert(newIdentifier.end(), newPositions.begin(), newPositions.end());
                    std::sort(newIdentifier.begin(), newIdentifier.end()); 
                }

                if (newMutationsMap.find(newIdentifier) != newMutationsMap.end()) {
                    newMutationsMap[newIdentifier] += 1;
                } else {
                    newMutationsMap[newIdentifier] = 1;
                }
            }
        }
    }

    for (const auto& newMutation : newMutationsMap) {
        if (mapTemp.find(newMutation.first) != mapTemp.end()) {
            mapTemp[newMutation.first] += newMutation.second;
        } else {
            mapTemp.insert({newMutation.first, newMutation.second});
        }
    }
}


void mergeMaps(FragmentMap& originalMap, const FragmentMap& tempMap) {
    for (const auto& pair : tempMap) {
        auto& identifier = pair.first;
        auto copyNumber = pair.second;

        if (originalMap.find(identifier) != originalMap.end()) {
            originalMap[identifier] += copyNumber;
        } else {
            originalMap.insert(pair);
        }
    }
}

std::string vectorToString(const std::vector<int>& v) {
    std::string result;
    for (auto i : v) {
        if (!result.empty()) result += ",";
        result += std::to_string(i);
    }
    return result.empty() ? "0" : result;
}

unsigned int deriveSeedForFragment(const std::string& fragmentName, unsigned int masterSeed) {
    std::hash<std::string> hasher;
    size_t hash = hasher(fragmentName);
    return static_cast<unsigned int>(hash ^ masterSeed);
}


void preGenerateFragmentSeeds(const std::unordered_map<std::string, int>& fragments, unsigned int masterSeed, std::unordered_map<std::string, unsigned int>& fragmentSeeds) {
    for (const auto& fragment : fragments) {
        fragmentSeeds[fragment.first] = deriveSeedForFragment(fragment.first, masterSeed);
    }
}


void simulatePCR(const std::unordered_map<std::string, int>& fragments, int numCycles, int midpointCycle, float kParameter, double polymeraseErrorRate, const std::unordered_map<std::string, unsigned int>& fragmentSeeds, const std::string& outputPath, float fracToSequencing) {
    int fragmentCount = 0;
    int fileCount = 0; 
    std::ofstream outputFile;
    std::ostringstream outputBuffer; 

    for (const auto& fragment : fragments) {
        if (fragmentCount % 100 == 0) {
            if (outputFile.is_open()) {
                outputFile << outputBuffer.str(); 
                outputFile.close(); 
                outputBuffer.str(""); 
                outputBuffer.clear(); 
            }
            outputFile.open(outputPath + "_" + std::to_string(++fileCount) + ".csv"); 
        }
        
        unsigned int fragmentSeed = fragmentSeeds.at(fragment.first);
        std::minstd_rand fragmentGenerator(fragmentSeed);
        
        FragmentMap F_map, R_map;
        F_map[{0}] = 1; 
        R_map[{0}] = 1;

        for (int cycle = 1; cycle <= numCycles; ++cycle) {
            float currentEfficiency = calculateEfficiency(cycle, midpointCycle, kParameter);
            FragmentMap F_map_temp, R_map_temp;
            
            recalibrateIdentifiers(R_map, F_map_temp, fragments.at(fragment.first));
            recalibrateIdentifiers(F_map, R_map_temp, fragments.at(fragment.first));
            applyBinomial(F_map_temp, currentEfficiency, fragmentGenerator);
            applyBinomial(R_map_temp, currentEfficiency, fragmentGenerator);
            simulateMutations(F_map_temp, polymeraseErrorRate, fragments.at(fragment.first), fragmentGenerator);
            simulateMutations(R_map_temp, polymeraseErrorRate, fragments.at(fragment.first), fragmentGenerator);
            mergeMaps(F_map, F_map_temp);
            mergeMaps(R_map, R_map_temp);
        }

        applyBinomial(F_map, fracToSequencing, fragmentGenerator);

        for (const auto& mapPair : F_map) {
            std::string identifierStr = (mapPair.first.empty() || mapPair.first[0] == 0) ? "No mutations" : vectorToString(mapPair.first);
            outputBuffer << fragment.first << ";1;" << fragments.at(fragment.first) << ";5-3;" << identifierStr << ";" << mapPair.second << "\n";
        }
        fragmentCount++;
    }
    if (outputFile.is_open()) {
        outputFile << outputBuffer.str(); 
        outputFile.close(); 
    }
}




int main(int argc, char* argv[]) {
    if (argc < 9) {
        std::cerr << "Usage: " << argv[0] << " <inputTablePath> <seed> <polymeraseErrorRate> <numCycles> <midpointCycle> <kParameter> <fracToSequencing> <outputPath>" << std::endl;
        return 1;
    }

    std::string inputTablePath = argv[1];
    int seed = std::stoi(argv[2]);
    double polymeraseErrorRate = std::stod(argv[3]);
    int numCycles = std::stoi(argv[4]);
    int midpointCycle = std::stoi(argv[5]);
    float kParameter = std::stof(argv[6]);
    float fracToSequencing = std::stof(argv[7]);
    std::string outputPath = argv[8]; 

    std::unordered_map<std::string, int> fragments;
    readInputTable(inputTablePath, fragments);

    std::unordered_map<std::string, unsigned int> fragmentSeeds;
    preGenerateFragmentSeeds(fragments, static_cast<unsigned int>(seed), fragmentSeeds);

    simulatePCR(fragments, numCycles, midpointCycle, kParameter, polymeraseErrorRate, fragmentSeeds, outputPath, fracToSequencing);

    return 0;
}

