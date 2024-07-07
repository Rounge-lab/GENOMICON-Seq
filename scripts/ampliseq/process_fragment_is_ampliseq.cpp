#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <sqlite3.h>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <string>

// Function declarations
std::vector<std::string> parseCsvLine(const std::string& line);
std::vector<std::string> split(const std::string &s, char delimiter);
std::string executeQuery(sqlite3* db, const std::string& query);
std::string queryDatabase(const std::string& databasePath, const std::string& chrName, int start, int end, int genomeLength);
void applyMutations(std::string& sequence, const std::string& mutations, int fragmentStart, int fragmentEnd, int genomeLength, const std::string& chrName, std::vector<std::string>& realMutationTracker, std::map<std::string, int>& mutationCount);
void mutatePcrSequence(std::string& sequence, const std::string& pcrMutations, std::mt19937& gen, std::map<std::string, int>& mutationTracking, int fragmentStart, int fragmentEnd, int genomeLength);
void writeFastaChunks(const std::string& outputBasePath, const std::string& fastaContent, size_t& fileCounter, std::string& buffer, size_t& currentSize);
std::string processFragment(const std::string &fragmentRow, const std::map<std::string, std::string> &sampleTablePaths, const std::string &sqlDatabasePath, std::string &lastChrName, std::ifstream &sampleTableFile, std::mt19937& gen, const std::unordered_map<std::string, int>& genomeLengths, std::map<std::string, int>& mutationTracking, std::vector<std::string>& realMutationTracker, std::map<std::string, int>& mutationCount);
int main(int argc, char* argv[]);

// Function definitions
std::vector<std::string> parseCsvLine(const std::string& line) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream stream(line);
    char delimiter = ',';
    char quote = '"';

    while (stream.good()) {
        getline(stream, token, delimiter);
        token.erase(0, token.find_first_not_of(' '));
        token.erase(token.find_last_not_of(' ') + 1);
        if (!token.empty() && token.front() == quote) {
            if (token.back() == quote) {
                token = token.substr(1, token.size() - 2);
            } else {
                std::string nextPart;
                while (getline(stream, nextPart, delimiter)) {
                    token += ',' + nextPart;
                    if (!nextPart.empty() && nextPart.back() == quote) {
                        token = token.substr(1, token.size() - 2);
                        break;
                    }
                }
            }
        }
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string executeQuery(sqlite3* db, const std::string& query) {
    sqlite3_stmt* stmt;
    std::string sequence;
    int rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
    if (rc == SQLITE_OK) {
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            sequence += reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        }
    } else {
        std::cerr << "SQL error during query: " << sqlite3_errmsg(db) << std::endl;
    }
    sqlite3_finalize(stmt);
    return sequence;
}

std::string queryDatabase(const std::string& databasePath, const std::string& chrName, int start, int end, int genomeLength) {
    sqlite3* db;
    int rc = sqlite3_open(databasePath.c_str(), &db);
    if (rc) {
        std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
        return "";
    }

    std::string query;
    std::string nonmutatedSequence;

    if (start > end && genomeLength > 0) {
        // For circular logic
        query = "SELECT nucleotide FROM " + chrName + " WHERE rowid BETWEEN " + std::to_string(start) + " AND " + std::to_string(genomeLength) + ";";
        nonmutatedSequence += executeQuery(db, query);

        std::string querySecondPart = "SELECT nucleotide FROM " + chrName + " WHERE rowid BETWEEN 1 AND " + std::to_string(end) + ";";
        nonmutatedSequence += executeQuery(db, querySecondPart);
    } else {
        // For linear logic
        query = "SELECT nucleotide FROM " + chrName + " WHERE rowid BETWEEN " + std::to_string(start) + " AND " + std::to_string(end) + ";";
        nonmutatedSequence = executeQuery(db, query);
    }

    sqlite3_close(db);

    return nonmutatedSequence;
}

void applyMutations(std::string& sequence, const std::string& mutations, int fragmentStart, int fragmentEnd, int genomeLength, const std::string& chrName, std::vector<std::string>& realMutationTracker, std::map<std::string, int>& mutationCount) {
    if (mutations == "No mutation") {
        return;
    }

    std::vector<std::string> mutationList = split(mutations, ',');
    for (const std::string& mut : mutationList) {
        int globalMutPosition = std::stoi(mut.substr(0, mut.size() - 1));
        char mutNucleotide = mut.back();

        bool withinRange = false;

        // Check if the mutation is within the range
        if (fragmentStart <= fragmentEnd) {
            // Linear genome
            withinRange = (globalMutPosition >= fragmentStart && globalMutPosition <= fragmentEnd);
        } else {
            // Circular genome
            withinRange = (globalMutPosition >= fragmentStart || globalMutPosition <= fragmentEnd);
        }

        if (withinRange) {
            // Convert global position to local position relative to fragment start
            int localMutPosition = globalMutPosition - fragmentStart;
            if (localMutPosition < 0) {
                localMutPosition += genomeLength;
            }

            // Ensure local position is within the fragment sequence length
            localMutPosition %= sequence.length();
            if (localMutPosition >= 0 && localMutPosition < sequence.length()) {
                // Apply mutation to the sequence
                sequence[localMutPosition] = mutNucleotide;
                // Record the mutation only if it is applied
                std::string mutationKey = chrName + "," + std::to_string(globalMutPosition) + mutNucleotide;
                realMutationTracker.push_back(mutationKey);
                mutationCount[mutationKey]++;
            }
        }
    }
}

void mutatePcrSequence(std::string& sequence, const std::string& pcrMutations, std::mt19937& gen, std::map<std::string, int>& mutationTracking, int fragmentStart, int fragmentEnd, int genomeLength) {
    if (pcrMutations == "No_mutations") {
        return;
    }

    std::vector<std::string> mutationPositions = split(pcrMutations, ',');
    std::unordered_map<char, std::vector<char>> nucleotideReplacements{
        {'A', {'C', 'G', 'T'}},
        {'C', {'A', 'G', 'T'}},
        {'G', {'A', 'C', 'T'}},
        {'T', {'A', 'C', 'G'}}
    };

    for (const std::string& posStr : mutationPositions) {
        int localPosition = std::stoi(posStr) - 1; // Adjusting PCR position from 1-based to 0-based
        if (localPosition < 0 || localPosition >= sequence.length()) {
            continue;
        }

        char originalNucleotide = sequence[localPosition];
        const auto& replacements = nucleotideReplacements[originalNucleotide];

        std::uniform_int_distribution<> distr(0, replacements.size() - 1);
        char newNucleotide = replacements[distr(gen)];

        sequence[localPosition] = newNucleotide; // Apply mutation

        // Calculate the genomic position of the mutation
        int genomicPosition;
        if (fragmentStart > fragmentEnd) { // Fragment spans the genome's origin
            int adjustedPosition = localPosition + fragmentStart;
            if (adjustedPosition > genomeLength) {
                adjustedPosition = adjustedPosition % genomeLength;
            }
            genomicPosition = adjustedPosition;
        } else {
            genomicPosition = localPosition + fragmentStart;
        }

        // Build a unique key for the mutation (genomic position and new nucleotide) and update the tracking map
        std::string mutationKey = std::to_string(genomicPosition) + newNucleotide;
        mutationTracking[mutationKey]++;
    }
}

void writeFastaChunks(const std::string& outputBasePath, const std::string& fastaContent, size_t& fileCounter, std::string& buffer, size_t& currentSize) {
    size_t chunkSizeLimit = 10485760; // 1MB in bytes
    buffer += fastaContent;
    currentSize += fastaContent.size();

    if (currentSize >= chunkSizeLimit) {
        std::ofstream outFile(outputBasePath + "-" + std::to_string(fileCounter++) + ".fasta");
        outFile << buffer;
        outFile.close();

        buffer.clear();
        currentSize = 0;
    }
}

std::string processFragment(const std::string &fragmentRow, const std::map<std::string, std::string> &sampleTablePaths, const std::string &sqlDatabasePath, std::string &lastChrName, std::ifstream &sampleTableFile, std::mt19937& gen, const std::unordered_map<std::string, int>& genomeLengths, std::map<std::string, int>& mutationTracking, std::vector<std::string>& realMutationTracker, std::map<std::string, int>& mutationCount) {
    std::vector<std::string> row = split(fragmentRow, ';');
    std::string fragmentName = row[0];
    std::string pcrMutations = row[1];
    std::string fragmentCoordinates = row[2];

    // Extract chr_name and original_genome_name
    std::vector<std::string> nameParts = split(fragmentName, '_');
    std::string chrName = nameParts[0];
    std::string originalGenomeName = chrName + "_" + nameParts[1];

    int genomeLength = 0;

    if (chrName != lastChrName) {
        lastChrName = chrName;
        sampleTableFile.close();
        sampleTableFile.open(sampleTablePaths.at(chrName));

        auto it = genomeLengths.find(chrName);
        if (it != genomeLengths.end()) {
            genomeLength = it->second;
        } else {
            it = genomeLengths.find(chrName);
            if (it != genomeLengths.end()) {
                genomeLength = it->second;
            } else {
                std::cerr << "Failed to fetch a valid genome length for " << chrName << ". Aborting fragment processing." << std::endl;
                return "";
            }
        }
    } else {
        if (genomeLength <= 0) {
            auto lengthIt = genomeLengths.find(chrName);
            if (lengthIt != genomeLengths.end()) {
                genomeLength = lengthIt->second;
            } else {
                std::cerr << "Unable to correct invalid genome length for " << chrName << ". Aborting fragment processing." << std::endl;
                return "";
            }
        }
    }

    std::string originalMutations = "No mutation";
    if (sampleTableFile.is_open()) {
        std::string sampleLine;
        std::getline(sampleTableFile, sampleLine);

        while (std::getline(sampleTableFile, sampleLine)) {
            std::vector<std::string> sampleRow = parseCsvLine(sampleLine);

            std::string genomeNameInTable = sampleRow[0];
            if (genomeNameInTable.front() == '"' && genomeNameInTable.back() == '"') {
                genomeNameInTable = genomeNameInTable.substr(1, genomeNameInTable.length() - 2);
            }

            if (genomeNameInTable == originalGenomeName) {
                if (sampleRow[1].front() == '"' && sampleRow[1].back() == '"') {
                    originalMutations = sampleRow[1].substr(1, sampleRow[1].length() - 2);
                } else {
                    originalMutations = sampleRow[1];
                }
                break;
            }
        }
        sampleTableFile.clear();
        sampleTableFile.seekg(0, std::ios::beg);
    }

    std::vector<std::string> coordinates = split(fragmentCoordinates, ':');

    int fragmentStart = 0;
    int fragmentEnd = 0;
    if (coordinates.size() >= 2) {
        try {
            fragmentStart = std::stoi(coordinates[0]);
            fragmentEnd = std::stoi(coordinates[1]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid argument. Unable to convert '" << coordinates[0] << "' or '" << coordinates[1] << "' to integers." << std::endl;
            return "";
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Out of range. The values '" << coordinates[0] << "' or '" << coordinates[1] << "' are out of integer range." << std::endl;
            return "";
        }
    } else {
        std::cerr << "Error: Incorrect format for fragment coordinates. Expected format 'start:end', got '" << fragmentCoordinates << "'." << std::endl;
        return "";
    }

    std::string databasePath = sqlDatabasePath + "/" + chrName + ".sqlite";
    std::string nonmutatedSequence = queryDatabase(databasePath, chrName, fragmentStart, fragmentEnd, genomeLength);

    applyMutations(nonmutatedSequence, originalMutations, fragmentStart, fragmentEnd, genomeLength, chrName, realMutationTracker, mutationCount);
    mutatePcrSequence(nonmutatedSequence, pcrMutations, gen, mutationTracking, fragmentStart, fragmentEnd, genomeLength);

    std::replace(pcrMutations.begin(), pcrMutations.end(), ' ', '_');
    std::replace(pcrMutations.begin(), pcrMutations.end(), ',', '_');

    return ">" + fragmentName + "_" + pcrMutations + "\n" + nonmutatedSequence + "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <pcr_table_temp> <sample_table_base_path> <SQLdatabase_path> <seed> <genome_lengths_table> <output_base_path>" << std::endl;
        return 1;
    }

    std::string pcrTablePath = argv[1];
    std::string sampleTableBasePath = argv[2];
    std::string sqlDatabasePath = argv[3];
    int seed = std::stoi(argv[4]);
    std::string genomeLengthsTable = argv[5];
    std::string outputBasePath = argv[6];

    std::map<std::string, int> mutationTracking;
    std::map<std::string, int> mutationCount; // New map to track mutation counts

    std::mt19937 gen(seed);
    std::ifstream pcrTableFile(pcrTablePath);
    if (!pcrTableFile.is_open()) {
        std::cerr << "Failed to open file: " << pcrTablePath << std::endl;
        return 1;
    }
    std::unordered_map<std::string, int> genomeLengths;
    std::ifstream genomeLengthsFile(genomeLengthsTable);
    if (!genomeLengthsFile.is_open()) {
        std::cerr << "Failed to open genome lengths file: " << genomeLengthsTable << std::endl;
        return 1;
    }

    std::string line;
    std::getline(genomeLengthsFile, line);

    while (std::getline(genomeLengthsFile, line)) {
        std::istringstream ss(line);
        std::string fastaName, lengthStr;
        std::getline(ss, fastaName, ',');
        std::getline(ss, lengthStr, ',');

        if (fastaName.front() == '"' && fastaName.back() == '"') {
            fastaName = fastaName.substr(1, fastaName.size() - 2);
        }

        try {
            int length = std::stoi(lengthStr);
            genomeLengths[fastaName] = length;
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error converting length to integer for " << fastaName << " with value: " << lengthStr << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Length value out of range for " << fastaName << " with value: " << lengthStr << std::endl;
        }
    }

    std::set<std::string> uniqueChromosomes;
    while (std::getline(pcrTableFile, line)) {
        std::stringstream ss(line);
        std::string fragmentName;
        std::getline(ss, fragmentName, ';');

        size_t underscorePos = fragmentName.find('_');
        if (underscorePos != std::string::npos) {
            uniqueChromosomes.insert(fragmentName.substr(0, underscorePos));
        }
    }

    std::map<std::string, std::string> sampleTablePaths;
    for (const std::string& chromosome : uniqueChromosomes) {
        sampleTablePaths[chromosome] = sampleTableBasePath + "/" + chromosome + "_sample_table.csv";
    }

    pcrTableFile.clear();
    pcrTableFile.seekg(0);

    std::string lastChrName = "";
    std::ifstream sampleTableFile;

    size_t fileCounter = 1;
    std::string buffer;
    size_t currentSize = 0;

    std::vector<std::string> realMutationTracker;

    while (std::getline(pcrTableFile, line)) {
        std::string fastaSequence = processFragment(line, sampleTablePaths, sqlDatabasePath, lastChrName, sampleTableFile, gen, genomeLengths, mutationTracking, realMutationTracker, mutationCount);
        writeFastaChunks(outputBasePath, fastaSequence, fileCounter, buffer, currentSize);
    }

    if (!buffer.empty()) {
        std::ofstream outFile(outputBasePath + "-" + std::to_string(fileCounter) + ".fasta");
        outFile << buffer;
        outFile.close();
    }

    if (!mutationTracking.empty()) {
        // Only proceed with file writing if there are mutations to track
        std::ofstream mutationTrackingFile(outputBasePath + "-mutation_tracking.csv");
        mutationTrackingFile << "pcr_mutations,copy_number\n"; // Header row
        for (const auto& pair : mutationTracking) {
            // Extract position and nucleotide from the key
            std::string posAndNuc = pair.first;
            int count = pair.second;

            // Write to file
            mutationTrackingFile << posAndNuc << "," << count << "\n";
        }
        mutationTrackingFile.close();
    }

    if (!realMutationTracker.empty()) {
        std::ofstream realMutationFile(outputBasePath + "-real_mutations.csv");
        realMutationFile << "chrName,inserted_mutation,count\n"; // Header row
        for (const auto& pair : mutationCount) {
            realMutationFile << pair.first << "," << pair.second << "\n";
        }
        realMutationFile.close();
    }

    pcrTableFile.close();
    sampleTableFile.close();

    return 0;
}
