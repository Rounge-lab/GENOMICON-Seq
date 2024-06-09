#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <string>
#include <sqlite3.h>

std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::pair<int, char>> parseAndFilterMutations(const std::string& mutations, int fragmentStart, int fragmentEnd) {
    if (mutations == "No_mutations") return {};

    size_t estimatedMutationCount = std::count(mutations.begin(), mutations.end(), ',') + 1;
    std::vector<std::string> mutationList;
    mutationList.reserve(estimatedMutationCount);
    mutationList = split(mutations, ',');

    std::vector<std::pair<int, char>> filteredMutations;
    filteredMutations.reserve(50);
    for (const auto& mut : mutationList) {
        int position = std::stoi(mut.substr(0, mut.size() - 1));
        char nucleotide = mut.back();
        if (position >= fragmentStart && position <= fragmentEnd) {
            filteredMutations.emplace_back(position - fragmentStart, nucleotide);
        }
    }
    return filteredMutations;
}

std::string queryDatabase(sqlite3* db, const std::string& chrName, int start, int end) {
    std::string query = "SELECT nucleotide FROM " + chrName + " WHERE rowid BETWEEN " + std::to_string(start) + " AND " + std::to_string(end) + ";";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);

    std::string nonmutatedSequence;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        nonmutatedSequence += reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
    }

    sqlite3_finalize(stmt);
    return nonmutatedSequence;
}


void applyMutations(std::string& sequence, const std::vector<std::pair<int, char>>& mutations) {
    for (const auto& mutation : mutations) {
        int pos = mutation.first;
        char nucleotide = mutation.second;
        if (pos >= 0 && pos < sequence.length()) {
            sequence[pos] = nucleotide;
        }
    }
}

void mutatePcrSequence(std::string& sequence, const std::string& pcrMutations, std::mt19937& gen) {
    if (pcrMutations == "No_mutations") {
        return;
    }

    std::vector<std::string> mutationPositions = split(pcrMutations, ',');
    std::unordered_map<char, std::vector<char>> nucleotideReplacements = {
        {'A', {'C', 'G', 'T'}},
        {'C', {'A', 'G', 'T'}},
        {'G', {'A', 'C', 'T'}},
        {'T', {'A', 'C', 'G'}}
    };

    for (const std::string& posStr : mutationPositions) {
        int position = std::stoi(posStr) - 1;
        if (position < 0 || position >= sequence.length()) {
            continue;
        }

        char originalNucleotide = sequence[position];
        const auto& replacements = nucleotideReplacements[originalNucleotide];
        std::uniform_int_distribution<> distr(0, replacements.size() - 1);
        char newNucleotide = replacements[distr(gen)];
        sequence[position] = newNucleotide;
    }
}

void writeFastaChunks(const std::string& outputBasePath, const std::string& fastaContent, size_t& fileCounter, std::string& buffer, size_t& currentSize) {
    size_t chunkSizeLimit = 20971520;  // 20MB in bytes
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

std::string processFragment(const std::string &fragmentRow, const std::map<std::string, std::string> &genomeToMutations,
                            const std::string &sqlDatabasePath, std::string &lastChrName, sqlite3* &db, std::mt19937& gen) {
    // Split the input row into constituent parts
    std::vector<std::string> row = split(fragmentRow, ';');
    std::string fragmentName = row[0];
    std::string pcrMutations = row[1];
    std::string fragmentCoordinates = row[3];

    // Parse the coordinates for the fragment
    std::vector<std::string> coordinates = split(fragmentCoordinates, ':');
    int fragmentStart = std::stoi(coordinates[0]);
    int fragmentEnd = std::stoi(coordinates[1]);

    // Extract chromosome and genome information from the fragment name
    std::vector<std::string> nameParts = split(fragmentName, '_');
    std::string chrName = nameParts[0];
    std::string originalGenomeName = nameParts[0] + "_" + nameParts[1];

    // Manage database connection based on chromosome changes
    if (chrName != lastChrName) {
        if (db != nullptr) {
            sqlite3_close(db);  // Close the current connection
            db = nullptr;
        }
        std::string dbPath = sqlDatabasePath + "/" + chrName + ".sqlite";
        sqlite3_open(dbPath.c_str(), &db);  // Open a new connection for the current chromosome
        lastChrName = chrName;
    }

    // Query the database for non-mutated sequence
    std::string nonmutatedSequence = queryDatabase(db, chrName, fragmentStart, fragmentEnd);

    // Filter mutations that apply to this fragment
    std::vector<std::pair<int, char>> filteredMutations;
    if (genomeToMutations.count(originalGenomeName)) {
        filteredMutations = parseAndFilterMutations(genomeToMutations.at(originalGenomeName), fragmentStart, fragmentEnd);
    }

    // Apply any genomic and PCR-induced mutations
    applyMutations(nonmutatedSequence, filteredMutations);
    mutatePcrSequence(nonmutatedSequence, pcrMutations, gen);

    // Format mutations for output
    std::replace(pcrMutations.begin(), pcrMutations.end(), ' ', '_');
    std::replace(pcrMutations.begin(), pcrMutations.end(), ',', '_');

    // Return the formatted sequence
    return ">" + fragmentName + "_" + pcrMutations + "\n" + nonmutatedSequence + "\n";
}


int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <pcr_table_path> <sample_table_base_path> <SQL_database_path> <seed> <output_base_path>\n";
        return 1;
    }

    std::string pcrTablePath = argv[1];
    std::string sampleTableBasePath = argv[2];
    std::string sqlDatabasePath = argv[3];
    int seed = std::stoi(argv[4]);
    std::string outputBasePath = argv[5];

    std::mt19937 gen(seed);
    std::ifstream pcrTableFile(pcrTablePath);
    if (!pcrTableFile.is_open()) {
        std::cerr << "Failed to open PCR table file: " << pcrTablePath << "\n";
        return 1;
    }

    std::set<std::string> neededGenomes;
    std::string line;
    while (std::getline(pcrTableFile, line)) {
        std::vector<std::string> parts = split(line, ';');
        std::vector<std::string> nameParts = split(parts[0], '_');
        neededGenomes.insert(nameParts[0] + "_" + nameParts[1]);  // Collecting chrX_MGX
    }

    std::map<std::string, std::string> genomeToMutations;
    for (const auto& genome : neededGenomes) {
        std::string chrPart = genome.substr(0, genome.find('_'));
        std::string sampleFilePath = sampleTableBasePath + "/" + chrPart + "_sample_table.csv";
        std::ifstream sampleFile(sampleFilePath);
        if (!sampleFile.is_open()) {
            std::cerr << "Failed to open sample table file: " << sampleFilePath << "\n";
            continue;
        }
        std::string sampleLine;
        while (std::getline(sampleFile, sampleLine)) {
            std::vector<std::string> sampleParts = split(sampleLine, ';');
            if (genome == sampleParts[0]) {
                genomeToMutations[genome] = sampleParts[1];
            }
        }
    }

    pcrTableFile.clear();
    pcrTableFile.seekg(0);
    std::string buffer;
    size_t currentSize = 0;
    size_t fileCounter = 1;
    std::string lastChrName = "";
    sqlite3* db = nullptr;  // Initialize the database pointer to null

    while (std::getline(pcrTableFile, line)) {
        std::string fastaSequence = processFragment(line, genomeToMutations, sqlDatabasePath, lastChrName, db, gen);
        writeFastaChunks(outputBasePath, fastaSequence, fileCounter, buffer, currentSize);
    }

    if (!buffer.empty()) {
        std::ofstream outFile(outputBasePath + "-" + std::to_string(fileCounter) + ".fasta");
        outFile << buffer;
        outFile.close();
    }

    if (db != nullptr) {  // Ensure the database connection is closed properly
        sqlite3_close(db);
    }

    pcrTableFile.close();

    return 0;
}
