#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <thread>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

void processChunk(int startLine, int endLine, const std::string& inputFileName, const std::string& outputFileName, 
                  const std::map<std::string, int>& newStartPositions, const std::map<std::string, int>& fastaLengths) {
    std::ifstream csv_file(inputFileName);
    if (!csv_file.is_open()) {
        std::cerr << "Error opening input file: " << inputFileName << std::endl;
        return;
    }

    std::ofstream bed_file(outputFileName, std::ios_base::out | std::ios_base::binary);
    if (!bed_file.is_open()) {
        std::cerr << "Error opening output file: " << outputFileName << std::endl;
        return;
    }

    boost::iostreams::filtering_ostream out;
    out.push(boost::iostreams::gzip_compressor());
    out.push(bed_file);

    std::string line;
    int currentLine = 0;
    while (std::getline(csv_file, line)) {
        currentLine++;
        if (currentLine < startLine || currentLine > endLine || currentLine == 1) continue; 

        std::vector<std::string> fields;
        boost::split(fields, line, boost::is_any_of(";"));

        std::vector<std::string> name_parts;
        boost::split(name_parts, fields[0], boost::is_any_of("_"));
        std::string first_column = name_parts[0];

        std::vector<std::string> coordinates;
        boost::split(coordinates, fields[1], boost::is_any_of(","));

        int adjustment = newStartPositions.at(fields[0]);
        int fastaLength = fastaLengths.at(first_column);

        for (int i = 0; i < coordinates.size(); i++) {
            std::vector<std::string> pair;
            boost::split(pair, coordinates[i], boost::is_any_of(":"));

            int start = std::stoi(pair[0]) + adjustment;
            int end = std::stoi(pair[1]) + adjustment;

            if (start > fastaLength && end > fastaLength) {
                start = start - fastaLength + 1;
                end = end - fastaLength + 1;
            } else if (end > fastaLength) {
                out << first_column << "\t" << start << "\t" << fastaLength << "\t" << fields[0] << "_f" << i+1 << "\n";
                start = 1;
                end = end - fastaLength + 1;
            }

            out << first_column << "\t" << start << "\t" << end << "\t" << fields[0] << "_f" << i+1 << "\n";
        }
    }

    csv_file.close();
    out.pop();
    bed_file.close();
}

void mergeFiles(const std::vector<std::string>& tempFiles, const std::string& finalOutput) {
    std::ofstream finalFile(finalOutput, std::ios_base::out | std::ios_base::binary);
    if (!finalFile.is_open()) {
        std::cerr << "Error opening final output file: " << finalOutput << std::endl;
        return;
    }

    for (const auto& tempFile : tempFiles) {
        std::ifstream inFile(tempFile, std::ios_base::in | std::ios_base::binary);
        if (!inFile.is_open()) {
            std::cerr << "Error opening temporary file for merging: " << tempFile << std::endl;
            continue;
        }
        finalFile << inFile.rdbuf();
        inFile.close();
        std::remove(tempFile.c_str()); 
    }

    finalFile.close();
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <input.csv> <new_start.csv> <fasta_lengths.csv> <output.bed.gz> <num_threads>\n";
        return 1;
    }

    std::map<std::string, int> newStartPositions;
    std::ifstream newStartFile(argv[2]);
    if (!newStartFile.is_open()) {
        std::cerr << "Error opening new start positions file: " << argv[2] << std::endl;
        return 1;
    }
    std::string newStartLine;
    std::getline(newStartFile, newStartLine); // Skip header
    while (std::getline(newStartFile, newStartLine)) {
        std::vector<std::string> newStartFields;
        boost::split(newStartFields, newStartLine, boost::is_any_of(","));
        std::string genomeName = newStartFields[0].substr(1, newStartFields[0].length() - 2);
        int newStart = std::stoi(newStartFields[1]);
        newStartPositions[genomeName] = newStart;
    }
    newStartFile.close();

    std::map<std::string, int> fastaLengths;
    std::ifstream fastaFile(argv[3]);
    if (!fastaFile.is_open()) {
        std::cerr << "Error opening fasta lengths file: " << argv[3] << std::endl;
        return 1;
    }
    std::string fastaLine;
    std::getline(fastaFile, fastaLine); // Skip header
    while (std::getline(fastaFile, fastaLine)) {
        std::vector<std::string> fastaFields;
        boost::split(fastaFields, fastaLine, boost::is_any_of(","));
        std::string fastaName = fastaFields[0].substr(1, fastaFields[0].length() - 2);
        int fastaLength = std::stoi(fastaFields[1]);
        fastaLengths[fastaName] = fastaLength;
    }
    fastaFile.close();

    int numThreads = std::stoi(argv[5]);
    std::vector<std::thread> threads;
    std::vector<std::string> tempFileNames;

    std::ifstream tempFile(argv[1]);
    if (!tempFile.is_open()) {
        std::cerr << "Error opening input CSV file: " << argv[1] << std::endl;
        return 1;
    }
    int totalLines = std::count(std::istreambuf_iterator<char>(tempFile), std::istreambuf_iterator<char>(), '\n');
    tempFile.close();

    int linesPerThread = totalLines / numThreads;
    int startLine = 2; // Start from line 2 to skip header

    for (int i = 0; i < numThreads; ++i) {
        std::string tempFileName = "temp_output_" + std::to_string(i) + ".bed.gz";
        tempFileNames.push_back(tempFileName);
        int endLine = (i == numThreads - 1) ? totalLines : startLine + linesPerThread - 1;
        threads.emplace_back(processChunk, startLine, endLine, argv[1], tempFileName, newStartPositions, fastaLengths);
        startLine = endLine + 1;
    }

    for (auto& thread : threads) {
        thread.join();
    }

    mergeFiles(tempFileNames, argv[4]);

    return 0;
}
