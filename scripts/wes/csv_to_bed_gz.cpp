#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <thread>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

void processFile(const std::string& inputFileName, const std::string& outputDirectory, 
                 const std::map<std::string, int>& fastaLengths) {
    std::ifstream csv_file(inputFileName);
    std::string outputFileName = outputDirectory + "/" + 
                                 std::filesystem::path(inputFileName).stem().string() + ".bed.gz";
    std::ofstream bed_file(outputFileName, std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_ostream out;
    out.push(boost::iostreams::gzip_compressor());
    out.push(bed_file);

    std::string line;
    while (std::getline(csv_file, line)) {
        std::vector<std::string> fields;
        boost::split(fields, line, boost::is_any_of(";"));
        if (fields.size() != 2) continue;

        std::string genomeNameFull = fields[0];
        boost::erase_all(genomeNameFull, "\"");
        std::string genomeName = genomeNameFull.substr(0, genomeNameFull.find("_"));

        if (fastaLengths.find(genomeName) == fastaLengths.end()) {
            std::cerr << "Error: Genome name not found in fastaLengths: " << genomeName << std::endl;
            continue;
        }

        int fastaLength = fastaLengths.at(genomeName);

        std::vector<std::string> coordinates;
        boost::split(coordinates, fields[1], boost::is_any_of(","));

        for (size_t i = 0; i < coordinates.size(); i++) {
            std::vector<std::string> pair;
            boost::split(pair, coordinates[i], boost::is_any_of(":"));
            if (pair.size() != 2) continue;

            int start = std::stoi(pair[0]) + 1;
            int end = std::stoi(pair[1]) + 1;

            out << genomeName << "\t" << start << "\t" << end << "\t" << genomeNameFull << "_f" << i + 1 << "\n";
        }
    }

    csv_file.close();
    out.pop();
    bed_file.close();

    try {
        if (std::filesystem::exists(inputFileName)) {
            std::filesystem::remove(inputFileName);
        } else {
            std::cerr << "Warning: File not found for deletion: " << inputFileName << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error deleting file: " << e.what() << '\n';
    }
}

std::vector<std::string> listBatchFiles(const std::string& directoryPath) {
    std::vector<std::string> batchFiles;
    for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
        std::string filename = entry.path().filename().string();
        if (filename.find("_batch") != std::string::npos && filename.find(".csv") != std::string::npos) {
            batchFiles.push_back(entry.path().string());
        }
    }
    return batchFiles;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <batch_files_directory> <fasta_lengths.csv> <output_directory> <num_threads>\n";
        return 1;
    }

    std::string batchFilesDirectory = argv[1];
    std::string fastaFile = argv[2];
    std::string outputDirectory = argv[3];
    int numThreads = std::stoi(argv[4]);

    std::vector<std::string> batchFiles = listBatchFiles(batchFilesDirectory);

    std::map<std::string, int> fastaLengths;
    std::ifstream fastaFileStream(fastaFile);
    if (!fastaFileStream.is_open()) {
        std::cerr << "Error opening fasta lengths file: " << fastaFile << std::endl;
        return 1;
    }

    std::string fastaLine;
    std::getline(fastaFileStream, fastaLine);
    while (std::getline(fastaFileStream, fastaLine)) {
        std::vector<std::string> fastaFields;
        boost::split(fastaFields, fastaLine, boost::is_any_of(","));
        std::string fastaName = fastaFields[0];
        int fastaLength = std::stoi(fastaFields[1]);
        fastaLengths[fastaName] = fastaLength;
    }
    fastaFileStream.close();

    std::vector<std::thread> threads;
    for (size_t i = 0; i < batchFiles.size(); ++i) {
        threads.emplace_back(processFile, batchFiles[i], outputDirectory, fastaLengths);
        if (threads.size() == numThreads || i == batchFiles.size() - 1) {
            for (auto& t : threads) {
                t.join();
            }
            threads.clear();
        }
    }

    return 0;
}
