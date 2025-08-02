#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

struct Primer {
    std::string ref;
    int start;
    int end;
    char type; // 'F' or 'R'
};

struct Fragment {
    std::string ref;  
    int start;
    int end;
    std::string name; 
};

std::map<std::string, int> readGenomeLengths(const std::string& filename) {
    std::map<std::string, int> genomeLengths;
    std::ifstream file(filename);
    std::string line, ref, lengthStr;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        std::istringstream iss(line);
        getline(iss, ref, ',');
        getline(iss, lengthStr);
        ref = ref.substr(1, ref.length() - 2); 
        genomeLengths[ref] = std::stoi(lengthStr);
    }
    return genomeLengths;
}

std::vector<Primer> readPrimers(const std::string& filename) {
    std::vector<Primer> primers;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Primer p;
        std::string primerName;
        iss >> p.ref >> p.start >> p.end >> primerName; 
        p.type = primerName.substr(primerName.find_last_of('-') + 1).back(); 
        primers.push_back(p);
    }
    return primers;
}

void processFragments(const std::string& filename, const std::vector<Primer>& primers, const std::map<std::string, int>& genomeLengths, const std::string& outputFilename, int validDistance) {
    std::ofstream outputFile(outputFilename);
    outputFile << "fragment_name;genome_coordinates;lengths" << std::endl; 

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        Fragment f;
        iss >> f.ref >> f.start >> f.end;
        iss.ignore(std::numeric_limits<std::streamsize>::max(), '\t'); 
        iss >> f.name; 

        bool isCircularFragment = f.start > f.end;
        int genomeLength = genomeLengths.at(f.ref);

        std::vector<std::pair<int, int>> validPairs; 
        std::vector<int> lengths; 

        if (isCircularFragment) {
            if (f.start == 0) f.start = genomeLength;
            if (f.end == 0) f.end = 1; 
        }

        Primer defaultF = {f.ref, f.start, f.start, 'F'};
        Primer defaultR = {f.ref, f.end, f.end, 'R'};
        bool hasF = false, hasR = false;

        for (const auto& primer : primers) {
            if (primer.ref == f.ref) {
                if (primer.type == 'F') hasF = true;
                if (primer.type == 'R') hasR = true;
            }
        }

        for (const auto& primerF : (hasF ? primers : std::vector<Primer>{defaultF})) {
            if (primerF.type != 'F' || primerF.ref != f.ref) continue;
            for (const auto& primerR : (hasR ? primers : std::vector<Primer>{defaultR})) {
                if (primerR.type != 'R' || primerR.ref != f.ref) continue;

                bool isFWithin = (primerF.end >= f.start && primerF.end <= f.end);
                bool isRWithin = (primerR.start >= f.start && primerR.start <= f.end);

                if (isFWithin && isRWithin) {
                    int fEnd = primerF.end;
                    int rStart = primerR.start;
                    if (isCircularFragment && fEnd > rStart) {
                        rStart += genomeLength; 
                    }

                    int distance = rStart - fEnd;
                    if (distance < 0) distance += genomeLength; 
                    if (distance >= validDistance) {
                        validPairs.emplace_back(fEnd % genomeLength ? fEnd % genomeLength : genomeLength, rStart % genomeLength ? rStart % genomeLength : genomeLength);
                        lengths.push_back(distance);
                    }
                }
            }
        }

        if (!validPairs.empty()) {
            outputFile << f.name << ";";
            for (size_t i = 0; i < validPairs.size(); ++i) {
                if (i > 0) outputFile << ",";
                outputFile << validPairs[i].first << ":" << validPairs[i].second;
            }
            outputFile << ";";
            for (size_t i = 0; i < lengths.size(); ++i) {
                if (i > 0) outputFile << ",";
                outputFile << lengths[i];
            }
            outputFile << std::endl;
        }
    }
}




int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <primers.bed> <fragments.bed.gz> <genome_length.csv> <output.csv> <valid_distance>\n";
        return 1;
    }

    std::string primerFile = argv[1];
    std::string fragmentFile = argv[2];
    std::string genomeLengthFile = argv[3];
    std::string outputFilename = argv[4];
    int validDistance = std::stoi(argv[5]);

    auto genomeLengths = readGenomeLengths(genomeLengthFile);
    auto primers = readPrimers(primerFile);
    processFragments(fragmentFile, primers, genomeLengths, outputFilename, validDistance);

    return 0;
}
