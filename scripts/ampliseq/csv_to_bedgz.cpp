// csv_to_bedgz.cpp
// Convert fragment_coordinates.csv + new_start.csv into a BED.GZ
// g++ -O3 -std=c++17 -pthread csv_to_bedgz.cpp -lboost_iostreams -o csv_to_bedgz

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <thread>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

inline std::string dequote(std::string s) {
    if (!s.empty() && (s.front()=='\"' || s.front()=='\'')) s.erase(0,1);
    if (!s.empty() && (s.back() =='\"' || s.back() =='\'')) s.pop_back();
    return s;
}
inline bool is_integer(const std::string& s) {
    static const std::regex re(R"(^\s*-?\d+\s*$)");
    return std::regex_match(s,re);
}

// ---------------------------------------------------------------------
// Worker: convert a CSV slice to a gz-compressed BED chunk
// ---------------------------------------------------------------------
void processChunk(int startLine, int endLine,
                  const std::string& inCSV, const std::string& outPart,
                  const std::map<std::string,int>& newStart,
                  const std::map<std::string,int>& fastaLen)
{
    std::ifstream csv(inCSV);
    if (!csv) { std::cerr<<"Cannot open "<<inCSV<<'\n'; return; }

    std::ofstream bed(outPart, std::ios::binary);
    boost::iostreams::filtering_ostream gzout;
    gzout.push(boost::iostreams::gzip_compressor());
    gzout.push(bed);

    std::string line; int cur=0;
    while (std::getline(csv,line)) {
        ++cur; if (cur<startLine || cur>endLine || cur==1) continue;   // skip header/range

        std::vector<std::string> fields;
        boost::split(fields,line,boost::is_any_of(";"));
        if (fields.size()<2) continue;

        std::string copy   = dequote(fields[0]);
        std::string ref    = copy.substr(0, copy.find('_'));

        auto ns=newStart.find(copy); auto fl=fastaLen.find(ref);
        if (ns==newStart.end() || fl==fastaLen.end()) continue;

        int shift=ns->second, gLen=fl->second;

        std::vector<std::string> coords;
        boost::split(coords,fields[1],boost::is_any_of(","));

        for (size_t i=0;i<coords.size();++i) {
            std::vector<std::string> p; boost::split(p,coords[i],boost::is_any_of(":"));
            if (p.size()!=2){ std::cerr<<"Warn: bad token '"<<coords[i]<<"'\n"; continue;}
            boost::trim(p[0]); boost::trim(p[1]);
            if (!is_integer(p[0])||!is_integer(p[1])){ std::cerr<<"Warn: bad ints '"<<coords[i]<<"'\n"; continue;}

            int start=std::stoi(p[0])+shift, end=std::stoi(p[1])+shift;

            if (start>gLen && end>gLen){ start-=gLen; end-=gLen; }
            else if (end>gLen){
                gzout<<ref<<'\t'<<start<<'\t'<<gLen<<'\t'<<copy<<"_f"<<i+1<<'\n';
                start=1; end-=gLen;
            }
            gzout<<ref<<'\t'<<start<<'\t'<<end<<'\t'<<copy<<"_f"<<i+1<<'\n';
        }
    }
    gzout.reset();
}

// ---------------------------------------------------------------------
// Merge part-files and clean up
// ---------------------------------------------------------------------
void mergeParts(const std::vector<std::string>& parts,const std::string& out){
    std::ofstream final(out, std::ios::binary);
    for (auto& p:parts){ std::ifstream in(p, std::ios::binary); final<<in.rdbuf(); std::remove(p.c_str()); }
}

// ---------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------
int main(int argc,char* argv[]){
    if (argc!=6){
        std::cerr<<"Usage: "<<argv[0]<<" <fragment_coordinates.csv> <new_start.csv> <fasta_lengths.csv> <output.bed.gz> <threads>\n";
        return 1;
    }
    std::string coordCSV=argv[1], newCSV=argv[2], lenCSV=argv[3], outBED=argv[4];
    int threads=std::stoi(argv[5]);

    // ---------- load new_start.csv (grab last column = numeric new_start)
    std::map<std::string,int> newStart;
    {
        std::ifstream f(newCSV); std::string l; std::getline(f,l);
        while (std::getline(f,l)){
            std::vector<std::string> c; boost::split(c,l,boost::is_any_of(","));
            if (c.empty()) continue;
            std::string key=dequote(c[0]), val=dequote(c.back());
            boost::trim(val); if (!is_integer(val)) continue;
            newStart[key]=std::stoi(val);
        }
    }

    // ---------- load fasta_lengths.csv
    std::map<std::string,int> fastaLen;
    {
        std::ifstream f(lenCSV); std::string l; std::getline(f,l);
        while (std::getline(f,l)){
            std::vector<std::string> c; boost::split(c,l,boost::is_any_of(","));
            if (c.size()<2) continue;
            std::string key=dequote(c[0]), val=dequote(c[1]);
            boost::trim(val); if (!is_integer(val)) continue;
            fastaLen[key]=std::stoi(val);
        }
    }

    // ---------- dispatch threads
    std::ifstream tmp(coordCSV);
    int total=std::count(std::istreambuf_iterator<char>(tmp),std::istreambuf_iterator<char>(),'\n');
    int per=total/threads, start=2;
    std::vector<std::thread> th; std::vector<std::string> parts;

    for(int i=0;i<threads;++i){
        int end=(i==threads-1)?total:start+per-1;
        std::string part="part_"+std::to_string(i)+".bed.gz"; parts.push_back(part);
        th.emplace_back(processChunk,start,end,coordCSV,part,std::cref(newStart),std::cref(fastaLen));
        start=end+1;
    }
    for(auto& t:th) t.join();
    mergeParts(parts,outBED);
    return 0;
}
