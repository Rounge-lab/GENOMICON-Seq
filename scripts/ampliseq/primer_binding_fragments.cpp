// GENOMICON-Seq candidate: strand-directed pairing, name-defined F/R roles.
// Compatible with v1.2.3 primer BED (0-based half-open) and fragment BED
// (1-based inclusive), output: fragment_name;genome_coordinates;lengths.
// Usage: primer_binding_fragments primers.bed fragments.bed.gz lengths.csv
//        output.csv MIN_LENGTH [CIRCULAR(Y|N)] [MAX_AMPLICON_LENGTH]
// Missing MAX_AMPLICON_LENGTH defaults to 1000 bp (paired PCR ONLY).
// MIN_LENGTH is applied to every generated product.
// F-only/R-only TaME-style fragment-boundary reactions do NOT use the maximum.
// The alignment strand of F selects its search direction; it never changes
// the F/R roles. R's mapping strand does not veto an experimentally validated pair.
// Compile: g++ -O2 -std=c++17 -Wall -Wextra -pedantic primer_binding_fragments_directional.cpp -lz -o primer_binding_fragments_directional

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>

using i64 = std::int64_t;

struct Primer {
    i64 start, end;              // 1-based inclusive, reference coordinates
    char role, strand;           // F/R from name; +/- from BED
};
struct Fragment {
    std::string ref, name;
    i64 start, end;              // 1-based inclusive, may wrap for circular
};
struct Bound {
    const Primer* primer;
    i64 start, end;              // unwrapped on physical fragment
};
struct PrimerPool {
    std::vector<Primer> hits;
    bool hasF = false, hasR = false;
};
using PrimerIndex = std::unordered_map<std::string, PrimerPool>;
using LengthIndex = std::unordered_map<std::string, i64>;

static std::string dequote(std::string s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    s = s.substr(a, b - a + 1);
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') ||
                          (s.front() == '\'' && s.back() == '\'')))
        s = s.substr(1, s.size() - 2);
    return s;
}

static LengthIndex readLengths(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open fasta lengths: " + path);
    LengthIndex lengths;
    std::string line;
    std::getline(in, line); // CSV header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        if (comma == std::string::npos) throw std::runtime_error("Bad fasta lengths row: " + line);
        const std::string ref = dequote(line.substr(0, comma));
        const i64 L = std::stoll(dequote(line.substr(comma + 1)));
        if (ref.empty() || L <= 0) throw std::runtime_error("Invalid fasta length: " + line);
        lengths[ref] = L;
    }
    if (lengths.empty()) throw std::runtime_error("No fasta lengths were loaded");
    return lengths;
}

static PrimerIndex readPrimers(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open primer BED: " + path);
    PrimerIndex byRef;
    std::string line;
    size_t n = 0;
    while (std::getline(in, line)) {
        ++n;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        std::string ref, name, score;
        i64 start, end;
        char strand;
        if (!(row >> ref >> start >> end >> name >> score >> strand))
            throw std::runtime_error("Primer BED must have 6 fields at line " + std::to_string(n));
        if ((strand != '+' && strand != '-') || start < 0 || end <= start)
            throw std::runtime_error("Invalid primer strand/coordinates at line " + std::to_string(n));
        if (name.size() < 3 || name[name.size()-2] != '-' ||
            (name.back() != 'F' && name.back() != 'R'))
            throw std::runtime_error("Primer name must end in -F or -R: " + name);
        auto& pool = byRef[ref];
        pool.hasF |= (name.back() == 'F');
        pool.hasR |= (name.back() == 'R');
        pool.hits.push_back({start + 1, end, name.back(), strand});
    }
    if (byRef.empty()) throw std::runtime_error("No primers in primer BED");
    for (auto& item : byRef)
        std::sort(item.second.hits.begin(), item.second.hits.end(),
                  [](const Primer& a, const Primer& b) { return a.start < b.start; });
    return byRef;
}

static Fragment parseFragment(const std::string& line) {
    std::istringstream row(line);
    Fragment f;
    if (!(row >> f.ref >> f.start >> f.end >> f.name))
        throw std::runtime_error("Invalid fragment BED row: " + line);
    return f;
}

static i64 wrap(i64 x, i64 L) {
    return ((x - 1) % L + L) % L + 1;
}

// Iterate the gzip file without loading millions of fragments into memory.
template <typename Consumer>
static void gzLines(const std::string& path, Consumer consume) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("Cannot open gzipped fragments: " + path);
    char buf[65536];
    std::string pending;
    try {
        while (true) {
            const int n = gzread(f, buf, sizeof(buf));
            if (n < 0) {
                int code = 0;
                const char* msg = gzerror(f, &code);
                throw std::runtime_error("Error reading gzipped fragments: " + std::string(msg));
            }
            if (n == 0) break;
            pending.append(buf, static_cast<size_t>(n));
            size_t from = 0;
            size_t pos;
            while ((pos = pending.find('\n', from)) != std::string::npos) {
                if (pos > from) consume(pending.substr(from, pos - from));
                from = pos + 1;
            }
            if (from) pending.erase(0, from);
        }
        if (!pending.empty()) consume(pending);
        gzclose(f);
    } catch (...) {
        gzclose(f);
        throw;
    }
}

static void process(
    const Fragment& f, const PrimerIndex& primers, const LengthIndex& lengths,
    bool circular, i64 minimum, i64 maximum, std::ofstream& out,
    std::uint64_t& fragmentCount, std::uint64_t& productCount
) {
    const auto lenIt = lengths.find(f.ref);
    if (lenIt == lengths.end()) throw std::runtime_error("Unknown reference: " + f.ref);
    const i64 L = lenIt->second;
    if (f.start < 1 || f.end < 1 || f.start > L || f.end > L)
        throw std::runtime_error("Fragment outside reference: " + f.name);
    if (!circular && f.start > f.end)
        throw std::runtime_error("Origin-spanning fragment on linear reference: " + f.name);
    const bool crossesOrigin = circular && f.start > f.end;
    const i64 arcStart = f.start;
    const i64 arcEnd = f.end + (crossesOrigin ? L : 0);
    const i64 physicalLength = arcEnd - arcStart + 1;
    const auto it = primers.find(f.ref);
    if (it == primers.end()) return;
    const auto& pool = it->second.hits;
    // Primer roles are precomputed per reference, avoiding an O(primer-count)
    // scan for every one of millions of chromosome fragments.
    const bool hasF = it->second.hasF;
    const bool hasR = it->second.hasR;

    std::vector<Bound> bound;
    auto gather = [&](i64 lo, i64 hi, i64 offset) {
        auto p = std::lower_bound(pool.begin(), pool.end(), lo,
            [](const Primer& item, i64 key) { return item.start < key; });
        for (; p != pool.end() && p->start <= hi; ++p) {
            if (p->end <= hi && p->start >= lo)
                bound.push_back({&*p, p->start + offset, p->end + offset});
        }
    };
    if (crossesOrigin) {
        gather(f.start, L, 0);
        gather(1, f.end, L);
    } else {
        gather(f.start, f.end, 0);
    }
    if (bound.empty()) return;

    std::vector<std::pair<i64, i64>> coords;
    std::vector<i64> sizes;
    auto append = [&](i64 left, i64 right, bool paired) {
        const i64 n = right - left + 1;
        if (left < arcStart || right > arcEnd || left > right ||
            n < minimum || n > physicalLength || (paired && n > maximum)) return;
        coords.emplace_back(wrap(left, L), wrap(right, L));
        sizes.push_back(n);
    };

    if (hasF && hasR) {
        for (const auto& forward : bound) {
            if (forward.primer->role != 'F') continue;
            for (const auto& reverse : bound) {
                if (reverse.primer->role != 'R') continue;
                if (forward.primer->strand == '+') {
                    // F extends toward increasing reference coordinates.
                    // The R binding site must be strictly to its right.
                    if (reverse.start > forward.end)
                        append(forward.start, reverse.end, true);
                } else {
                    // F extends toward decreasing reference coordinates.
                    // The R binding site must be strictly to its left.
                    if (reverse.end < forward.start)
                        append(reverse.start, forward.end, true);
                }
                // R mapping strand deliberately does not veto the pair:
                // primer roles and validation are supplied by the user.
            }
        }
    } else {
        // Preserve GENOMICON's TaME-style single-primer rules:
        // F-only: from F binding END to physical fragment's right boundary.
        // R-only: from physical fragment's left boundary to R binding START.
        // These special rules follow primer ROLE, not the mapped strand.
        // The paired-PCR maximum is not applied to these products.
        for (const auto& p : bound) {
            if (p.primer->role == 'F') append(p.end, arcEnd, false);
            else append(arcStart, p.start, false);
        }
    }

    if (coords.empty()) return;
    ++fragmentCount;
    productCount += coords.size();
    out << f.name << ';';
    for (size_t i = 0; i < coords.size(); ++i) {
        if (i) out << ',';
        out << coords[i].first << ':' << coords[i].second;
    }
    out << ';';
    for (size_t i = 0; i < sizes.size(); ++i) {
        if (i) out << ',';
        out << sizes[i];
    }
    out << '\n';
}

int main(int argc, char** argv) {
    if (argc < 6 || argc > 8) {
        std::cerr << "Usage: " << argv[0]
                  << " primers.bed fragments.bed.gz fasta_lengths.csv output.csv"
                  << " MIN_LENGTH [Y|N] [MAX_AMPLICON_LENGTH=1000]\n";
        return 1;
    }
    try {
        const i64 minimum = std::stoll(argv[5]);
        const std::string flag = argc >= 7 ? argv[6] : "N";
        const i64 maximum = argc == 8 ? std::stoll(argv[7]) : 1000;
        if (flag != "Y" && flag != "N")
            throw std::runtime_error("CIRCULAR must be Y or N");
        if (minimum < 1 || maximum < 1)
            throw std::runtime_error("Minimum and maximum lengths must be positive");
        const auto lengths = readLengths(argv[3]);
        const auto primers = readPrimers(argv[1]);
        std::ofstream out(argv[4]);
        if (!out) throw std::runtime_error("Cannot write PCR output: " + std::string(argv[4]));
        out << "fragment_name;genome_coordinates;lengths\n";

        std::optional<Fragment> last;
        std::uint64_t fragmentsWithProducts = 0, products = 0;
        auto flush = [&]() {
            if (!last) return;
            process(*last, primers, lengths, flag == "Y", minimum, maximum,
                    out, fragmentsWithProducts, products);
            last.reset();
        };
        gzLines(argv[2], [&](const std::string& line) {
            Fragment current = parseFragment(line);
            if (last && last->ref == current.ref && last->name == current.name) {
                if (flag == "Y" && last->end == lengths.at(current.ref) &&
                    current.start == 1) {
                    // GENOMICON writes [start,L] followed by [1,end]
                    // under the SAME fragment identifier.
                    last->end = current.end;
                    return;
                }
                throw std::runtime_error("Unexpected duplicate fragment ID: " + current.name);
            }
            flush();
            last = std::move(current);
        });
        flush();
        out.close();
        if (products == 0)
            throw std::runtime_error("No eligible amplicons found; check primer sites, F direction, length bounds and fragments");
        std::cerr << "Primer binding: " << fragmentsWithProducts
                  << " fragments yielded " << products << " candidate products"
                  << "; circular=" << flag << "; max paired=" << maximum << " bp\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR primer_binding_fragments_directional: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
