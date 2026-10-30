// Heuristic priority only. Every undecided position still requires fallback.
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Seed { int query, offset; };
struct Hit { int score, length; std::string reference; };

int base(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

bool better(const Hit& a, const Hit& b) {
    if (a.score != b.score) return a.score > b.score;
    if (a.length != b.length) return a.length < b.length;
    return a.reference < b.reference;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::ifstream input(argv[1]);
    if (!input) return 3;
    std::vector<std::string> ids, targets;
    std::vector<Seed> postings[65536];
    std::string id, seq;
    while (input >> id >> seq) {
        if (seq.size() != 21) return 4;
        std::string desired;
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            int b = base(*it);
            if (b < 0) return 5;
            desired += "TGCA"[b];
        }
        int qi = ids.size();
        ids.push_back(id);
        targets.push_back(desired);
        for (int offset = 0; offset <= 13; ++offset) {
            unsigned key = 0;
            for (int j = 0; j < 8; ++j)
                key = (key << 2) | base(desired[offset + j]);
            postings[key].push_back({qi, offset});
        }
    }
    std::vector<std::vector<Hit>> best(ids.size());
    while (std::cin >> id >> seq) {
        std::unordered_map<int, int> scores;
        unsigned key = 0;
        int valid = 0;
        for (int pos = 0; pos < static_cast<int>(seq.size()); ++pos) {
            int b = base(seq[pos]);
            if (b < 0) { key = 0; valid = 0; continue; }
            key = ((key << 2) | b) & 65535;
            if (++valid < 8) continue;
            for (const auto& seed : postings[key]) {
                int start = pos - 7 - seed.offset;
                if (start < 0 || start + 21 > static_cast<int>(seq.size()))
                    continue;
                const auto& desired = targets[seed.query];
                int paired = 0, run = 0, longest = 0, gc = 0;
                for (int j = 0; j < 21; ++j) {
                    char t = seq[start + j], d = desired[j];
                    bool exact = t == d;
                    bool gu = (d == 'C' && t == 'T') ||
                              (d == 'A' && t == 'G');
                    paired += exact || gu;
                    run = exact ? run + 1 : 0;
                    longest = std::max(longest, run);
                    gc += exact && (d == 'G' || d == 'C');
                }
                if (paired < 16 && longest < 16) continue;
                int score = paired * 100 + longest * 2 + gc;
                auto found = scores.find(seed.query);
                if (found == scores.end() || found->second < score)
                    scores[seed.query] = score;
            }
        }
        for (const auto& entry : scores) {
            auto& hits = best[entry.first];
            Hit hit{entry.second, static_cast<int>(seq.size()), id};
            if (hits.size() >= 24 && !better(hit, hits.back())) continue;
            hits.push_back(hit);
            std::sort(hits.begin(), hits.end(), better);
            if (hits.size() > 24) hits.pop_back();
        }
    }
    for (unsigned qi = 0; qi < ids.size(); ++qi)
        for (const auto& hit : best[qi])
            std::cout << ids[qi] << '\t' << hit.reference << '\t'
                      << hit.score << '\n';
    return 0;
}
