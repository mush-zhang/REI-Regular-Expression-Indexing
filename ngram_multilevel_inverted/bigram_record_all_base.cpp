#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <regex>
#include <chrono>
#include <numeric>
#include <random>
#include <tuple>
#include <queue>
#include <cmath>
#include <limits>
#include <ctime>
#include <unordered_map>
#include <unordered_set>

#include <re2/re2.h>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

const static int NUM_INDEXED = 2;
// Define base_generator as a Mersenne Twister. This is needed only to make the
// code a bit less verbose.
typedef boost::mt19937 base_generator;

// ------------------------------------------------------------------------------------------------------
//  String Splitting
// ------------------------------------------------------------------------------------------------------

bool bracket_capturing(const std::string &line, std::size_t pos) {
    auto temp_pos = pos;
    int num_escaped = 0;
    // count number of escape characters before the (
    while (temp_pos > 0) {
        if (line.at(--temp_pos) != '\\') 
            break;
        num_escaped++;
    }
    return num_escaped % 2 == 0;
}

std::tuple<std::string, std::string, std::string> split_regex(const std::string &line) {
    std::size_t pos = 0;
    std::string prefix = line;
    std::string regex = "";
    std::string suffix = "";
    while ((pos = line.find("(", pos)) != std::string::npos) {
        if (bracket_capturing(line, pos)) {
            std::size_t reg_start_pos = pos;
            std::size_t reg_end_pos = pos+1;
            while ((pos = line.find(")", pos+1)) != std::string::npos) {
                if (bracket_capturing(line, pos)) {
                    reg_end_pos = pos;
                }
                pos++;
            }
            prefix = line.substr(0, reg_start_pos);
            regex = line.substr(reg_start_pos, reg_end_pos+1-reg_start_pos);
            suffix = line.substr(reg_end_pos+1);
            
            break;
        }
        pos++;
    }
    prefix.erase(std::remove(prefix.begin(), prefix.end(), '\\'), prefix.end()); // remove escape (trick method. should only remove single \ and turning \\ to \)
    suffix.erase(std::remove(suffix.begin(), suffix.end(), '\\'), suffix.end()); // remove escape (trick method. should only remove single \ and turning \\ to \)
    return std::make_tuple(prefix, regex, suffix);
}

std::tuple<std::vector<std::string>, std::vector<std::string>, bool> split_regex_multi(const std::string &line) {
    std::size_t pos = 0;
    std::size_t prev_pos = 0;
    int pos2 = 0;
    std::vector<std::string> const_strings;
    std::vector<std::string> regexes;
    std::string prefix = line;
    std::string suffix = "";
    bool prefix_first = true;
    while ((pos = line.find("(", pos)) != std::string::npos) {
        if (bracket_capturing(line, pos)) {
            prefix = line.substr(prev_pos, pos-prev_pos);
            if (prev_pos == 0 && pos == 0) prefix_first = false;
            pos2 = pos;
            while ((pos2 = line.find(")", pos2)) != std::string::npos) {
                if (bracket_capturing(line, pos2)) {
                    suffix = line.substr(pos, pos2 + 1 -pos);
                    break;
                }
                pos2++;
            }
            
            if (pos2 == std::string::npos && suffix.empty()) {
                suffix = line.substr(pos);
            }
            pos2++; 
            regexes.push_back(suffix);
            if (!prefix.empty()) {
                prefix.erase(std::remove(prefix.begin(), prefix.end(), '\\'), prefix.end()); // remove escape (trick method. should only remove single \ and turning \\ to \)
                const_strings.push_back(prefix);
            }
            
            prev_pos = pos2;
            pos = pos2;
        } else {
            pos++;
        }
    }

    prefix = line.substr(pos2);
    if (!prefix.empty()) {
        prefix.erase(std::remove(prefix.begin(), prefix.end(), '\\'), prefix.end()); // remove escape (trick method. should only remove single \ and turning \\ to \)
        const_strings.push_back(prefix);
    }
    return std::tuple<std::vector<std::string>, std::vector<std::string>, bool>(const_strings, regexes, prefix_first);
}

// ------------------------------------------------------------------------------------------------------
//  BLARE & RE2 Baseline
// ------------------------------------------------------------------------------------------------------

std::pair<double, int> Re2DirectMatch(const std::vector<std::string> & lines, const std::string & reg_string) {
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::string sm;
    RE2 reg(reg_string);
    for (const auto & line : lines) {
        count += RE2::PartialMatch(line, reg, &sm);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    return std::make_pair(elapsed_seconds.count(), count);
}

bool MultiMatchSingle (const std::string & line, std::vector<std::shared_ptr<RE2>> & c_regs, std::shared_ptr<RE2> & reg0, const std::vector<std::string> prefixes, const std::vector<std::string> & regs, bool prefix_first, std::vector<size_t> & prev_prefix_pos) {
    std::string sm;
    bool match = false;
    if (regs.empty()) {
        return line.find(prefixes[0]) != std::string::npos;
    } else if (prefixes.empty()) {
        return RE2::PartialMatch(line, *reg0, &sm);
    } else {
        size_t pos = 0;
        size_t curr_prefix_pos = 0;
        size_t prefix_idx = 0;
        size_t reg_idx = 0;
        MATCH_LOOP_SINGLE:
            for (; prefix_idx < prefixes.size(); ) {
                // find pos of prefix before reg
                if ((curr_prefix_pos = line.find(prefixes[prefix_idx], pos)) == std::string::npos) {
                    if (prefix_idx == 0 || prev_prefix_pos[prefix_idx] == 0 || reg_idx == 0)
                        goto CONTINUE_OUTER_SINGLE;
                    else {
                        prefix_idx--;
                        reg_idx--;
                        pos = prev_prefix_pos[prefix_idx]+1;
                        continue;
                    }
                }  
                pos = curr_prefix_pos + prefixes[prefix_idx].size();
                prev_prefix_pos[prefix_idx] = curr_prefix_pos;
                if (prefix_idx == prev_prefix_pos.size()-1 || prev_prefix_pos[prefix_idx+1] >= pos) {
                    break;
                }
                prefix_idx++;
            }
            for (; reg_idx < prev_prefix_pos.size()-1; reg_idx++) {
                size_t prev_prefix_end_pos = prev_prefix_pos[reg_idx] + prefixes[reg_idx].size();
                auto curr = line.substr(prev_prefix_end_pos, prev_prefix_pos[reg_idx+1] - prev_prefix_end_pos);
                if (!RE2::FullMatch(curr, *(c_regs[reg_idx]), &sm)){
                    pos = prev_prefix_pos[reg_idx+1]+1;
                    prefix_idx = reg_idx+1;
                    goto MATCH_LOOP_SINGLE;
                }
            }
            if (prefixes.size() == regs.size() && prefix_first) {
                auto curr = line.substr(pos);
                re2::StringPiece input(curr);
                if (!RE2::Consume(&input, *(c_regs.back()), &sm)) {
                    prefix_idx = prev_prefix_pos.size() -1;
                    pos = prev_prefix_pos[prefix_idx]+1;
                    goto MATCH_LOOP_SINGLE;
                } 
            }
            if (!prefix_first) {
                auto curr = line.substr(0,prev_prefix_pos[0]);
                if (!RE2::PartialMatch(curr, *reg0, &sm)){
                    prefix_idx = 0;
                    pos = prev_prefix_pos[prefix_idx]+1;
                    goto MATCH_LOOP_SINGLE;
                }
            }
            match = true;                
        CONTINUE_OUTER_SINGLE:;
        std::fill(prev_prefix_pos.begin(), prev_prefix_pos.end(), 0);
    }
    return match;
}

bool SplitMatchSingle (const std::string & line, RE2 & reg, const std::tuple<std::string, std::string, std::string> & r) {
    bool match = false;
    std::string sm;
    auto prefix = std::get<0>(r);
    auto suffix = std::get<2>(r);
    if (std::get<1>(r).empty()) {
        match = line.find(prefix) != std::string::npos;
    } else if (prefix.empty()) {
        if (suffix.empty()) {
            match = RE2::PartialMatch(line, reg, &sm);
        } else {
            std::size_t pos = 0;
            while ((pos = line.find(suffix, pos)) != std::string::npos) {
                std::string curr_in = line.substr(0, pos); 
                re2::StringPiece input(curr_in);
                while (RE2::Consume(&input, reg, &sm)) {
                    if (input.ToString().empty()) {
                        match = true;
                        goto NEXT_LINE2;
                    } 
                }
                pos++;
            }
            NEXT_LINE2:;
        }
    } else if (suffix.empty()) {
        std::size_t pos = 0;
        while ((pos = line.find(prefix, pos)) != std::string::npos) {
            // for accuracy, should use line = line.substr(pos+1) next time
            std::string curr_in = line.substr(pos + prefix.length()); 
            re2::StringPiece input(curr_in);
            if (RE2::Consume(&input, reg, &sm)) {
                match = true;
                break;
            }
            pos++;
        }                       
    } else {
        std::size_t pos = 0;
        while ((pos = line.find(prefix, pos)) != std::string::npos) {
            std::size_t reg_start_pos = pos + prefix.length();
            std::size_t reg_end_pos = reg_start_pos;
            while ((reg_end_pos = line.find(suffix, reg_end_pos)) != std::string::npos) {
                if (RE2::FullMatch(line.substr(reg_start_pos, reg_end_pos - reg_start_pos ), reg, &sm)) {
                    match = true;
                    goto NEXT_LINE;
                }
                reg_end_pos++;
            }
            pos++;
        } 
        NEXT_LINE:;  
    }

    return match;
}

bool FullMatchSingle (const std::string & line, RE2 & reg) {
    std::string sm;
    return RE2::PartialMatch(line, reg, &sm);
}

// argmax returns the index of maximum element in vector v.
template<class T>
size_t argmax(const std::vector<T>& v){
  return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
}

// argmin returns the index of minimum element in vector v.
template<class T>
size_t argmin(const std::vector<T>& v){
  return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

std::tuple<double, int, unsigned int> Blare (const std::vector<std::string> & lines, const std::string & reg_string) {
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    unsigned int idx = 0;

    // pre_compile regs
    RE2 reg_full(reg_string);

    auto r = split_regex(reg_string);
    RE2 reg_suffix(std::get<1>(r));


    // Assumes prefix at first
    auto r_multi = split_regex_multi(reg_string);
    std::vector<std::string> prefixes = std::get<0>(r_multi);
    std::vector<std::string> regs_temp = std::get<1>(r_multi);
    auto regs = std::get<1>(r_multi);
    bool prefix_first = std::get<2>(r_multi);
    std::shared_ptr<RE2> reg0;
    if (!prefix_first) {
        reg0 = std::shared_ptr<RE2>(new RE2(regs_temp[0]+"$"));
        regs_temp.erase(regs_temp.begin());
    }
    std::vector<std::shared_ptr<RE2>> c_regs;
    for (const auto & reg : regs_temp) {
        std::shared_ptr<RE2> c_reg(new RE2(reg));
        c_regs.push_back(c_reg);
    }
    std::vector<size_t> prev_prefix_pos(prefixes.size(), 0);

    size_t skip_size = 100;
    size_t iteration_num = 10;
    size_t sample_size = lines.size() / 100000;
    if (sample_size < 200) {
        sample_size = 200;
        skip_size = 10;
    }
    else if (sample_size > 10000) 
        sample_size = 10000;

    base_generator gen(static_cast<std::uint32_t>(std::time(0)));
    boost::random::uniform_int_distribution<> dist{0, 2};

    // auto pred_results = std::vector<int>(iteration_num);
    std::vector<int> pred_results{0, 0, 0};
    auto arm_num = 3;

    for (; idx < skip_size; idx++) {
        // switch(b(rng)) {
        switch(dist(gen)) {
            case 0: count += SplitMatchSingle(lines[idx], reg_suffix, r); break;
            case 1: count += MultiMatchSingle(lines[idx], c_regs, reg0, prefixes, regs, prefix_first, prev_prefix_pos); break;
            case 2: count += FullMatchSingle(lines[idx], reg_full); break;
        }
    }

    for (size_t j = 0; j < iteration_num; j++) {
        // Number of trials per bandit
        auto trials = std::vector<unsigned int>(arm_num);
        // Average time per bandit
        auto ave_elapsed = std::vector<double>(arm_num);
        // Number of wins per bandit
        auto wins = std::vector<unsigned int>(arm_num);
        // Beta distributions of the priors for each bandit
        std::vector<boost::random::beta_distribution<>> prior_dists;
        // Initialize the prior distributions with alpha=1 beta=1
        for (size_t i = 0; i < arm_num; i++) {
            prior_dists.push_back(boost::random::beta_distribution<>(1, 1));
        }
        for (size_t k = 0; k < sample_size; k++, idx++) {
            std::vector<double> priors;
            // Sample a random value from each prior distribution.
            for (auto& dist : prior_dists) {
                priors.push_back(dist(gen));
            }
            // Select the bandit that has the highest sampled value from the prior
            size_t chosen_bandit = argmax(priors);
            
            // Pull the lever of the chosen bandit
            auto single_start = std::chrono::high_resolution_clock::now();
            switch(chosen_bandit) {
                case 0: count += SplitMatchSingle(lines[idx], reg_suffix, r); break;
                case 1: count += MultiMatchSingle(lines[idx], c_regs, reg0, prefixes, regs, prefix_first, prev_prefix_pos); break;
                case 2: count += FullMatchSingle(lines[idx], reg_full); break;
            }
            
            auto single_end = std::chrono::high_resolution_clock::now();
            auto single_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(single_end - single_start).count(); 
            trials[chosen_bandit] += 1;

            ave_elapsed[chosen_bandit] = (ave_elapsed[chosen_bandit] * (trials[chosen_bandit]-1) + single_elapsed) / ( 1.0 * trials[chosen_bandit]);
            wins[chosen_bandit] += argmin(ave_elapsed) == chosen_bandit;

            if (wins[chosen_bandit] >= sample_size / arm_num) 
                break;

            // Update the prior distribution of the chosen bandit
            auto alpha = 1 + wins[chosen_bandit];
            auto beta = 1 + trials[chosen_bandit] - wins[chosen_bandit];
            prior_dists[chosen_bandit] = boost::random::beta_distribution<>(alpha, beta);
        }
        auto winning_strategy = argmax(std::vector<double>{(wins[0]*1.0)/trials[0], (wins[1]*1.0)/trials[1], (wins[2]*1.0)/trials[2]});
        // pred_results[j] = winning_strategy;
        pred_results[winning_strategy]++;
        if (pred_results[winning_strategy] >= iteration_num / arm_num) 
            break;
    }
    auto chosen_bandit = argmax(pred_results);
    std::string sm;
    switch(chosen_bandit) {
        case 0: {
            auto prefix = std::get<0>(r);
            auto suffix = std::get<2>(r);
            if (std::get<1>(r).empty()) {
                for (; idx < lines.size(); idx++) {
                    count += lines[idx].find(prefix) != std::string::npos;
                }
            } else if (prefix.empty()) {
                if (suffix.empty()) {
                    for (; idx < lines.size(); idx++) {
                        count += RE2::PartialMatch(lines[idx], reg_suffix, &sm);
                    }
                } else {
                    for (; idx < lines.size(); idx++) {
                        std::size_t pos = 0;
                        while ((pos = lines[idx].find(suffix, pos)) != std::string::npos) {
                            std::string curr_in = lines[idx].substr(0, pos); 
                            re2::StringPiece input(curr_in);
                            while (RE2::Consume(&input, reg_suffix, &sm)) {
                                if (input.ToString().empty()) {
                                    count++;
                                    goto NEXT_LINE2;
                                } 
                            }
                            pos++;
                        }
                        NEXT_LINE2:;
                    }
                }
            } else if (suffix.empty()) {
                for (; idx < lines.size(); idx++) {
                    std::size_t pos = 0;
                    while ((pos = lines[idx].find(prefix, pos)) != std::string::npos) {
                        std::string curr_in = lines[idx].substr(pos + prefix.length()); 
                        re2::StringPiece input(curr_in);
                        if (RE2::Consume(&input, reg_suffix, &sm)) {
                            count++;
                            break;
                        }
                        pos++;
                    }
                }  
            } else {
                for (; idx < lines.size(); idx++) {
                    std::size_t pos = 0;
                    while ((pos = lines[idx].find(prefix, pos)) != std::string::npos) {
                        std::size_t reg_start_pos = pos + prefix.length();
                        std::size_t reg_end_pos = reg_start_pos;
                        while ((reg_end_pos = lines[idx].find(suffix, reg_end_pos)) != std::string::npos) {
                            if (RE2::FullMatch(lines[idx].substr(reg_start_pos, reg_end_pos - reg_start_pos ), reg_suffix, &sm)) {
                                count++;
                                goto NEXT_LINE;
                            }
                            reg_end_pos++;
                        }
                        pos++;
                    } 
                    NEXT_LINE:;
                }
                
            }
            break;
        }
        case 1: {
            if (regs.empty()) {
                for (; idx < lines.size(); idx++) {
                    count += lines[idx].find(prefixes[0]) != std::string::npos;
                }
            } else if (prefixes.empty()) {
                RE2 reg(regs[0]);
                for (; idx < lines.size(); idx++) {
                    count += RE2::PartialMatch(lines[idx], reg, &sm);
                }
            } else {
                std::vector<size_t> prev_prefix_pos(prefixes.size(), 0);
                for (; idx < lines.size(); idx++) {
                    auto line = lines[idx];
                    size_t pos = 0;
                    size_t curr_prefix_pos = 0;
                    size_t prefix_idx = 0;
                    size_t reg_idx = 0;
                    MATCH_LOOP_BLARE:
                        for (; prefix_idx < prefixes.size(); ) {
                            // find pos of prefix before reg
                            if ((curr_prefix_pos = line.find(prefixes[prefix_idx], pos)) == std::string::npos) {
                                if (prefix_idx == 0 || prev_prefix_pos[prefix_idx] == 0 || reg_idx == 0)
                                    goto CONTINUE_OUTER_BLARE;
                                else {
                                    prefix_idx--;
                                    reg_idx--;
                                    pos = prev_prefix_pos[prefix_idx]+1;
                                    continue;
                                }
                            }  
                            pos = curr_prefix_pos + prefixes[prefix_idx].size();
                            prev_prefix_pos[prefix_idx] = curr_prefix_pos;
                            if (prefix_idx == prev_prefix_pos.size()-1 || prev_prefix_pos[prefix_idx+1] >= pos) {
                                break;
                            }
                            prefix_idx++;
                        }
                        for (; reg_idx < prev_prefix_pos.size()-1; reg_idx++) {
                            size_t prev_prefix_end_pos = prev_prefix_pos[reg_idx] + prefixes[reg_idx].size();
                            auto curr = line.substr(prev_prefix_end_pos, prev_prefix_pos[reg_idx+1] - prev_prefix_end_pos);
                            if (!RE2::FullMatch(curr, *(c_regs[reg_idx]), &sm)){
                                pos = prev_prefix_pos[reg_idx+1]+1;
                                prefix_idx = reg_idx+1;
                                goto MATCH_LOOP_BLARE;
                            }
                        }
                        if (prefixes.size() == regs.size() && prefix_first) {
                            auto curr = line.substr(pos);
                            re2::StringPiece input(curr);
                            if (!RE2::Consume(&input, *(c_regs.back()), &sm)) {
                                prefix_idx = prev_prefix_pos.size() -1;
                                pos = prev_prefix_pos[prefix_idx]+1;
                                goto MATCH_LOOP_BLARE;
                            } 
                        }
                        if (!prefix_first) {
                            auto curr = line.substr(0,prev_prefix_pos[0]);
                            if (!RE2::PartialMatch(curr, *reg0, &sm)){
                                prefix_idx = 0;
                                pos = prev_prefix_pos[prefix_idx]+1;
                                goto MATCH_LOOP_BLARE;
                            }
                        }
                        count++;            
                    CONTINUE_OUTER_BLARE:;
                    std::fill(prev_prefix_pos.begin(), prev_prefix_pos.end(), 0);
                }
            }
            break;
        }
        case 2: {
            for (; idx < lines.size(); idx++) {
                count += RE2::PartialMatch(lines[idx], reg_full, &sm);
            } 
            break;
        }
    }

    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    return std::make_tuple(elapsed_seconds.count(), count, chosen_bandit);
}

// ------------------------------------------------------------------------------------------------------
//  Tuple Hash Functions
//  Hashing strategy of tuples: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3876.pdf
// ------------------------------------------------------------------------------------------------------

template <typename T>
void hash_combine (std::size_t& seed, const T& val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// auxiliary generic functions to create a hash value using a seed
template <typename T, typename... Types>
void hash_combine (std::size_t& seed, const T& val, const Types&... args) {
    hash_combine(seed,val);
    hash_combine(seed,args...);
}

// optional auxiliary generic functions to support hash_val() without arguments
void hash_combine (std::size_t& seed) {}

//  generic function to create a hash value out of a heterogeneous list of arguments
template <typename... Types>
std::size_t hash_val (const Types&... args) {
    std::size_t seed = 0;
    hash_combine(seed, args...);
    return seed;
}

struct hash_pair {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const
    {
        return hash_val(p.first, p.second);
    }
};

// ------------------------------------------------------------------------------------------------------
//  Input Data
// ------------------------------------------------------------------------------------------------------

std::vector<std::string> read_traffic() {
    std::string line;
    std::vector<std::string> lines;
    std::string data_file = "../US_Accidents_Dec21_updated.csv";
    std::ifstream data_in(data_file);
    if (!data_in.is_open()) {
        std::cerr << "Could not open the file - '" << data_file << "'" << std::endl;
        return lines;
    }
    while (getline(data_in, line)){
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::stringstream streamData(line);
        std::string s; 
        int i = 0;
        while (getline(streamData, s, ',')) {
            if (i++ == 9){
                lines.push_back(s);
                break;
            }
        }
    }
    data_in.close();
    return lines;
}

std::vector<std::string> read_db_x() {
    std::string line;

    std::vector<std::string> lines;
    std::string path = "../extracted";
    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        std::string data_file = entry.path();
        std::ifstream data_in(data_file);
        if (!data_in.is_open()) {
            std::cerr << "Could not open the file - '" << data_file << "'" << std::endl;
            return lines;
        }
        //  101876733 
        while (getline(data_in, line)){
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            if (line.size() > 2)
                lines.push_back(line);
        }
        data_in.close();
    }
    return lines;
}

// ------------------------------------------------------------------------------------------------------
//  N-gram Extraction & Index Building
// ------------------------------------------------------------------------------------------------------

std::unordered_set<std::pair<char, char>, hash_pair> make_unique_bigrams(const std::string& s) {
    std::unordered_set<std::pair<char, char>, hash_pair> bigrams;
    if (s.size() >= 2) {
        for(auto it = s.cbegin(); it != std::prev(s.cend()); ++it) {
            bigrams.insert(std::make_pair(*it , *std::next(it)));
        }
    }
    return bigrams;
}

void insert_unique_bigrams_to_set(std::unordered_set<std::pair<char, char>, hash_pair> & bigrams, const std::string & s) {
    if (s.size() >= 2) {
        for(auto it = s.cbegin(); it != std::prev(s.cend()); ++it) {
            bigrams.insert(std::make_pair(*it , *std::next(it)));
        }
    }
}

std::unordered_map<std::pair<char, char>, size_t, hash_pair> top_query_ngrams2(const std::vector<std::string> & regexes, int k) {

    // count number of occurence for all ngrams
    std::unordered_map<std::pair<char, char>, size_t, hash_pair> ngram_counts;
    std::vector<std::pair<char, char>> keys;
    
    for (const auto & reg_string : regexes) {
        auto curr_literals = std::get<0>(split_regex_multi(reg_string));
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams;
        for (const auto & lit : curr_literals) {
            insert_unique_bigrams_to_set(curr_grams, lit);
        }
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_counts.find(line_gram);
            if (it != ngram_counts.end()) {
                ngram_counts[line_gram]++;
            } else {
                ngram_counts.insert({line_gram, 1});
                keys.push_back(line_gram);
            }
        }
    }

    // sort by number of occurence
    std::sort(
        keys.begin(), 
        keys.end(),
        [&ngram_counts](std::pair<char, char> const& a, std::pair<char, char> const& b) {
            if (a == b) {
                return false;
            }
            if (ngram_counts[a] > ngram_counts[b]) {
                return true;
            }
            else if (ngram_counts[a] < ngram_counts[b]) {
                return false;
            }
            return a < b;
    });

    std::unordered_map<std::pair<char, char>, size_t, hash_pair> result;
    for (auto i = 0; i < NUM_INDEXED; i++) {
        result.insert({keys[i], i});
    }

    return result;
}

std::unordered_map<std::pair<char, char>, size_t, hash_pair> top_query_ngrams(const std::vector<std::string> & regex_lits, int k) {

    // count number of occurence for all ngrams
    std::unordered_map<std::pair<char, char>, size_t, hash_pair> ngram_counts;
    std::vector<std::pair<char, char>> keys;
    for (const auto & line : regex_lits) {
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_counts.find(line_gram);
            if (it != ngram_counts.end()) {
                ngram_counts[line_gram]++;
            } else {
                ngram_counts.insert({line_gram, 1});
                keys.push_back(line_gram);
            }
        }
    }

    // sort by number of occurence
    std::sort(
        keys.begin(), 
        keys.end(),
        [&ngram_counts](std::pair<char, char> const& a, std::pair<char, char> const& b) {
            if (a == b) {
                return false;
            }
            if (ngram_counts[a] > ngram_counts[b]) {
                return true;
            }
            else if (ngram_counts[a] < ngram_counts[b]) {
                return false;
            }
            return a < b;
    });

    std::unordered_map<std::pair<char, char>, size_t, hash_pair> result;
    for (auto i = 0; i < NUM_INDEXED; i++) {
        result.insert({keys[i], i});
    }

    return result;
}

// return a vector of m-bits array, each corresponds to existence of some bigrams in the log
std::vector<std::bitset<NUM_INDEXED>> build_index_single(const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::vector<std::string> & lines) {

    std::vector<std::bitset<NUM_INDEXED>> result;
    for (const auto & line : lines) {
        // init a bitset of all 0's
        std::bitset<NUM_INDEXED> curr_bitarr;
        auto curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            auto it = gram_to_idx.find(line_gram);
            if (it != gram_to_idx.end()) {
                // bigram exits, set bit at the index
                curr_bitarr.set(it->second);
            }
        }
        result.push_back(curr_bitarr);
    }
    return result;
}

std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> build_inverted_index_single(const std::vector<std::string> & regex_lits, const std::vector<std::string> & lines) {
    // create an entry for each ngram in the regex
    std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted;
    for (const auto & line : regex_lits) {
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            if (ngram_inverted.find(line_gram) == ngram_inverted.end()) {
                // ngram_inverted.emplace(std::pair<char, char>, std::vector<unsigned int>());
                ngram_inverted[line_gram] = std::vector<unsigned int>(); // Initialize key with null vector
            }
        }
    }

    for (size_t idx = 0; idx < lines.size(); idx++) {
        auto line = lines[idx];
        auto curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_inverted.find(line_gram);
            if (it != ngram_inverted.end()) {
                // bigram exits, add line idx to the inverted index
                ngram_inverted[line_gram].push_back(idx);
            }
        }
    }
    return ngram_inverted;
}

std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> build_topk_inverted_index_single(const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::vector<std::string> & lines) {
    // create an entry for each ngram
    std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted;
    for (const auto & keyval: gram_to_idx) {
        ngram_inverted[keyval.first] = std::vector<unsigned int>(); // Initialize key with empty vector
    }

    for (size_t idx = 0; idx < lines.size(); idx++) {
        auto line = lines[idx];
        auto curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_inverted.find(line_gram);
            if (it != ngram_inverted.end()) {
                // bigram exits, add line idx to the inverted index
                ngram_inverted[line_gram].push_back(idx);
            }
        }
    }
    return ngram_inverted;
}

// return a vector of m-bits array, each corresponds to existence of some bigrams in the log
std::vector<std::bitset<NUM_INDEXED>> build_index(const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::vector<std::string> & lines, int granularity) {
    if (granularity == 1) {
        return build_index_single(gram_to_idx, lines);
    }
    std::vector<std::bitset<NUM_INDEXED>> result;  
    size_t line_idx = 0;
    for (size_t group_idx = 0; group_idx < std::ceil(lines.size()/double(granularity)); group_idx++) {
        std::bitset<NUM_INDEXED> curr_bitarr;
        // get all bigrams for group of log lines
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams;
        for (size_t subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); subline_idx++, line_idx++) {
            insert_unique_bigrams_to_set(curr_grams, lines[line_idx]);
        }
        for (const auto & line_gram: curr_grams) {
            auto it = gram_to_idx.find(line_gram);
            if (it != gram_to_idx.end()) {
                // bigram exits, set bit at the index
                curr_bitarr.set(it->second);
            }
        }
        result.push_back(curr_bitarr);
    }
    return result;
}

std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> build_inverted_index(const std::vector<std::string> & regex_lits, const std::vector<std::string> & lines, int granularity) {
    if (granularity == 1) {
        return build_inverted_index_single(regex_lits, lines);
    }
    // create an entry for each ngram in the regex
    std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted;
    for (const auto & line : regex_lits) {
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams = make_unique_bigrams(line);
        for (const auto & line_gram: curr_grams) {
            if (ngram_inverted.find(line_gram) == ngram_inverted.end()) {
                ngram_inverted[line_gram] = std::vector<unsigned int>(); // Initialize key with null vector
            }
        }
    }

    size_t line_idx = 0;
    for (size_t group_idx = 0; group_idx < std::ceil(lines.size()/double(granularity)); group_idx++) {
        // get all bigrams for group of log lines
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams;
        for (size_t subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); subline_idx++, line_idx++) {
            insert_unique_bigrams_to_set(curr_grams, lines[line_idx]);
        }
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_inverted.find(line_gram);
            if (it != ngram_inverted.end()) {
                // bigram exits, add line idx to the inverted index
                ngram_inverted[line_gram].push_back(group_idx);
            }
        }
    }
    return ngram_inverted;
}

std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> build_topk_inverted_index(const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::vector<std::string> & lines, int granularity) {
    if (granularity == 1) {
        return build_topk_inverted_index_single(gram_to_idx, lines);
    }
    // create an entry for each ngram
    std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted;
    for (const auto & keyval: gram_to_idx) {
        ngram_inverted[keyval.first] = std::vector<unsigned int>(); // Initialize key with empty vector
    }
    
    size_t line_idx = 0;
    for (size_t group_idx = 0; group_idx < std::ceil(lines.size()/double(granularity)); group_idx++) {
        // get all bigrams for group of log lines
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams;
        for (size_t subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); subline_idx++, line_idx++) {
            insert_unique_bigrams_to_set(curr_grams, lines[line_idx]);
        }
        for (const auto & line_gram: curr_grams) {
            auto it = ngram_inverted.find(line_gram);
            if (it != ngram_inverted.end()) {
                // bigram exits, add line idx to the inverted index
                ngram_inverted[line_gram].push_back(group_idx);
            }
        }
    }
    return ngram_inverted;
}

// ------------------------------------------------------------------------------------------------------
//  Match with N-gram Index Filtering
// ------------------------------------------------------------------------------------------------------

std::bitset<NUM_INDEXED> GetBitMask(const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::string & reg_string) {
    // get the bit mask for the regex
    std::vector<std::string> reg_lits = std::get<0>(split_regex_multi(reg_string));
    std::unordered_set<std::pair<char, char>, hash_pair> reg_bigrams{};
    for (const auto & lit : reg_lits) {
        reg_bigrams.merge(make_unique_bigrams(lit));
    }
    std::bitset<NUM_INDEXED> bitmask;
    bitmask.set();
    // start with all 1's
    for (const auto & reg_gram: reg_bigrams) {
        auto it = gram_to_idx.find(reg_gram);
        if (it != gram_to_idx.end()) {
            // bigram exits, unset bit at the index
            bitmask.reset(it->second);
        }
    }
    // return a bitmask where the positions we want to tests are 0's, and others are 1's.
    return bitmask;
}

std::pair<double, int> Re2FilterMatchSingle(const std::vector<std::string> & lines, const std::vector<std::bitset<NUM_INDEXED>> & idxs, const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::string & reg_string) {
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::string sm;
    RE2 reg(reg_string);
    auto bitmask = GetBitMask(gram_to_idx, reg_string);
    for (auto i = 0; i < lines.size(); i++) {
        // check if the indexed ngrams in regular expressions are present in the current log line
        if ((idxs[i] | bitmask).all()) {
            count += RE2::PartialMatch(lines[i], reg, &sm);
        } 
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    return std::make_pair(elapsed_seconds.count(), count);
}

std::pair<double, int> Re2FilterMatch(const std::vector<std::string> & lines, const std::vector<std::bitset<NUM_INDEXED>> & idxs, const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::string & reg_string, int granularity) {
    if (granularity == 1) {
        return Re2FilterMatchSingle(lines, idxs, gram_to_idx, reg_string);
    }
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::string sm;
    RE2 reg(reg_string);
    auto bitmask = GetBitMask(gram_to_idx, reg_string);

    for (size_t group_idx = 0; group_idx < std::ceil(lines.size()/double(granularity)); group_idx++) {
        if ((idxs[group_idx] | bitmask).all()) {
            auto line_idx = group_idx*granularity;
            for (auto subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); line_idx++, subline_idx++) {
                count += RE2::PartialMatch(lines[line_idx], reg, &sm);
            }
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    // std::cout << "FilterTime: " << elapsed_seconds.count() << " Count " << count << std::endl;
    return std::make_pair(elapsed_seconds.count(), count);
}

std::pair<double, int> Re2InvertedFilterMatch2Single(const std::vector<std::string> & lines, std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted, const std::string & reg_string) {
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::string sm;
    RE2 reg(reg_string);
    auto curr_literals = std::get<0>(split_regex_multi(reg_string));
    std::unordered_set<std::pair<char, char>, hash_pair> unique_bigrams{};
    for (const auto & curr_lit : curr_literals) {
        unique_bigrams.merge(make_unique_bigrams(curr_lit));
    } 

    // find shortest
    std::vector<std::pair<char, char>> checking_ngram{};
    std::pair<char, char> shortest_ngram;
    int shortest_size = -1;
    bool first = true;

    for (const auto & ngram : unique_bigrams) {
        if (ngram_inverted.find(ngram) != ngram_inverted.end()) {
            int curr_size = ngram_inverted[ngram].size();
            if (first) {
                shortest_size = curr_size;
                shortest_ngram = ngram;
                first = false;
                if (curr_size = 0) {
                    break;
                }
            } else if (curr_size < shortest_size) {
                shortest_size = curr_size;
                if (curr_size = 0) {
                    break;
                }
                checking_ngram.push_back(shortest_ngram);
                shortest_ngram = ngram;
            } else {
                checking_ngram.push_back(ngram);
            }
        }
    }

    if (shortest_size == -1) {
        for (auto i = 0; i < lines.size(); i++) {
            count += RE2::PartialMatch(lines[i], reg, &sm);
        }
    } else if (shortest_size > 0) {
        std::vector<size_t> pos(checking_ngram.size(), 0);
        std::vector<unsigned int> candidates = ngram_inverted[shortest_ngram];
        std::vector<unsigned int> final_candidates = std::vector<unsigned int>();
        bool filtered = false;
        for (const auto & cand_line_idx : candidates) {
            for (size_t i = 0; i < checking_ngram.size(); i++) {
                auto curr_ngram = checking_ngram[i];
                auto it = std::lower_bound(ngram_inverted[curr_ngram].begin() + pos[i], ngram_inverted[curr_ngram].end(), cand_line_idx);
                pos[i] = it - ngram_inverted[curr_ngram].begin();
                if (it == ngram_inverted[curr_ngram].end()) {
                    // all unchecked idxs would be greater than the current one
                    // and none of them can be found
                    goto ALL_CANDIDATES_FOUND;
                }
                else if ((it != candidates.end()) && (*it != cand_line_idx)) {
                    // this idx does not satisfy, continue with next index in the candidates
                    goto NEXT_CANDIDATE;
                }
            }
            // all ngram has this index
            final_candidates.push_back(cand_line_idx);
            NEXT_CANDIDATE:;
        }
        ALL_CANDIDATES_FOUND:;

        // match on the indexes
        for (auto & i : candidates) {
            count += RE2::PartialMatch(lines[i], reg, &sm);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    // std::cout << "Time2: " << elapsed_seconds.count() << " Count2 " << count << std::endl;
    return std::make_pair(elapsed_seconds.count(), count);
}

std::pair<double, int> Re2InvertedFilterMatch2(const std::vector<std::string> & lines, std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> ngram_inverted, const std::string & reg_string, int granularity) {
    if (granularity == 1) {
        return Re2InvertedFilterMatch2Single(lines, ngram_inverted, reg_string);
    }    
    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    std::string sm;
    RE2 reg(reg_string);
    auto curr_literals = std::get<0>(split_regex_multi(reg_string));
    std::unordered_set<std::pair<char, char>, hash_pair> unique_bigrams{};
    for (const auto & curr_lit : curr_literals) {
        unique_bigrams.merge(make_unique_bigrams(curr_lit));
    } 

    // find shortest
    std::vector<std::pair<char, char>> checking_ngram{};
    std::pair<char, char> shortest_ngram;
    int shortest_size = -1;
    bool first = true;

    for (const auto & ngram : unique_bigrams) {
        if (ngram_inverted.find(ngram) != ngram_inverted.end()) {
            int curr_size = ngram_inverted[ngram].size();
            if (first) {
                shortest_size = curr_size;
                shortest_ngram = ngram;
                first = false;
                if (curr_size = 0) {
                    break;
                }
            } else if (curr_size < shortest_size) {
                shortest_size = curr_size;
                if (curr_size = 0) {
                    break;
                }
                checking_ngram.push_back(shortest_ngram);
                shortest_ngram = ngram;
            } else {
                checking_ngram.push_back(ngram);
            }
        }
    }

    if (shortest_size == -1) {
        for (auto i = 0; i < lines.size(); i++) {
            count += RE2::PartialMatch(lines[i], reg, &sm);
        }
    } else if (shortest_size > 0) {
        std::vector<size_t> pos(checking_ngram.size(), 0);
        std::vector<unsigned int> candidates = ngram_inverted[shortest_ngram];
        std::vector<unsigned int> final_candidates = std::vector<unsigned int>();
        bool filtered = false;
        for (const auto & cand_line_idx : candidates) {
            for (size_t i = 0; i < checking_ngram.size(); i++) {
                auto curr_ngram = checking_ngram[i];
                auto it = std::lower_bound(ngram_inverted[curr_ngram].begin() + pos[i], ngram_inverted[curr_ngram].end(), cand_line_idx);
                pos[i] = it - ngram_inverted[curr_ngram].begin();
                if (it == ngram_inverted[curr_ngram].end()) {
                    // all unchecked idxs would be greater than the current one
                    // and none of them can be found
                    goto ALL_CANDIDATES_FOUND;
                }
                else if ((it != candidates.end()) && (*it != cand_line_idx)) {
                    // this idx does not satisfy, continue with next index in the candidates
                    goto NEXT_CANDIDATE;
                }
            }
            // all ngram has this index
            final_candidates.push_back(cand_line_idx);
            NEXT_CANDIDATE:;
        }
        ALL_CANDIDATES_FOUND:;

        // match on the indexes
        for (auto & i : candidates) {
            auto line_idx = i*granularity;
            for (auto subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); line_idx++, subline_idx++) {
                count += RE2::PartialMatch(lines[line_idx], reg, &sm);
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    // std::cout << "Time2: " << elapsed_seconds.count() << " Count2 " << count << std::endl;
    return std::make_pair(elapsed_seconds.count(), count);
}

// ------------------------------------------------------------------------------------------------------
//  Auxiliary Information Collection
// ------------------------------------------------------------------------------------------------------

int StatsAfterFiltering(const std::vector<std::string> & lines, const std::vector<std::bitset<NUM_INDEXED>> & idxs, const std::unordered_map<std::pair<char, char>, size_t, hash_pair> & gram_to_idx, const std::string & reg_string, int granularity) {
    auto bitmask = GetBitMask(gram_to_idx, reg_string);
    int count = 0;
    for (size_t group_idx = 0; group_idx < std::ceil(lines.size()/double(granularity)); group_idx++) {
        if ((idxs[group_idx] | bitmask).all()) {
            auto line_idx = group_idx*granularity;
            for (auto subline_idx = 0; subline_idx < granularity && line_idx < lines.size(); line_idx++, subline_idx++) {
                count++;
            }
        }

    }
    return count;
}

template <typename T>
unsigned long long calculate_vector_size(const std::vector<T> & vec) {
    unsigned long long base_size = sizeof(std::vector<T>);
    unsigned long long content_size = sizeof(T) * vec.size();
    return base_size + content_size;
}

unsigned long long calculate_map_size(const std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> & inv_idx, bool max) {
    unsigned long long entrySize = sizeof(std::pair<char, char>) + sizeof(void*);
    unsigned long long bucketSize = sizeof(void*);
    unsigned long long adminSize = 3 * sizeof(void*) + sizeof(size_t);
    unsigned long long totalSize = 0;
    if (max) {
        totalSize += adminSize + inv_idx.size() * entrySize + inv_idx.max_bucket_count() * bucketSize;
    } else {
        totalSize += adminSize + inv_idx.size() * entrySize + inv_idx.bucket_count() * bucketSize;
    } 

    unsigned long long contentSize = 0;
    for (const auto& kv : inv_idx) {
        contentSize += calculate_vector_size(kv.second);
    }
    totalSize += contentSize;

    return totalSize;
}

// ------------------------------------------------------------------------------------------------------
//  Main
// ------------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Incorrect number of arguments. " << std::endl;
        std::cerr << "\tRun with: ./bigram_record_all_[NUM_INDEXED].o [expr_idx] [granularity]" << std::endl;
        return EXIT_FAILURE;
    }
    
    const auto dataset_name = "db_x";

    std::string line;
    // read all regexes
    std::vector<std::string> regexes;
    std::vector<std::string> regex_lits;

    std::string reg_file = "../sampled_db_x_regex.csv";
    std::ifstream reg_in(reg_file);
    if (!reg_in.is_open()) {
        std::cerr << "Could not open the file - '" << reg_file << "'" << std::endl;
        return EXIT_FAILURE;
    }
    while (getline(reg_in, line)){
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        regexes.push_back(line);
    }
    reg_in.close();
    
    auto lines = read_db_x();
    
    std::cout << "Number of Logs: " << lines.size() << std::endl;
    
    int granularity = atoi(argv[2]);
    
    // row level index building time
    for (const auto & reg_string : regexes) {
        auto curr_literals = std::get<0>(split_regex_multi(reg_string));
        regex_lits.insert(regex_lits.end(), curr_literals.begin(), curr_literals.end());
    }
    std::unordered_map<std::pair<char, char>, size_t, hash_pair> gram_to_idx = top_query_ngrams2(regexes, NUM_INDEXED);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto top_idxs = build_index(gram_to_idx, lines, granularity);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto bit_index_building_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto bit_index_size = calculate_vector_size(top_idxs);
    std::cout << "Bitvector Index building: " << bit_index_building_time.count() << " seconds" << std::endl;
    std::cout << "Bitvector Index Size: " <<  bit_index_size << std::endl;

    start = std::chrono::high_resolution_clock::now();
    // std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> inverted_index = build_inverted_index(regex_lits, lines);
    std::unordered_map<std::pair<char, char>, std::vector<unsigned int>, hash_pair> inverted_index = build_topk_inverted_index(gram_to_idx, lines, granularity);
    end = std::chrono::high_resolution_clock::now();
    auto inverted_index_building_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    
    auto inverted_index_size_fit = calculate_map_size(inverted_index, false);
    auto inverted_index_size_max = calculate_map_size(inverted_index, true);
    std::cout << "Inverted Index building: " << inverted_index_building_time.count() << " seconds" << std::endl;
    std::cout << "Inverted Index Size: " << inverted_index_size_fit << " " << inverted_index_size_max << std::endl;
    
    // return EXIT_FAILURE;

    std::ofstream ind_out{"ngram_record_all_script/index_building.csv", std::ios::app};
    ind_out << "bigram," << NUM_INDEXED << "," << granularity << "," << bit_index_building_time.count() << "," << bit_index_size << "," << inverted_index_building_time.count() << "," << inverted_index_size_fit << "," << inverted_index_size_max << std::endl;
    ind_out.close();

    int num_repeat = 10;

    const std::filesystem::path dir_path = "ngram_record_all_script/new_bigram_" + std::to_string(NUM_INDEXED) + "_" + argv[2];
    if (!std::filesystem::exists(dir_path)) {
        std::filesystem::create_directory(dir_path);
    }

    std::vector<std::ofstream> expr_out_files;
    for (size_t i = 0; i < num_repeat; i++) {
        std::ofstream expr_r_file{dir_path/("bigram_" + std::string(std::string(argv[1])) + "_" + std::to_string(i) + ".csv"), std::ofstream::out};
        expr_r_file << "regex\tfilter_time\tinverted_time\tdirect_time" << std::endl;
        expr_out_files.push_back(std::move(expr_r_file));
    }

    std::ofstream r_file{dir_path/("summary" +  std::string(std::string(argv[1])) +".csv"), std::ofstream::out};
    r_file << "regex\tfilter_time\tmid_filter_time\tinverted_time\tmid_inverted_time\tdirect_time\tmiddirect_time\tfilter\tnum_after_filter\tmatch_num_filter\tmatch_num_inverted\tmatch_num_direct" << std::endl;

    for (const std::string & r : regexes) {
        auto bitmask = GetBitMask(gram_to_idx, r);
        std::vector<double> elapsed_time_direct;
        std::vector<double> elapsed_time_filter;
        std::vector<double> elapsed_time_inverted;
        int match_count_inverted;
        int match_count_direct;
        int match_count_filter;
        for (int i = 0; i < num_repeat; i++) {
            std::pair<double, int> result_filter = Re2FilterMatch(lines, top_idxs, gram_to_idx, r, granularity);
            std::tuple<double, int, unsigned int> result_blare = Blare(lines, r);
            std::pair<double, int> result_direct = std::make_pair(std::get<0>(result_blare), std::get<1>(result_blare));
            auto result_inverted = Re2InvertedFilterMatch2(lines, inverted_index, r, granularity);

            elapsed_time_direct.push_back(result_direct.first);
            elapsed_time_filter.push_back(result_filter.first);
            elapsed_time_inverted.push_back(result_inverted.first);
            expr_out_files[i] << r << "\t" << result_filter.first << "\t" << result_inverted.first << "\t" << result_direct.first << std::endl;

            match_count_inverted = result_inverted.second;
            match_count_direct = result_direct.second;
            match_count_filter = result_filter.second;
        }

        std::sort(elapsed_time_direct.begin(), elapsed_time_direct.end());
        auto ave_direct = std::accumulate(elapsed_time_direct.begin(), elapsed_time_direct.end(), 0.0) / elapsed_time_direct.size();
        auto mid_ave_direct = std::accumulate(elapsed_time_direct.begin()+1, elapsed_time_direct.end()-1, 0.0) / (num_repeat-2);

        std::sort(elapsed_time_filter.begin(), elapsed_time_filter.end());
        auto ave_filter = std::accumulate(elapsed_time_filter.begin(), elapsed_time_filter.end(), 0.0) / elapsed_time_filter.size();        
        auto mid_ave_filter = std::accumulate(elapsed_time_filter.begin()+1, elapsed_time_filter.end()-1, 0.0) / (num_repeat-2);
        
        std::sort(elapsed_time_inverted.begin(), elapsed_time_inverted.end());
        auto ave_invert = std::accumulate(elapsed_time_inverted.begin(), elapsed_time_inverted.end(), 0.0) / elapsed_time_inverted.size();        
        auto mid_ave_invert = std::accumulate(elapsed_time_inverted.begin()+1, elapsed_time_inverted.end()-1, 0.0) / (num_repeat-2);

        r_file << r << "\t";
        r_file << ave_filter << "\t" << mid_ave_filter << "\t";
        r_file << ave_invert << "\t" << mid_ave_invert << "\t";
        r_file << ave_direct << "\t" << mid_ave_direct << "\t";
        r_file << bitmask.to_string() << "\t" << StatsAfterFiltering(lines, top_idxs, gram_to_idx, r, granularity) << "\t";
        r_file << match_count_filter << "\t" << match_count_inverted << "\t" << match_count_direct << std::endl;
    }
    r_file.close();

    for (auto & expr_out : expr_out_files) {
        expr_out.close();
    }

}
// g++ -O3 -std=c++17   -Ofast -march=native -mfma -mavx -fomit-frame-pointer -ffp-contract=fast -flto -DARMA_NO_DEBUG -pthread bigram_record_all_base.cpp  -L/usr/local/lib/ -lre2 -lstdc++fs -o bigram_record_all_base.o