#include <cstdlib>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

// ------------------------------------------------------------------------------------------------------
//  String Splitting & Manipulating
// ------------------------------------------------------------------------------------------------------

std::string get_printable_char_helper(char c) {
    if (!std::isprint(c)) {
        return " ";
    }
    return std::string(1, c);
}

std::string get_printable_bigram(const std::pair<char, char> & b) {
    return get_printable_char_helper(b.first) + get_printable_char_helper(b.second);
}

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
        // while (getline(data_in, line)){
        while (lines.size() < 1000000 && getline(data_in, line)){ // subset
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

std::unordered_map<std::pair<char, char>, size_t, hash_pair> top_query_ngrams2(const std::vector<std::string> & regexes) {

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
    for (auto i = 0; i < keys.size(); i++) {
        result.insert({keys[i], i});
    }

    return result;
}

// return a map where key is all bigrams in regex, and value is number of regexes with the bigram
std::pair<std::unordered_map<std::pair<char, char>, double, hash_pair>, std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair>> bigram_regex_count(const std::vector<std::string> & regexes) {
    std::unordered_map<std::pair<char, char>, double, hash_pair> ngram_counts;
    std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair> ngram_exist_regexes;
    // for (const auto & reg_string : regexes) {
    for (size_t idx = 0; idx < regexes.size(); idx++) {
        auto reg_string = regexes[idx];
        auto curr_literals = std::get<0>(split_regex_multi(reg_string));
        std::unordered_set<std::pair<char, char>, hash_pair> curr_grams;
        for (const auto & lit : curr_literals) {
            insert_unique_bigrams_to_set(curr_grams, lit);
        }
        for (const auto & regex_gram: curr_grams) {

            auto it = ngram_counts.find(regex_gram);
            if (it == ngram_counts.end()) {
                ngram_counts.insert({regex_gram, 0});
                std::vector<size_t> vect;
                ngram_exist_regexes.insert({regex_gram, std::move(vect)});
            }
            ngram_counts[regex_gram]++;
            ngram_exist_regexes[regex_gram].push_back(idx);
        }
    }
    return std::make_pair(ngram_counts, ngram_exist_regexes);
}

// return a map where key is all bigrams in regex, and value is number of log lines without this q garm
std::pair<std::unordered_map<std::pair<char, char>, double, hash_pair>, std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair>> bigram_line_count(
    const std::vector<std::string> & lines, const std::unordered_map<std::pair<char, char>, double, hash_pair> & regex_grams) 
{
    std::unordered_map<std::pair<char, char>, double, hash_pair> ngram_not_exist_counts;
    std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair> ngram_exist_lines;
    for (size_t idx = 0; idx < lines.size(); idx++) {
        auto line = lines[idx];
        auto curr_grams = make_unique_bigrams(line);

        // check if bigrams in regexes are in the current line
        for (const auto & kv : regex_grams) {
            auto regex_gram = kv.first;
            if (ngram_not_exist_counts.find(regex_gram) == ngram_not_exist_counts.end()) {
                ngram_not_exist_counts.insert({regex_gram, 0});
                std::vector<size_t> vect;
                ngram_exist_lines.insert({regex_gram, std::move(vect)});
            }
            // if it is not in the line
            if (curr_grams.find(regex_gram) == curr_grams.end()) {
                ngram_not_exist_counts[regex_gram]++;
            } else {
                ngram_exist_lines[regex_gram].push_back(idx);
            }
        }
    }
    return std::make_pair(ngram_not_exist_counts, ngram_exist_lines);
}

// ------------------------------------------------------------------------------------------------------
//  Individual Benefit Calculation
// ------------------------------------------------------------------------------------------------------

std::pair<std::unordered_map<std::pair<char, char>, double, hash_pair>, std::vector<std::pair<char, char>>> individual_benefit(
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & line_not_exist_grams, 
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & regex_grams) {

    std::unordered_map<std::pair<char, char>, double, hash_pair> indi_benefit;
    std::vector<std::pair<char, char>> benefit_ranking;

    for (const auto & keyval : regex_grams) {
        auto bigram_key = keyval.first;
        auto num_r_w_key = keyval.second;

        auto num_l_wo_key = line_not_exist_grams.at(bigram_key);
        indi_benefit.insert({bigram_key, num_r_w_key*num_l_wo_key});
        benefit_ranking.push_back(bigram_key);
    }

    // sort by number of occurence
    std::sort(
        benefit_ranking.begin(), 
        benefit_ranking.end(),
        [&indi_benefit](std::pair<char, char> const& a, std::pair<char, char> const& b) {
            if (a == b) {
                return false;
            }
            if (indi_benefit[a] > indi_benefit[b]) {
                return true;
            }
            else if (indi_benefit[a] < indi_benefit[b]) {
                return false;
            }
            return a < b;
    });
    return std::make_pair(indi_benefit, benefit_ranking);
}

std::vector<std::pair<char, char>> individual_benefit_percentage(
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & line_not_exist_grams, size_t line_count, 
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & regex_grams, size_t regex_count) 
{    
    std::ofstream ind_out{"individual_bigram_benefit_db_x_subset.csv", std::ios::out};
    ind_out << "bigram\tnum_g\tbenefit_single\tabs_benefit\tbenefit_perc" << std::endl;

    std::unordered_map<std::pair<char, char>, double, hash_pair> indi_benefit;
    std::vector<std::pair<char, char>> benefit_ranking;
    for (const auto & keyval : regex_grams) {
        auto bigram_key = keyval.first;
        auto num_r_w_key = keyval.second;
        auto perc_r_w_key = num_r_w_key/regex_count;
        auto num_l_wo_key = line_not_exist_grams.at(bigram_key);
        auto perc_l_wo_key = num_l_wo_key/line_count;
        indi_benefit.insert({bigram_key, perc_r_w_key*perc_l_wo_key});
        ind_out << get_printable_bigram(bigram_key) << "\t" << num_r_w_key << "\t" << num_l_wo_key << "\t" << num_r_w_key*num_l_wo_key << "\t" << perc_r_w_key*perc_l_wo_key << std::endl;
        benefit_ranking.push_back(bigram_key);        
    }

    ind_out.close();
    // sort by number of occurence
    std::sort(
        benefit_ranking.begin(), 
        benefit_ranking.end(),
        [&indi_benefit](std::pair<char, char> const& a, std::pair<char, char> const& b) {
            if (a == b) {
                return false;
            }
            if (indi_benefit[a] > indi_benefit[b]) {
                return true;
            }
            else if (indi_benefit[a] < indi_benefit[b]) {
                return false;
            }
            return a < b;
    });
    return benefit_ranking;
}

// ------------------------------------------------------------------------------------------------------
//  incremental Benefit Calculation
// ------------------------------------------------------------------------------------------------------

std::pair<double, double> num_queries_helper(const std::vector<size_t> & vect1, const std::vector<size_t> & vect2, size_t regex_count) {
    double num_g1_g2 = 0;
    double num_g2_not_g1 = 0;
    size_t i = 0;
    size_t j = 0;
    for (size_t reg_idx = 0; reg_idx < regex_count; reg_idx++) {
        bool g1 = false;
        while (reg_idx > vect1[i]) {
            i++;
        }
        if (reg_idx == vect1[i]) {
            // current index in g1
            g1 = true;
        } 

        bool g2 = false;
        while (reg_idx > vect2[j]) {
            j++;
        }
        if (reg_idx == vect2[j]) {
            g2 = true;
        }
        
        if (g1 && g2) {
            num_g1_g2++;
            i++;
            j++;
        } else if ((!g1) && g2){
            num_g2_not_g1++;
        }
    }
    return std::make_pair(num_g1_g2, num_g2_not_g1);
}

void incremental_bigram_line(
    const std::vector<std::string> & lines, 
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & line_not_exist_grams, 
    const std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair> & ngram_exist_lines, size_t line_count,
    const std::unordered_map<std::pair<char, char>, double, hash_pair> & regex_grams, 
    const std::unordered_map<std::pair<char, char>, std::vector<size_t>, hash_pair> & ngram_exist_regexes, size_t regex_count,
    const std::vector<std::pair<char, char>> & benefit_ranking) 
{
    // key is condition bigram1
    // value is a map where 
    //    key is bigram2 used after condition filters, value is count of filtered out by bigram2
    std::unordered_map<std::pair<char, char>, std::unordered_map<std::pair<char, char>, double, hash_pair>, hash_pair> result;

    std::ofstream ind_out{"incremental_bigram_benefit_db_x_subset.csv", std::ios::out};
    ind_out << "bigram\tcondition_bigram\tincremental_benefit_single\tnum_g1_g2\tbenefit_g2_single\tnum_g2_not_g1\tnum_incremental_benefit\tg1_remaining_base\tincremental_benefit_perc" << std::endl;

    long double total_time = 0;
    size_t gram2_count = 0;
    auto start = std::chrono::high_resolution_clock::now();

    // for each bigram2 in regex
    for (const auto & bigram_key2 : benefit_ranking) {
        if (gram2_count++ > 200)
            break;
        auto g2_filtered = line_not_exist_grams.at(bigram_key2);

        if (g2_filtered == 0) {
            // filter nothing, don't bother
            continue;
        }

        size_t gram1_count = 0;

        // for each bigram1 in regex
        for (const auto & bigram_key1 : benefit_ranking) {
            if (gram1_count++ > 200) 
                break;
            if (bigram_key2.first == bigram_key1.first && bigram_key2.second == bigram_key1.second) {
                continue;
            }
            double incre_count = 0;

            auto g1_filtered = line_not_exist_grams.at(bigram_key1);
            auto [num_g1_g2, num_g2_not_g1] = num_queries_helper(ngram_exist_regexes.at(bigram_key1), ngram_exist_regexes.at(bigram_key2), regex_count);

            if (g1_filtered == 0) {
                // g1 filtered nothing
                incre_count = g2_filtered;
            } else if (num_g1_g2 == 0) {
                // no intersection, assign random number without matching (we will time it to 0)
                incre_count = -1;
            } else {
                // calculate incremental benefit
                for (const auto & idx : ngram_exist_lines.at(bigram_key1)) {
                    auto curr_grams = make_unique_bigrams(lines[idx]);
                    if (curr_grams.find(bigram_key2) == curr_grams.end()) {
                        incre_count++;
                    }
                }
            }

            auto abs_inc_benefit = incre_count* num_g1_g2 +g2_filtered* num_g2_not_g1;
            auto reg_w_g1 = regex_grams.at(bigram_key1);
            auto cond_base = (line_count - g1_filtered)*reg_w_g1 + line_count*(regex_count - reg_w_g1);

            if (abs_inc_benefit == 0) {
                // no incremental match just return 0
                continue;
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
            total_time += elapsed_seconds;

            ind_out << get_printable_bigram(bigram_key2) << "\t" << get_printable_bigram(bigram_key1) << "\t";
            ind_out << incre_count << "\t" << num_g1_g2 << "\t";
            ind_out << g2_filtered << "\t" << num_g2_not_g1 << "\t";
            ind_out << abs_inc_benefit << "\t" << cond_base << "\t" << abs_inc_benefit/cond_base << std::endl;

            start = std::chrono::high_resolution_clock::now();
        }
    }

    std::cout << "Incremental time part: " << total_time << std::endl;
    std::cout << "-- Add benefit time to get the total incremental time" << std::endl;
}

// ------------------------------------------------------------------------------------------------------
//  Main
// ------------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    std::string line;
    // read all regexes
    std::vector<std::string> regexes;

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
    auto start = std::chrono::high_resolution_clock::now();
    std::unordered_map<std::pair<char, char>, size_t, hash_pair> gram_to_idx = top_query_ngrams2(regexes);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Frequency time: " << elapsed_seconds << std::endl;

    start = std::chrono::high_resolution_clock::now();
    auto [bigram_regex_num, bigram_regex_idx] = bigram_regex_count(regexes);
    auto [bigram_filtered_line_num, bigram_remain_line_idx] = bigram_line_count(lines, bigram_regex_num);
    auto indi_perc_ranking = individual_benefit_percentage(bigram_filtered_line_num, lines.size(), bigram_regex_num, regexes.size());
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Benefit time: " << elapsed_seconds << std::endl;

    incremental_bigram_line(lines, bigram_filtered_line_num, bigram_remain_line_idx, lines.size(), bigram_regex_num, bigram_regex_idx, regexes.size(), indi_perc_ranking);
}
// g++ -O3 -std=c++17   -Ofast -march=native -mfma -mavx -fomit-frame-pointer -ffp-contract=fast -flto -DARMA_NO_DEBUG -pthread benefit-200.cpp  -L/usr/local/lib/ -lre2 -lstdc++fs -o benefit-200.o