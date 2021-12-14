/*
    This program generates the list of
    conformations, ordered by descendent
    values of energy.
*/


#include <set>
#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>

#include "solution.h"
#include "energy_functions.h"


std::set<char> hp_model = {'H', 'P'};
std::set<char> mj_model = {'Y', 'C', 'K', 'E', 'L', 'V', 'R', 'H',
                           'P', 'M', 'A', 'W', 'S', 'I', 'T', 'Q', 
                           'F', 'D', 'G', 'N'};


bool check_if_hp(std::string ps_seq) {

    std::set<char> alphabet, diff;

    for (int i=0; i<ps_seq.size(); i++) {
        alphabet.insert(ps_seq[i]);
    }

    if (std::includes(hp_model.begin(), hp_model.end(),
                      alphabet.begin(), alphabet.end())) {
        return true;
    } else {
        return false;
    }

}

bool check_if_mj(std::string ps_seq) {

    std::set<char> alphabet, diff;

    for (int i=0; i<ps_seq.size(); i++) {
        alphabet.insert(ps_seq[i]);
    }

    if (std::includes(mj_model.begin(), mj_model.end(),
                      alphabet.begin(), alphabet.end())) {
        return true;
    } else {
        return false;
    }

}

bool check_only_hp(std::string ps_seq) {

    bool could_be_hp=true;
    std::set<char> alphabet, diff;

    for (int i=0; i<ps_seq.size(); i++) {
        alphabet.insert(ps_seq[i]);
    }

    for (auto it=mj_model.begin(); it!=mj_model.end(); it++) {
        if (*it == 'H' or *it == 'P') continue;
        if (alphabet.find(*it) != alphabet.end()) {
            could_be_hp = false;
            break;
        }
    }

    return could_be_hp;

}

std::list<std::pair<int, int>> make_contact_list(std::string line) {

    int number;
    std::list<std::pair<int, int>> contact_list;
    std::vector<int> number_list;
    std::stringstream stream(line);

    while (stream >> number) {
        number_list.push_back(number);
    }

    int n_targets = number_list.size(); // always even
    for (int i=0; i<n_targets/2; i++) {
        std::pair<int, int> pair;
        pair.first = number_list[i];
        pair.second = number_list[i + n_targets/2];
        contact_list.push_back(pair);
    }

    return contact_list;

}

int main(int argc, char* argv[]) {

    // Process input and check

    if (argc < 2) {
        std::cout << "Usage: \n";
        std::cout << "  energy-calc   file_prefix  sequence\n";
        exit(1);
    }

    std::string file_prefix(argv[1]);
    std::string ps_seq(argv[2]);
    std::string default_enc;

    if (argc >= 4) {
        std::string temp(argv[3]);
        default_enc = temp;
    }

    bool is_hp = check_if_hp(ps_seq);
    bool is_mj = check_if_mj(ps_seq);

    if (is_hp and is_mj) {

        bool only_hp = check_only_hp(ps_seq);
    
        if (!only_hp) {
            is_hp = false;
        } else {
            std::cout << "Warning: introduced chain is ambigous, so HP model was selected.\n";
            is_mj = false;
        }

    }

    if (!(is_hp or is_mj)) {

        std::cout << "Fatal error: unrecognised protein model.\n";
        exit(1);

    }

    if (argc >= 4) {
        if (default_enc == "HP") {
            is_hp = true;
            is_mj = false;
        } else if (default_enc == "MJ") {
            is_hp = false;
            is_mj = true;
        }
    }

    // Get metadata

    std::ifstream meta_file;
    meta_file.open(file_prefix + ".meta");
    if (!meta_file) {
        std::cout << "Fatal error: Could not read " << file_prefix << ".meta\n";
        exit(1);
    }

    unsigned long long n_confs;
    meta_file >> n_confs;
    meta_file.close();

    // Test conformations

    std::ifstream contact_file;
    contact_file.open(file_prefix + ".contacts");
    if (!contact_file) {
        std::cout << "Fatal error: Could not read " << file_prefix << ".contacts\n";
        exit(1);
    }
    double* energies = new double[n_confs];

    for (unsigned long long i=0; i<n_confs; i++) {

        std::list<std::pair<int, int>> contact_list;
        std::string line;
        std::getline(contact_file, line);

        contact_list = make_contact_list(line);

        for (auto it = contact_list.begin(); it != contact_list.end(); it++){

            char aa1 = ps_seq[it->first];
            char aa2 = ps_seq[it->second];            

            if (is_hp) {
                energies[i] += hp_energy(aa1, aa2);
            } else if(is_mj) {
                energies[i] += mj_energy(aa1, aa2);
            } else {
                std::cout << "Something went wrong... neither MJ nor HP!\n";
            }
        }
        
    }

    contact_file.close();

    // Read in turn sequences

    std::ifstream seq_file;
    seq_file.open(file_prefix + ".seq");
    if (!seq_file) {
        std::cout << "Fatal error: Could not read " << file_prefix << ".seq\n";
        exit(1);
    }

    std::list<Solution> solutions;
    for (unsigned long long i=0; i<n_confs; i++) {

        std::string seq;
        seq_file >> seq;
        std::cout << seq << "\n";

        Solution sol;
        sol.seq = seq;
        sol.energy = energies[i];

        solutions.push_back(sol);

    }

    seq_file.close();

    // Write everything to file

    solutions.sort(cmp_solution);

    std::ofstream out_file;
    out_file.open(ps_seq + "_" + file_prefix + ".sol");

    for (auto it = solutions.begin(); it != solutions.end(); it++) {
        out_file << it->seq << "  " << it->energy << "\n";
    }

    out_file.close();

}
