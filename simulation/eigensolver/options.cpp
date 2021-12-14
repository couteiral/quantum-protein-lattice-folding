/*

    This object is used to pass
    the options to the main program.

    Carlos Outeiral
    University of Oxford
    May 2019

*/

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "options.h"
#include "petscmat.h"
#include "hamiltonian.h"


const char program_help[] = "\n"
"    This is a software utility to analyse arbitrary Pauli\n"
"    hamiltonians used for adiabatic quantum computation. It\n"
"    incorporates routines for obtaining the lowest eigenvalues\n"
"    of a hamiltonian (for example to study the spectral gap).\n"
"\n"
"\n"
"      USAGE: \n"
"         ./eigensolver  [flags]\n"
"\n"
"      FLAGS: \n"
"\n"
"         --ham=\n"
"            Selects a file to read the hamiltonian from.\n"
"\n"
"         --coef=\n"
"            Fraction of the annealing process.\n"
"\n"
"         [--out=]\n"
"            Output file. Otherwise it will use default \"output.dat\"\n"
"\n"
"         [--report-time=]\n"
"            Reports times to a file. If not given, timings will not be\n"
"            reported."
"\n"
"         [--help]\n"
"            Prints this help message. Ignores other flags.\n";


bool check_cmd_option(std::string args, std::string key) {

    if (args.find(key) == std::string::npos) {
        return false;
    } else {
        return true;
    }

}

std::string get_cmd_option(std::string args, std::string key) {

    assert(check_cmd_option(args, key));

    std::size_t p_start, p_end, len;
    p_start = args.find(key + '=');
    p_end = args.find("--", p_start+1);
    len = p_end - (p_start + key.length()+1);

    std::string result= args.substr(p_start + key.length()+1, len);
    remove_trailing_whitespaces(result);

    return result;

}

std::string process_args(int argc, char* argv[]) {

    std::string args(argv[1]);

    for (int i=2; i<argc; i++) {
        args.append(argv[i]).append(" ");
    }

    return args;

}

void remove_trailing_whitespaces(std::string &str) {

    size_t ws_start=0, ws_end=0;

    for (int i=0; i<str.length(); i++) {
        if (str[i] == ' ') {
            ws_start++;
        } else {
            break;
        }
    }

    for (int i=str.length()-1; i>=0; i++) {
        if (str[i] == ' ') {
            ws_end++;
        } else {
            break;
        }
    }

    if (ws_start == ws_end and ws_end != 0) {
        str = std::string("\0");
    } else {
        str.erase(str.length()-ws_end, ws_end);
        str.erase(0, ws_start);
    }

}

Options::Options(int argc, char* argv[]) {

    /*
        First, we process the command-line arguments
    */

    if (argc == 1) {
        std::cout << "ERROR: Not enough arguments provided. Printing help...\n\n";
        std::cout << program_help;
        exit(1);
    }

    std::string args = process_args(argc, argv);

    if (check_cmd_option(args, "--help")) {
        std::cout << program_help;
        exit(1);
    }


    if (!check_cmd_option(args, "--coef")) {
        std::cout << "ERROR: Coefficient not specified, use flag --coef.\n";
        exit(1);
    }
    coef = std::stod( get_cmd_option(args, "--coef").c_str() );
    if (std::abs(coef) < 0.0 or std::abs(coef) > 1.0) {
        std::cout << "ERROR: Coefficient must be in range [0, 1]\n";
        exit(1);
    }

    if (!check_cmd_option(args, "--ham")) {
        std::cout << "ERROR: Hamiltonian file not specified, use flag --ham.\n";
        exit(1);
    }

    ham_file = get_cmd_option(args, "--ham");

    if (check_cmd_option(args, "--out")) {
        out_file = get_cmd_option(args, "--out");
    } else {
        out_file = "output.dat";
    }

    if (check_cmd_option(args, "--report-time")) {
        report = true;
        perf_file = get_cmd_option(args, "--report-time");
    } else {
        report = false;
    }

    /*
        Then, we read in data from files
    */

    std::ifstream fpt(ham_file);
    if (!fpt) {
        std::cout << "ERROR: Hamiltonian file" << ham_file << "does not exist.\n";
        exit(1);
    }

    std::string line;
    getline(fpt, line);
    std::stringstream strst(line);

    if (!(strst >> n_qbits)) {
        std::cout << "ERROR: Problem reading header in " << ham_file << "\n";
        fpt.close();
        exit(1);
    }
    size = 1ULL << n_qbits;
    if (!(strst >> n_terms)) {
        std::cout << "ERROR: Problem reading header in " << ham_file << "\n";
        fpt.close();
        exit(1);
    }
    
    ham = new Hamiltonian(n_qbits, n_terms);

    for (int i=0; i<n_terms; i++) {
       
        getline(fpt, line);
        std::stringstream strst(line);

        PetscScalar weight;
        strst >> weight;
       
        int qbit;
        std::vector<int> qbit_list;

        while ( strst >> qbit ) {
            if (qbit == 0) {
                qbit_list.push_back(0);
            } else {
                qbit_list.push_back(1);
            }
        }
        if (qbit_list.size() < n_qbits) {
            std::cout << "ERROR: Could not read hamiltonian file.\n";
            fpt.close();
            exit(1);
        }

        ham->set_term(i, weight, qbit_list);

    }

    fpt.close();

}

Options::~Options() {
    delete ham;
}

Hamiltonian* Options::get_hamiltonian() {
    return ham;
}

std::string Options::get_ham_file() {
    return ham_file;
}

std::string Options::get_out_file() {
    return out_file;
}

std::string Options::get_perf_file() {
    return perf_file;
}

bool Options::get_report() {
    return report;
}

PetscScalar Options::get_coef() {
    return coef;
}

int Options::get_n_terms() {
    return n_terms;
}

int Options::get_n_qbits() {
    return n_qbits;
}

unsigned long long Options::get_size() {
    return size;
}
