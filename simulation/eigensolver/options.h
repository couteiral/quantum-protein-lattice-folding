#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include "petscmat.h"
#include "hamiltonian.h"

class Options {

    private:

        Hamiltonian* ham;

        std::string ham_file;
        std::string out_file;
        std::string perf_file;

        int n_qbits;
        int n_terms;
        bool report;
        PetscScalar coef;
        unsigned long long size;

    public:

        Options(int argc, char* argv[]);
        ~Options();

        Hamiltonian* get_hamiltonian();
        std::string get_ham_file();
        std::string get_out_file();
        std::string get_perf_file();
        bool get_report();
        PetscScalar get_coef();
        int get_n_terms();
        int get_n_qbits();
        unsigned long long get_size();

};

bool check_cmd_option(std::string args, std::string key);
std::string get_cmd_option(std::string args, std::string key);
void remove_trailing_whitespaces(std::string &str);

#endif
