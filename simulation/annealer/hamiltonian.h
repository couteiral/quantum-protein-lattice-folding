#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>
#include "petscmat.h"

class Hamiltonian {

    private:

        int n_qbits;
        int n_terms;
        int* qbits_array;
        PetscScalar* weights_array;

    public:

        Hamiltonian(int n_qbits, int n_terms);
        ~Hamiltonian();

        void set_term(int i, PetscScalar weight, std::vector<int> qbits);
        PetscScalar get_term_weight(int i);
        std::vector<int> get_term_qbits(int i);

};

#endif
