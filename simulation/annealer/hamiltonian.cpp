/*

    This class is used to store the
    hamiltonian structure read from
    the file.

    Carlos Outeiral
    University of Oxford
    June 2019

*/

#include <vector>
#include <assert.h>
#include "hamiltonian.h"


Hamiltonian::Hamiltonian(int n_qbits, int n_terms) 
        : n_qbits(n_qbits), n_terms(n_terms) {

    qbits_array = new int[n_terms * n_qbits]();
    weights_array = new PetscScalar[n_terms]();

}

Hamiltonian::~Hamiltonian() {

    delete[] qbits_array;
    delete[] weights_array;

}

void Hamiltonian::set_term(int i, PetscScalar weight, std::vector<int> qbits) {

    assert(i < n_terms);

    weights_array[i] = weight;
    for (int j=0; j<n_qbits; j++) {
        qbits_array[i*n_qbits + j] = qbits[j];
    }

}

PetscScalar Hamiltonian::get_term_weight(int i) {

    assert(i < n_terms);

    return weights_array[i];

}

std::vector<int> Hamiltonian::get_term_qbits(int i) {

    assert(i < n_terms);

    std::vector<int> qbits;

    for (int j=0; j<n_qbits; j++) {
        qbits.push_back( qbits_array[i*n_qbits + j] ); 
    }

    return qbits;
}
