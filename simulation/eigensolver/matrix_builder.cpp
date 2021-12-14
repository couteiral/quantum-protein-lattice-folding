/*
    This routine generates hamiltonian 
    matrices with minimum expense in
    memory.
*/

#include <vector>
#include <fstream>
#include <assert.h>
#include <iostream>
#include "petscmat.h" 
#include "options.h"
#include "matrix_builder.h"


MatrixBuilder::MatrixBuilder(Options* opts) : opts(opts) {

    coef    = opts->get_coef();
    size    = opts->get_size();
    n_terms = opts->get_n_terms();
    n_qbits = opts->get_n_qbits();

    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetType(A, MATMPIAIJ);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, size, size);
    MatSetUp(A);

    // First set the diagonal to be non-zero
    for (PetscInt i=0; i<size; i++) {
        PetscScalar zero(0.0, 0.0);
        MatSetValues(A, 1, &i, 1, &i, &zero, INSERT_VALUES);
    }
    MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);

    // If you read when coef = 0, this creates zero entries in
    // the sparse matrix, which is a waste of space
    if (std::abs(coef) > 0.0) {

        Hamiltonian* ham = opts->get_hamiltonian();
        for (int i=0; i<n_terms; i++) {

            PetscScalar weight = ham->get_term_weight(i);
            std::vector<int> Z_array = ham->get_term_qbits(i);
            add_Z_chain(A, coef * weight, Z_array);

        }
    }

    if (std::abs(coef) < 1.0) {

        add_I(A, (1.0-coef)*(0.5*n_qbits), n_qbits-1);

        MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);

        for (int i=0; i<n_qbits; i++) {
            add_X(A, -0.5 * (1.0-coef), i, n_qbits-1);
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);

}

Mat MatrixBuilder::get_matrix() {
    return A;
}

void add_X(Mat A, PetscScalar coef, int pos, int dim) {
    /* 
        Adds a matrix of the form
            coef * I \otimes I \otimes ... \otimes X \otimes ... I \otimes I
        where the position of the X matrix is indicated by
        pos, and the number of matrices in the tensor product
        is indicated by dim.
    */

    assert(dim >= 0);
    assert(pos >= 0);
    assert(pos <= dim);

    PetscInt N = (pos != 0)   ? pos     : 0; // Matrices before X
    PetscInt M = (pos != dim) ? dim-pos : 0; // Matrices after X
    PetscInt N_it = 1ULL << N;
    PetscInt M_it = 1ULL << M;

    /*
        The following section of code employs a simple
        algorithm to build the matrix, just by recognising
        where the non-zero entries are supposed to be.

        We note that:
            - The first N identity matrices create a identity
              matrix of size 2^N
            - The subsequent X matrix creates a 2^(N+1) matrix
              composed of 2^N blocks, each of which is a X matrix
            - The application of the next M identity matrices,
              create 2^N blocks of dimension 2^(M+1), which are
              at the same time composed of two "off-diagonal" blocks
              which are 2^M identity matrices

        Therefore, a matrix of the desired form can be built if
        we just generate 2^N blocks with this form, which can be
        done with negligible cost in memory, with the code below.
    */

    for (PetscInt n=0; n<N_it; n++) {

        PetscInt pt = n * 2 * M_it;

        for (PetscInt m=0; m<M_it; m++) {
            PetscInt pt1 = pt+m;
            PetscInt pt2 = pt+m+M_it; 
            MatSetValues(A, 1, &pt2, 1, &pt1, &coef, INSERT_VALUES);
            MatSetValues(A, 1, &pt1, 1, &pt2, &coef, INSERT_VALUES);
        }

    }

}

void add_Z_chain(Mat A, PetscScalar coef, std::vector<int> Z_chain) {
    /*
        Adds a chain of matrices that includes several Z 
        matrices.
    */

    assert(Z_chain.size() > 0);

    
    /*
        The following section of code employs a simple
        algorithm to build the matrix, just by recognising
        where the non-zero entries are supposed to be.

        A matrix of this form is purely diagonal, formed of
        1s and -1s. We can build it trivially just by
        recognising the symmetries that arise by application
        of the identity and Z matrices.

        This algorithm starts from the end of the chain, and
        progressively builds the matrix. The ith matrix in
        the chain will take the first 2^(N-1-i) numbers and
        replicate them below; if it is a Z, this numbers will
        be preceded by a -1, otherwise it will be a copy.

        Example: I \otimes Z \otimes I

            First:       1  0
                         0  1

            Second:      1  0
                         0  1
                              -1  0
                               0 -1
            Third:       1  0
                         0  1
                              -1  0
                               0 -1
                                     1  0
                                     0  1
                                          -1  0
                                           0 -1
        

    */

    int n = 0;
    PetscInt tot_len;
    PetscScalar* diag;
    tot_len = 1ULL << Z_chain.size();
    PetscMalloc1(tot_len, &diag);

    for (auto it=Z_chain.begin(); it!=Z_chain.end(); it++) {

        if (n == 0) {

            if (*it == 0) { // Start with identity matrix
                diag[0] =  1.0 * coef;
                diag[1] =  1.0 * coef;
            } else if (*it == 1) { // Start with Z matrix
                diag[0] =  1.0 * coef;
                diag[1] = -1.0 * coef;
            } else {
               std::cout << "Fatal error with Z_chain. Value not 0 or 1\n";
               exit(1);
            }

        } else {

            PetscInt pt = 1ULL << n;;

            for (PetscInt i=0; i<pt; i++) { 
                if (*it == 0) {
                    diag[pt+i] =        diag[i];
                } else if (*it == 1) {
                    diag[pt+i] = -1.0 * diag[i];
                } else {
                    std::cout << "Fatal error with Z_chain. Value not 0 or 1\n";
                }
            }

        }

        n++;

    }

    for (PetscInt i=0; i<tot_len; i++) {
        MatSetValues(A, 1, &i, 1, &i, &diag[i], ADD_VALUES);
    }

    PetscFree(diag);

}

void add_I(Mat A, PetscScalar coef, int dim) {
    /* 
        Adds an identity matrix of dimension dim
    */

    assert(dim >= 0);
    PetscInt N_it = 1ULL << (dim+1);

    for (PetscInt n=0; n<N_it; n++) {

        MatSetValues(A, 1, &n, 1, &n, &coef, ADD_VALUES);

    }

}

