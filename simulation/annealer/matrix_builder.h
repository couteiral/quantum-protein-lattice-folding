#ifndef MATRIX_BUILDER_H
#define MATRIX_BUILDER_H

#include <vector>
#include <cassert>
#include "petscmat.h" 
#include "options.h"


class MatrixBuilder {

    public:

        Mat A;
        Mat H_start, H_end;
        std::string file;
        int n_qbits;
        int n_terms;
        PetscInt size;
        Options* opts;

        MatrixBuilder(Options* opts, Mat A);
        ~MatrixBuilder();

};

void add_Z(Mat A, PetscScalar coef, int pos, int dim);
void add_X(Mat A, PetscScalar coef, int pos, int dim);
void add_I(Mat A, PetscScalar coef, int dim);
void add_Z_chain(Mat A, PetscScalar coef, std::vector<int> Z_chain);
PetscErrorCode MatPreallocate(Mat A, PetscInt dim);

int read_nqbits(std::string file);

#endif
