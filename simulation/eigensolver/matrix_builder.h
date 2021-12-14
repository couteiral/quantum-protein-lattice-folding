#ifndef MATRIX_BUILDER_H
#define MATRIX_BUILDER_H

#include <vector>
#include <cassert>
#include "petscmat.h" 
#include "options.h"


class MatrixBuilder {

    private:

        Mat A;
        std::string file;
        int n_qbits;
        int n_terms;
        PetscInt size;
        PetscScalar coef;
        Options* opts;

    public:

        MatrixBuilder(Options* opts);
        Mat get_matrix();

};

void add_X(Mat A, PetscScalar coef, int pos, int dim);
void add_I(Mat A, PetscScalar coef, int dim);
void add_Z_chain(Mat A, PetscScalar coef, std::vector<int> Z_chain);

#endif
