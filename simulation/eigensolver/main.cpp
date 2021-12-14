/*

    This program is a utility to find
    the lowest eigenvalues of arbitrary
    Pauli hamiltonians

    Carlos Outeiral Rubiera
    University of Oxford
    May 2019

*/

#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <slepceps.h>

#include "main.h"
#include "timer.h"
#include "options.h"
#include "matrix_builder.h"


static char help[] = "This program finds the lowest eigenvalues of Pauli hamiltonians. Use --help for more information.\n\n";

int main(int argc, char* argv[]) {

    SlepcInitialize(&argc, &argv, (char*)0, help);

    Timer timer;
    timer.tick_zero();

    Options opts(argc, argv);

    Eigensolver eigsolver(&opts, &timer);
    eigsolver.run();

    DataProcessor dataproc(&eigsolver, &opts);
    dataproc.process_output();

    timer.tick_finished();

    if (opts.get_report()) {
        timer.report_time(opts.get_perf_file());
    }

    SlepcFinalize();

}

Eigensolver::Eigensolver(Options* opts, Timer* timer) : opts(opts), timer(timer) {
}

void Eigensolver::run() {

    Mat         A;
    EPS         eps;
    PetscScalar kr, ki;
    PetscInt    nconv;

    MatrixBuilder matbuild(opts);
    A = matbuild.get_matrix();

    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, A, NULL);
    EPSSetProblemType(eps, EPS_HEP);
    EPSSetDimensions(eps, 10 * opts->get_n_qbits(), PETSC_DEFAULT, PETSC_DEFAULT);
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetTolerances(eps, 1.0E-12, PETSC_DEFAULT);
    EPSSetFromOptions(eps);

    timer->tick_start();
    EPSSolve(eps);
    timer->tick_solved();

    EPSGetConverged(eps, &nconv);

    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id == 0) {
        for (int i=0; i<nconv; i++) {

            EPSGetEigenvalue(eps, i, &kr, &ki);
            PetscScalar xr = kr;
            sols.push_back(xr);
        }
    }
    
    EPSDestroy(&eps);

}

std::vector<PetscScalar> Eigensolver::get_sols() {
    return sols;
}

DataProcessor::DataProcessor(Eigensolver* eigsol, Options* opts) 
        : opts(opts) {

    sols = eigsol->get_sols();

}

void DataProcessor::process_output() {

    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id == 0) {

        std::ofstream fpt;
        fpt.open(opts->get_out_file().c_str());

        for (auto it=sols.begin(); it!=sols.end(); it++) {
            fpt << (*it).real() << "\n";
        }

        fpt.close();

    }

}
