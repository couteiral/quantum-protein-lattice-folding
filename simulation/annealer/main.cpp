/*

    This program is a utility to simulate
    quantum annealing with Pauli hamiltonians.

    Carlos Outeiral Rubiera
    University of Oxford
    May 2019

*/

#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <petscts.h>

#include "main.h"
#include "timer.h"
#include "options.h"
#include "matrix_builder.h"


static char help[] = "This program computes the adiabatic evolution of Pauli hamiltonians. Use --help for more information.\n\n";

int main(int argc, char* argv[]) {

    /*

        This function contains the main body of the program:

            - Reads the input parameters and files (through the
              Options object, which is then used as a bundle of
              variables that is passed to functions)

            - Sets up the PETSc TS environment.

            - Sets up the wavefunction vector.

            - Sets up the initial and final hamiltonian matrices,
              which will be combined during the integration to
              create the instantaneous hamiltonian.

            - Reports the results to an output file.

    */

    Options opts(argc, argv);
    PetscInitialize(&argc, &argv, (char*)0, help);

    Timer timer;
    timer.tick_zero();

    // Basic program variables
    TS ts;
    Mat H;
    Vec psi;
    PetscReal tstep = 0.0001;

    // Variables provided in Options object
    PetscInt n_qbits = opts.get_n_qbits();
    PetscInt size = opts.get_size();
    PetscScalar time = opts.get_time();
    PetscInt max_steps = (PetscInt) std::abs(time) / tstep;
    std::string out_file = opts.get_out_file();

    // Set up PETSc TS
    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetType(ts, TSRK);
    //TSRKSetType(ts, TSRK4);
    TSRKSetType(ts, TSRK5F);
    TSSetProblemType(ts, TS_LINEAR);
    TSSetTimeStep(ts, tstep);
    TSSetMaxSteps(ts, max_steps);
    TSSetMaxTime(ts, std::abs(opts.get_time()));
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetTolerances(ts, 1.0e-10, NULL, PETSC_DECIDE, NULL);
    TSSetFromOptions(ts);

    // Uncomment to disable adaptivity
    //TSAdapt adapt;
    //TSGetAdapt(ts, &adapt);
    //TSAdaptSetType(adapt, TSADAPTNONE);

    // Set up wavefunction vector
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size, &psi);
    PetscScalar  val = 1/sqrt(size);
    VecSet(psi, val);
    VecAssemblyBegin(psi);
    VecAssemblyEnd(psi);

    // Allocate hamiltonian matrix
    MatCreate(PETSC_COMM_WORLD, &H);
    MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, size, size);
    MatSetType(H, MATMPIAIJ);
    MatSetFromOptions(H);
    MatSetUp(H);
    MatPreallocate(H, n_qbits); // this is a custom function

    // Compute H_start and H_end and pass it to hamiltonian
    // computation function
    MatrixBuilder* matbuild = new MatrixBuilder(&opts, H);
    TSSetRHSFunction(ts, NULL, TSComputeRHSFunctionLinear, matbuild);
    TSSetRHSJacobian(ts, H, H, ComputeHamiltonian, matbuild);
    TSSetPreStep(ts, Normalize);

    // Start integrator
    timer.tick_start();
    TSSolve(ts, psi);
    timer.tick_solved();

    // Report data to file
    Vec conj;
    VecDuplicate(psi, &conj);
    VecCopy(psi, conj);
    VecConjugate(conj);
    VecPointwiseMult(conj, psi, conj);

    PetscViewer view;
    PetscViewerCreate(PETSC_COMM_WORLD, &view);
    PetscViewerSetType(view, PETSCVIEWERASCII);
    PetscViewerFileSetMode(view, FILE_MODE_WRITE);
    PetscViewerFileSetName(view, out_file.c_str());
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, out_file.c_str(), &view);
    VecView(conj, view);
    PetscViewerDestroy(&view);

    // Destroy everything
    TSDestroy(&ts);
    MatDestroy(&H);
    VecDestroy(&psi);
    delete matbuild;
    PetscFinalize();
    timer.tick_finished();

    // Report timings
    if (opts.get_report()) {
        timer.report_time(opts.get_perf_file());
    }

}

PetscErrorCode ComputeHamiltonian(TS ts, PetscReal t, Vec X, Mat AA, Mat BB, void *ctx) {

    /*

        This function computes the instantaneous hamiltonian. The
        MatrixBuilder object contains the matrices H_start and H_end,
        and this routine combines them into another matrix to compute
        the evolution.

    */

    MatrixBuilder* matbuild = (MatrixBuilder*) ctx;

    PetscScalar coef = t / matbuild->opts->get_time();
    if (std::abs(coef) < 0.0 or std::abs(coef) > 1.0) {
        std::cout << "RUNTIME ERROR: coefficient " << coef << " out of range [0, 1]\n";
        exit(1);
    }

    MatZeroEntries(AA);
    MatAXPY(AA, PETSC_i * (1.0-coef), matbuild->H_start, SAME_NONZERO_PATTERN);
    MatAXPY(AA, PETSC_i * coef, matbuild->H_end, SUBSET_NONZERO_PATTERN);

    return 0;

}

PetscErrorCode Normalize(TS ts) {

    /*

        Renormalises the wavefunction vector every 1,000 steps,
        for numerical reasons. 

    */

    PetscInt step;
    TSGetStepNumber(ts, &step);

    if (step % 10 == 0) {

        Vec vec;
        PetscReal norm;

        TSGetSolution(ts, &vec);
        VecNorm(vec, NORM_2, &norm);
        VecScale(vec, (PetscScalar) 1.0/norm);

    }

    return 0;

}
