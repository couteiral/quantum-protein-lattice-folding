#ifndef MAIN_H
#define MAIN_H

#include <petscts.h>
#include <list>

#include "timer.h"
#include "options.h"


int main(int argc, char* argv[]);
PetscErrorCode ComputeHamiltonian(TS ts, PetscReal t, Vec X, Mat AA, Mat BB, void *ctx);
PetscErrorCode Normalize(TS ts);

#endif
