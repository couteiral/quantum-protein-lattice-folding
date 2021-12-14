#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <list>
#include "timer.h"
#include "options.h"
#include "matrix_builder.h"


class Eigensolver {

    private:

        Timer*   timer;
        Options* opts;

        std::vector<PetscScalar> sols;

    public:

        Eigensolver(Options* opts, Timer* timer);
        void run();
        std::vector<PetscScalar> get_sols();

};

class DataProcessor {

    private:

        Options* opts;
        std::vector<PetscScalar> sols;

    public:

        DataProcessor(Eigensolver* eigsol, Options* opts);
        void process_output();

};

#endif
