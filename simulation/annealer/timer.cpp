#include <mpi.h>
#include <string>
#include <fstream>
#include <iostream>
#include "timer.h"


void Timer::tick_zero() {
    t_began = MPI_Wtime();
}

void Timer::tick_start() {
    t_start = MPI_Wtime();
}

void Timer::tick_solved() {
    t_solved = MPI_Wtime();
}

void Timer::tick_finished() {
    t_finished = MPI_Wtime();
}

double Timer::setup_time() {
    return t_start - t_began;
}

double Timer::solution_time() {
    return t_solved - t_start;
}

double Timer::setdown_time() {
    return t_finished - t_solved;
}

double Timer::total_time() {
    return t_finished - t_began;
}

void Timer::report_time(std::string perf_file) {

    std::ofstream fpt;
    fpt.open(perf_file);

    fpt << "SETUP    TIME = " << this->setup_time() << "\n";
    fpt << "SOLUTION TIME = " << this->solution_time() << "\n";
    fpt << "SETDOWN  TIME = " << this->setdown_time() << "\n";
    fpt << "TOTAL    TIME = " << this->total_time() << "\n";

    fpt.close();

}
