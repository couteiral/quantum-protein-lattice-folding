#ifndef TIMER_H
#define TIMER_H

class Timer {

    private:

        double t_began;
        double t_start;
        double t_solved;
        double t_finished;

    public:

        void tick_zero();
        void tick_start();
        void tick_solved();
        void tick_finished();
        void report_time(std::string perf_file);

        double setup_time();
        double solution_time();
        double setdown_time();
        double total_time();

};

#endif
